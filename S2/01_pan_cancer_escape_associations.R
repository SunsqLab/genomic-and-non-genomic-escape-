# S2: Co-occurrence of immune escape features
# In RStudio, set the working directory to this file's directory, edit the paths,
# then Source this file or run it from top to bottom. The functions are sourced below.
# From another directory: source("path/to/01_pan_cancer_escape_associations.R", chdir = TRUE).
# Inputs keep the ORIGINAL RDS objects. No new TSV schema or preprocessing is required.
# DataCenter.rds must contain ID and Path; Path points to each original cohort folder.
# Edit the two query paths to match your existing folders; do not move the data.
DataCenter.file <- "/path/to/DataCenter.rds"
PanCancer.path <- "/path/to/PanCancer_TCGA.dataset"
escape.query <- "Results/BioImmune/ImmuneEscape"
b2m.query <- "Results/BioGenomics/B2M_biallelic_inactivation_matrix.rds"
functions.dir <- normalizePath("functions", mustWork = TRUE)

library(magrittr)
library(ggplot2)
library(corrplot)
library(dplyr)

Dir.output="results"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive = TRUE) }
Dir.output <- normalizePath(Dir.output, mustWork = TRUE)

profile.escape <- readRDS(file.path(PanCancer.path, escape.query, "ProfilingImmuneEvading.rds"))

source(file.path(functions.dir, "GetInfor.PatientCenter.R"))

patient.center=readRDS(file.path(PanCancer.path, "PatientCenter", "PatientCenter.rds"))
# Primary samples only; use raw continuous scores for non-genomic features.
profile.escape = lapply(profile.escape, function(x){
            x[GetInfor.PatientCenter(patient.center, SampleID=x$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
        })

profile.escape.df=Reduce(function(x, y) dplyr::full_join(x, y, by = "SampleID"), profile.escape)
profile.escape.df$all.APM.mut=ifelse(profile.escape.df$all.APM.mut>0, TRUE, FALSE)
profile.escape.df$HLA.mut=ifelse(profile.escape.df$HLA.mut>0, TRUE, FALSE)
profile.escape.df=Filter(function(x) !all(is.na(x)), profile.escape.df)
profile.escape.df <- list(PanCancer = profile.escape.df)

source(file.path(functions.dir, "DatasetLabels.R"))
source(file.path(functions.dir, "CombineData.XTeam.R"))
ID.xteam <- "PanCancer_TCGA.dataset"

ID.xteams <- get(ID.xteam)
B2M.loh <- CombineData.XTeam(ID.xteams, file.query = b2m.query)

B2M.loh.df <- dplyr::bind_rows(lapply(B2M.loh, function(x) {
    data.frame(
        SampleID = colnames(x),
        B2M.LOH = x["B2M", ] == 1,
        stringsAsFactors = FALSE
    )
}))
profile.escape.df <- lapply(profile.escape.df, \(x) { m <- match(x$SampleID, B2M.loh.df$SampleID); x$B2M.LOH[!is.na(m)] <- B2M.loh.df$B2M.LOH[m[!is.na(m)]]; x })

features = list(
    APMalt = c("HLA.LOH", "B2M.LOH", "HLA.mut", "B2M.mut", "all.APM.mut"),
    ActMalt = c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
    Checkalt = c("CD274.deepAmp"),
    Checkexp = c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    Supprcell = c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    Supprsig = c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    Actdown = c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown = c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo = c("mean.neo.exprs")
)
all.features = unlist(features)
names(all.features) <- rep(c("#c82621", "#fa8b69ff", "#F6C141", "#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"), lengths(features))
genome.features <- unlist(features[c("APMalt", "ActMalt", "Checkalt")],use.names = FALSE)
non.genomic <- unlist(features[c("Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo")],use.names = FALSE)

# Genomic pairs: preserve the source direction-selected one-sided Fisher tests.
# HLA.mut and B2M.mut are nested in all.APM.mut and excluded from these pairs.
source(file.path(functions.dir, "EscapePairFisher.R"))
temp.genome.features <- genome.features[!genome.features %in% c("HLA.mut", "B2M.mut")]
stat.data=EscapePairFisher(data=profile.escape.df, escape.feature=temp.genome.features)
data=stat.data      %>%
    dplyr::mutate(log10.OR=log10(odds_ratio),
        log10.FDR=log10(p_adj)*(-1),
        color=factor(ifelse(odds_ratio>1, "OR > 1", "OR < 1"), levels=c("OR > 1", "OR < 1")),
        significance=factor(ifelse(log10.FDR>2, "strong", ifelse(log10.FDR>1.301, "medium", "weak")), levels=c("strong", "medium", "weak")),
        label=paste(geneA, geneB, sep=':'))
data$label[data$significance=='weak']=''

p=ggplot(data, aes(x = log10.OR, y = log10.FDR, fill = color)) +
        geom_point(aes(size = significance), shape = 21) +
        scale_fill_manual(values = c("OR > 1"="#E41A1C", "OR < 1"="#377eb8")) +
        scale_size_manual(values=c("strong"=4, "medium"=2, "weak"=0.1), labels = c("strong" = "FDR < 0.01", "medium" = "0.01 < FDR < 0.05", "weak" = "FDR > 0.05"))+
        ggrepel::geom_text_repel(aes(label=label), size=2.5, max.overlaps=50, segment.color = "#d6ccc2", segment.size=0.3)+
        theme_bw() +
        xlim(c(-1, 3.5))+
        labs(y='-log10 (p.value)', x='log10 (Odds Ratio)')
ggsave(file.path(Dir.output, '1.genomic_escape_associations.pan_cancer.pdf'), p, width=7, height=5)

source(file.path(functions.dir, "EscapePairSpearman.R"))
non.genomic <- unlist(features[c("Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo")],use.names = FALSE)
names(non.genomic) <- names(all.features)[match(non.genomic, all.features)]
# Panel-wide complete cases; BH includes both matrix directions, as in the source.
stat.data=EscapePairSpearman(data=profile.escape.df, non.genomic)
data=stat.data      %>%
    dplyr::mutate(significance=factor(ifelse(p_adj<0.01, "strong", ifelse(p_adj<0.05, "medium", "weak")), levels=c("strong", "medium", "weak")),
        feaA=factor(feaA, levels=unique(feaA)[order(match(unique(feaA), non.genomic))]),
        feaB=factor(feaB, levels=unique(feaB)[order(match(unique(feaB), non.genomic))]))

data.cor <- stat.data %>%
  dplyr::select(feaA, feaB, cor) %>%
  tidyr::pivot_wider(names_from = feaB, values_from = cor) %>%
  tibble::column_to_rownames("feaA") %>%
  as.matrix()

data.p <- stat.data %>%
  dplyr::select(feaA, feaB, p_adj) %>%
  tidyr::pivot_wider(names_from = feaB, values_from = p_adj) %>%
  tibble::column_to_rownames("feaA") %>%
  as.matrix()

var.order <- intersect(colnames(data.cor), rownames(data.cor))
data.cor  <- data.cor[var.order, var.order]
data.p    <- data.p[var.order, var.order]

tl.col <- names(non.genomic)[match(var.order, non.genomic)]
col.RWB <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(200)

pdf(file.path(Dir.output, '2.non_genomic_escape_associations.pan_cancer.pdf'), width = 16, height = 12)
p = corrplot(data.cor,method = "square",col = col.RWB, outline = "grey",
  order = "original",diag = TRUE,tl.cex = 1.2,tl.col = tl.col, tl.srt = 45,
  addgrid.col = "grey",type = "lower",
  p.mat = data.p,
  sig.level = 0.05,
  insig = "blank", na.label = " "
)
dev.off()

source(file.path(functions.dir, "EscapePairWilcoxon.R"))

# Pairwise complete cases; at least 5 samples per group; BH within this cohort.
stat.data <- EscapePairWilcoxon(data = profile.escape.df, disc.feature = genome.features, cont.feature = non.genomic)
stat.data <- na.omit(stat.data)
data <- stat.data %>%
    dplyr::mutate(
        log2.fold_change = log2(fold_change),
        log10.p_adj = log10(p_adj) * (-1),
        color = factor(ifelse(fold_change > 1, "fold change > 1", "fold change < 1"), levels = c("fold change > 1", "fold change < 1")),
        significance = factor(ifelse(p_adj < 0.01, "strong", ifelse(p_adj < 0.05, "medium", "weak")), levels = c("strong", "medium", "weak"))
    )
data$label[data$p_adj > 0.05] <- ""

data <- data %>%
    dplyr::mutate(
        sig_label = ifelse(significance == "strong", "**", ifelse(significance == "medium", "*", "")),
        feaA = factor(feaA, levels = genome.features),
        feaB=factor(feaB, levels=rev(non.genomic)),
        feaA.class = factor(rep(names(features), lengths(features))[match(feaA, unlist(features))],levels=c("APMalt","ActMalt","Checkalt")),
        feaB.class = factor(rep(names(features), lengths(features))[match(feaB, unlist(features))],levels=c("Checkexp","Supprcell","Supprsig","Actdown","APMdown","Neo")),
    )

feature_col_map <- setNames(names(all.features), all.features)
p <- ggplot(data, aes(x = feaA, y = feaB, fill = log2.fold_change)) +
    geom_tile(color = "grey85", size = 0.3) +
    geom_text(aes(label = sig_label), size = 3, color = "black") +
    scale_fill_gradient2(
        low = "#4B71B3", mid = "white", high = "#D66B6B",
        midpoint = 0, limits = c(-3, 2), breaks = -3:2,
        oob = scales::squish, na.value = "white", name = "log2FC",
        guide = guide_colorbar(
            direction = "horizontal",
            title.position = "bottom",
            title.hjust = 0.5,
            label.position = "bottom",
            barwidth = grid::unit(4.8, "cm"),
            barheight = grid::unit(0.35, "cm"),
            ticks = TRUE
        )
    ) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
        axis.text.y = ggtext::element_markdown(size = 10),
        axis.title = element_blank(), panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        strip.background = element_blank(), strip.text = element_blank(),
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels = function(x) {paste0("<span style='color:", feature_col_map[x], "'>", x, "</span>")}) +
    facet_grid(cols = vars(feaA.class), scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL)
ggsave(file.path(Dir.output, "3.genomic_non_genomic_associations.pan_cancer.pdf"), p, width = 5, height = 9)

# Original Sections 4-6: feature consensus clustering and module summaries.
# These analyses were already present in the original main script.
non.genomic.consensus.k <- "auto"
non.genomic.consensus.k.range <- 2:8
non.genomic.consensus.pac.tolerance <- 0.02
non.genomic.consensus.reps <- 500
non.genomic.consensus.p.item <- 0.8
non.genomic.consensus.p.feature <- 0.8
non.genomic.consensus.seed <- 20260605
non.genomic.consensus.min.finite <- 30
non.genomic.consensus.file.prefix <- "4.non_genomic_consensus_clustering.pan_cancer"
non.genomic.consensus.output.dir <- file.path(Dir.output, non.genomic.consensus.file.prefix)
if (!dir.exists(non.genomic.consensus.output.dir)) {
    dir.create(non.genomic.consensus.output.dir, recursive = TRUE)
}

feature.color <- setNames(names(all.features), all.features)
feature.category <- rep(names(features), lengths(features))
names(feature.category) <- unlist(features, use.names = FALSE)
category.colors <- vapply(names(features), function(category) {
    current.color <- unique(unname(feature.color[features[[category]]]))
    current.color <- current.color[!is.na(current.color)]
    if (length(current.color) == 0) "#BDBDBD" else current.color[1]
}, character(1))
names(category.colors) <- names(features)

manual.k <- NA_integer_
if (!identical(non.genomic.consensus.k, "auto")) {
    manual.k <- suppressWarnings(as.integer(non.genomic.consensus.k))
    if (is.na(manual.k) || length(manual.k) != 1 || manual.k < 2) {
        stop("non.genomic.consensus.k must be 'auto' or an integer >= 2.")
    }
}

non.genomic.available <- non.genomic[non.genomic %in% colnames(profile.escape.df$PanCancer)]
non.genomic.missing <- setdiff(non.genomic, non.genomic.available)
if (length(non.genomic.missing) > 0) {
    warning("Skipping non-genomic features absent from profile.escape.df$PanCancer: ",
            paste(non.genomic.missing, collapse = ", "))
}

non.genomic.consensus.df <- profile.escape.df$PanCancer[, c("SampleID", non.genomic.available), drop = FALSE]
non.genomic.consensus.matrix <- data.frame(
    lapply(non.genomic.consensus.df[, non.genomic.available, drop = FALSE], function(x) {
        suppressWarnings(as.numeric(as.character(x)))
    }),
    check.names = FALSE
)
rownames(non.genomic.consensus.matrix) <- make.unique(as.character(non.genomic.consensus.df$SampleID))

non.genomic.consensus.normalization <- "raw_continuous_values_no_cancer_zscore"
message("Section 4 uses raw continuous non-genomic feature values without cancer-type normalization.")

feature.qc <- data.frame(
    feature = colnames(non.genomic.consensus.matrix),
    original_category = unname(feature.category[colnames(non.genomic.consensus.matrix)]),
    n_finite = vapply(non.genomic.consensus.matrix, function(x) sum(is.finite(x)), integer(1)),
    n_missing = vapply(non.genomic.consensus.matrix, function(x) sum(!is.finite(x)), integer(1)),
    missing_rate = vapply(non.genomic.consensus.matrix, function(x) mean(!is.finite(x)), numeric(1)),
    sd = vapply(non.genomic.consensus.matrix, function(x) stats::sd(x, na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
)
feature.qc$normalization <- non.genomic.consensus.normalization
feature.qc$keep_for_consensus <- feature.qc$n_finite >= non.genomic.consensus.min.finite &
    is.finite(feature.qc$sd) & feature.qc$sd > 0
utils::write.table(
    feature.qc,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".feature_qc.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

valid.feature <- setNames(feature.qc$keep_for_consensus, feature.qc$feature)
if (!all(valid.feature)) {
    warning("Skipping non-genomic features with too few finite values or zero variance: ",
            paste(names(valid.feature)[!valid.feature], collapse = ", "))
}
non.genomic.consensus.matrix <- non.genomic.consensus.matrix[, valid.feature, drop = FALSE]
if ("mean.neo.exprs" %in% colnames(non.genomic.consensus.matrix)) {
    message("mean.neo.exprs is retained for consensus clustering.")
} else {
    warning("mean.neo.exprs is excluded from consensus clustering; check finite-value counts and variance in feature_qc.tsv.")
}

if (ncol(non.genomic.consensus.matrix) < 3) {
    stop("Fewer than 3 valid non-genomic features remain for consensus clustering.")
}

missing.filter.summary <- data.frame(
    n_samples_in_matrix = nrow(non.genomic.consensus.matrix),
    n_features_after_feature_filter = ncol(non.genomic.consensus.matrix),
    n_complete_case_samples_all_features = sum(apply(non.genomic.consensus.matrix, 1, function(x) all(is.finite(x)))),
    normalization = non.genomic.consensus.normalization,
    distance_missing_strategy = "pairwise.complete.obs",
    stringsAsFactors = FALSE
)
utils::write.table(
    missing.filter.summary,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".missing_filter_summary.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

ConsensusSpearmanPairwiseDistance <- function(x) {
    cor.mat <- suppressWarnings(stats::cor(t(x), method = "spearman", use = "pairwise.complete.obs"))
    cor.mat[!is.finite(cor.mat)] <- 0
    cor.mat <- (cor.mat + t(cor.mat)) / 2
    diag(cor.mat) <- 1
    cor.mat[cor.mat > 1] <- 1
    cor.mat[cor.mat < -1] <- -1
    stats::as.dist(1 - cor.mat)
}
assign("ConsensusSpearmanPairwiseDistance", ConsensusSpearmanPairwiseDistance, envir = .GlobalEnv)

max.allowed.k <- ncol(non.genomic.consensus.matrix) - 1
non.genomic.consensus.k.range <- sort(unique(non.genomic.consensus.k.range))
non.genomic.consensus.k.range <- non.genomic.consensus.k.range[
    non.genomic.consensus.k.range >= 2 & non.genomic.consensus.k.range <= max.allowed.k
]
if (length(non.genomic.consensus.k.range) == 0) {
    stop("No k in non.genomic.consensus.k.range is compatible with the number of valid features.")
}
if (!is.na(manual.k) && manual.k > max.allowed.k) {
    stop("The manually specified k exceeds the number of valid non-genomic features minus 1.")
}

non.genomic.consensus.k.eval <- non.genomic.consensus.k.range
if (!is.na(manual.k)) {
    non.genomic.consensus.k.eval <- sort(unique(c(non.genomic.consensus.k.eval, manual.k)))
}
max.k <- max(non.genomic.consensus.k.eval)

non.genomic.consensus.result <- local({
    consensus.tmp.plot.name <- paste0(".ConsensusClusterPlus_tmp_", Sys.getpid())
    consensus.tmp.plot.dir <- file.path(non.genomic.consensus.output.dir, consensus.tmp.plot.name)
    old.wd <- setwd(non.genomic.consensus.output.dir)
    on.exit({
        setwd(old.wd)
        unlink(consensus.tmp.plot.dir, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    ConsensusClusterPlus::ConsensusClusterPlus(
        d = as.matrix(non.genomic.consensus.matrix),
        maxK = max.k,
        reps = non.genomic.consensus.reps,
        pItem = non.genomic.consensus.p.item,
        pFeature = non.genomic.consensus.p.feature,
        clusterAlg = "hc",
        distance = "ConsensusSpearmanPairwiseDistance",
        innerLinkage = "average",
        finalLinkage = "average",
        seed = non.genomic.consensus.seed,
        title = consensus.tmp.plot.name,
        plot = "pdf",
        writeTable = FALSE
    )
})

calc_pac <- function(consensus.matrix, lower = 0.1, upper = 0.9) {
    x <- consensus.matrix[upper.tri(consensus.matrix)]
    mean(x > lower & x < upper, na.rm = TRUE)
}

summarise_consensus_k <- function(k) {
    res <- non.genomic.consensus.result[[k]]
    consensus.matrix <- res$consensusMatrix
    consensus.class <- res$consensusClass
    upper.value <- consensus.matrix[upper.tri(consensus.matrix)]
    same.module <- outer(consensus.class, consensus.class, FUN = "==")[upper.tri(consensus.matrix)]
    data.frame(
        k = k,
        PAC = calc_pac(consensus.matrix),
        mean_within_consensus = mean(upper.value[same.module], na.rm = TRUE),
        mean_between_consensus = mean(upper.value[!same.module], na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

k.selection <- do.call(rbind, lapply(non.genomic.consensus.k.eval, summarise_consensus_k))
k.selection$within_minus_between <- k.selection$mean_within_consensus - k.selection$mean_between_consensus

if (identical(non.genomic.consensus.k, "auto")) {

    best.pac <- min(k.selection$PAC, na.rm = TRUE)
    selected.k <- min(k.selection$k[k.selection$PAC <= best.pac + non.genomic.consensus.pac.tolerance])
    selection.mode <- "auto_by_min_PAC_choose_smallest"
} else {
    selected.k <- manual.k
    selection.mode <- "manual"
}

k.selection$selected <- k.selection$k == selected.k
utils::write.table(
    k.selection,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".k_selection.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

p.k <- ggplot2::ggplot(k.selection, ggplot2::aes(x = k, y = PAC)) +
    ggplot2::geom_line(color = "grey40", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(fill = selected), shape = 21, size = 3, color = "black") +
    ggplot2::scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#c82621"), guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::labs(
        x = "Number of consensus modules (k)",
        y = "PAC score",
        title = paste0("Selected k = ", selected.k, " (", selection.mode, ")")
    )
ggplot2::ggsave(
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".PAC_k_selection.pdf")),
    p.k, width = 6, height = 4
)

selected.result <- non.genomic.consensus.result[[selected.k]]
selected.class <- selected.result$consensusClass
selected.features <- names(selected.class)
input.feature.names <- colnames(non.genomic.consensus.matrix)
if (is.null(selected.features) || length(selected.features) != length(input.feature.names) ||
    !all(selected.features %in% input.feature.names)) {
    selected.features <- input.feature.names
    names(selected.class) <- selected.features
}
feature.order <- selected.features
if (!is.null(selected.result$consensusTree)) {
    tree.labels <- selected.result$consensusTree$labels[selected.result$consensusTree$order]
    if (length(tree.labels) == length(selected.features) && all(tree.labels %in% selected.features)) {
        feature.order <- tree.labels
    } else {
        tree.index <- suppressWarnings(as.integer(tree.labels))
        if (length(tree.index) == length(selected.features) &&
            all(is.finite(tree.index)) &&
            all(tree.index >= 1 & tree.index <= length(selected.features))) {
            feature.order <- selected.features[tree.index]
        } else if (length(selected.result$consensusTree$order) == length(selected.features)) {
            feature.order <- selected.features[selected.result$consensusTree$order]
        }
    }
}

module.levels <- paste0("M", sort(unique(as.integer(selected.class))))
selected.module <- paste0("M", as.integer(selected.class))
selected.module <- factor(selected.module, levels = module.levels)
module.size <- table(selected.module)

feature.module.table <- data.frame(
    feature = selected.features,
    original_category = unname(feature.category[selected.features]),
    original_category_color = unname(feature.color[selected.features]),
    consensus_module = as.character(selected.module),
    module_size = as.integer(module.size[as.character(selected.module)]),
    feature_order = match(selected.features, feature.order),
    selected_k = selected.k,
    selection_mode = selection.mode,
    PAC_at_selected_k = k.selection$PAC[match(selected.k, k.selection$k)],
    stringsAsFactors = FALSE
)
feature.module.table <- feature.module.table[order(feature.module.table$consensus_module, feature.module.table$feature_order), ]

module.summary <- do.call(rbind, lapply(split(feature.module.table, feature.module.table$consensus_module), function(x) {
    data.frame(
        consensus_module = unique(x$consensus_module),
        module_size = nrow(x),
        original_categories = paste(unique(x$original_category), collapse = ";"),
        features = paste(x$feature, collapse = ";"),
        stringsAsFactors = FALSE
    )
}))

utils::write.table(
    feature.module.table,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".feature_module_appendix.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)
utils::write.table(
    module.summary,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".module_summary.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

consensus.matrix.selected <- selected.result$consensusMatrix
if (is.null(rownames(consensus.matrix.selected)) ||
    is.null(colnames(consensus.matrix.selected)) ||
    !all(feature.order %in% rownames(consensus.matrix.selected)) ||
    !all(feature.order %in% colnames(consensus.matrix.selected))) {
    dimnames(consensus.matrix.selected) <- list(selected.features, selected.features)
}
consensus.matrix.selected <- consensus.matrix.selected[feature.order, feature.order]
annotation.df <- data.frame(
    ConsensusModule = feature.module.table$consensus_module[match(feature.order, feature.module.table$feature)],
    OriginalCategory = feature.module.table$original_category[match(feature.order, feature.module.table$feature)],
    stringsAsFactors = FALSE
)
rownames(annotation.df) <- feature.order
annotation.df$ConsensusModule[is.na(annotation.df$ConsensusModule) | annotation.df$ConsensusModule == ""] <- "Unknown"
annotation.df$OriginalCategory[is.na(annotation.df$OriginalCategory) | annotation.df$OriginalCategory == ""] <- "Unknown"

module.levels.plot <- sort(unique(annotation.df$ConsensusModule))
module.colors <- setNames(grDevices::hcl.colors(length(module.levels.plot), palette = "Dark 3"), module.levels.plot)
category.levels.plot <- names(features)[names(features) %in% unique(as.character(annotation.df$OriginalCategory))]
category.levels.plot <- c(category.levels.plot, setdiff(unique(as.character(annotation.df$OriginalCategory)), category.levels.plot))
annotation.df$OriginalCategory <- factor(as.character(annotation.df$OriginalCategory), levels = category.levels.plot)
category.colors.plot <- category.colors[intersect(category.levels.plot, names(category.colors))]
unknown.category <- setdiff(category.levels.plot, names(category.colors.plot))
if (length(unknown.category) > 0) {
    category.colors.plot <- c(category.colors.plot, setNames(rep("#BDBDBD", length(unknown.category)), unknown.category))
}
category.colors.plot <- category.colors.plot[category.levels.plot]
annotation.colors <- list(
    ConsensusModule = module.colors[module.levels.plot],
    OriginalCategory = category.colors.plot
)

if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(
        consensus.matrix.selected,
        color = grDevices::colorRampPalette(c("#F7F7F7", "#2166AC"))(101),
        breaks = seq(0, 1, length.out = 102),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        annotation_row = annotation.df,
        annotation_col = annotation.df,
        annotation_colors = annotation.colors,
        border_color = NA,
        fontsize_row = 7,
        fontsize_col = 7,
        main = paste0("Non-genomic escape features consensus matrix, k = ", selected.k),
        filename = file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".selected_k_consensus_heatmap.pdf")),
        width = 12,
        height = 10
    )

    spearman.cor.selected <- suppressWarnings(stats::cor(
        non.genomic.consensus.matrix[, feature.order, drop = FALSE],
        method = "spearman",
        use = "pairwise.complete.obs"
    ))
    spearman.cor.selected[!is.finite(spearman.cor.selected)] <- 0
    spearman.cor.selected <- (spearman.cor.selected + t(spearman.cor.selected)) / 2
    diag(spearman.cor.selected) <- 1
    spearman.cor.selected <- spearman.cor.selected[feature.order, feature.order]
    pheatmap::pheatmap(
        spearman.cor.selected,
        color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
        breaks = seq(-1, 1, length.out = 102),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        annotation_row = annotation.df,
        annotation_col = annotation.df,
        annotation_colors = annotation.colors,
        border_color = NA,
        fontsize_row = 7,
        fontsize_col = 7,
        main = "Spearman correlation among clustered non-genomic features",
        filename = file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".selected_features_spearman_heatmap.pdf")),
        width = 12,
        height = 10
    )
} else {
    warning("pheatmap is not installed; skipping selected-k consensus and Spearman heatmaps.")
}

message("Selected non-genomic consensus clustering k: ", selected.k, " (", selection.mode, ")")
print(module.summary)

calc_module_eigenscore <- function(x) {
    x <- as.data.frame(x, check.names = FALSE)
    n.sample <- nrow(x)
    x <- data.frame(lapply(x, function(v) suppressWarnings(as.numeric(as.character(v)))),
                    check.names = FALSE)

    keep <- vapply(x, function(v) {
        sum(is.finite(v)) >= non.genomic.consensus.min.finite &&
            is.finite(stats::sd(v, na.rm = TRUE)) &&
            stats::sd(v, na.rm = TRUE) > 0
    }, logical(1))
    x <- x[, keep, drop = FALSE]
    if (ncol(x) == 0) {
        return(rep(NA_real_, n.sample))
    }
    if (ncol(x) == 1) {
        score <- as.numeric(scale(x[[1]]))
        return(score)
    }

    x.scaled <- scale(as.matrix(x))
    x.scaled[!is.finite(x.scaled)] <- NA
    for (j in seq_len(ncol(x.scaled))) {
        missing.idx <- !is.finite(x.scaled[, j])
        if (any(missing.idx)) {
            x.scaled[missing.idx, j] <- stats::median(x.scaled[, j], na.rm = TRUE)
        }
    }

    pc <- stats::prcomp(x.scaled, center = FALSE, scale. = FALSE)
    score <- pc$x[, 1]
    module.mean <- rowMeans(x.scaled, na.rm = TRUE)
    score.cor <- suppressWarnings(stats::cor(score, module.mean, use = "pairwise.complete.obs"))
    if (is.finite(score.cor) && score.cor < 0) {
        score <- -score
    }
    as.numeric(scale(score))
}

module.feature.list <- split(feature.module.table$feature, feature.module.table$consensus_module)
module.feature.list <- module.feature.list[module.levels]

module.eigenscore.matrix <- do.call(cbind, lapply(names(module.feature.list), function(module.id) {
    module.features <- module.feature.list[[module.id]]
    module.features <- module.features[module.features %in% colnames(non.genomic.consensus.matrix)]
    calc_module_eigenscore(non.genomic.consensus.matrix[, module.features, drop = FALSE])
}))
colnames(module.eigenscore.matrix) <- names(module.feature.list)
rownames(module.eigenscore.matrix) <- rownames(non.genomic.consensus.matrix)

module.eigenscore.df <- data.frame(
    SampleID = rownames(module.eigenscore.matrix),
    module.eigenscore.matrix,
    check.names = FALSE
)
utils::write.table(
    module.eigenscore.df,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".module_eigenscores.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

module.eigenscore.cor <- suppressWarnings(stats::cor(
    module.eigenscore.matrix,
    method = "spearman",
    use = "pairwise.complete.obs"
))
module.eigenscore.cor[!is.finite(module.eigenscore.cor)] <- NA_real_
diag(module.eigenscore.cor) <- 1

module.eigenscore.p <- matrix(NA_real_, nrow = ncol(module.eigenscore.matrix), ncol = ncol(module.eigenscore.matrix))
rownames(module.eigenscore.p) <- colnames(module.eigenscore.matrix)
colnames(module.eigenscore.p) <- colnames(module.eigenscore.matrix)

module.pair.table <- do.call(rbind, lapply(utils::combn(colnames(module.eigenscore.matrix), 2, simplify = FALSE), function(pair) {
    x <- module.eigenscore.matrix[, pair[1]]
    y <- module.eigenscore.matrix[, pair[2]]
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 3) {
        rho <- NA_real_
        p.value <- NA_real_
    } else {
        test <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
        rho <- unname(test$estimate)
        p.value <- test$p.value
    }
    data.frame(
        module_a = pair[1],
        module_b = pair[2],
        n_pairwise_samples = sum(keep),
        spearman_rho = rho,
        p_value = p.value,
        abs_rho = abs(rho),
        stringsAsFactors = FALSE
    )
}))
module.pair.table$FDR <- stats::p.adjust(module.pair.table$p_value, method = "BH")
for (i in seq_len(nrow(module.pair.table))) {
    module.eigenscore.p[module.pair.table$module_a[i], module.pair.table$module_b[i]] <- module.pair.table$p_value[i]
    module.eigenscore.p[module.pair.table$module_b[i], module.pair.table$module_a[i]] <- module.pair.table$p_value[i]
}
diag(module.eigenscore.p) <- 0

module.relationship.summary <- data.frame(
    selected_k = selected.k,
    n_modules = ncol(module.eigenscore.matrix),
    n_module_pairs = nrow(module.pair.table),
    max_abs_rho = max(module.pair.table$abs_rho, na.rm = TRUE),
    median_abs_rho = stats::median(module.pair.table$abs_rho, na.rm = TRUE),
    n_pairs_abs_rho_ge_0.3 = sum(module.pair.table$abs_rho >= 0.3, na.rm = TRUE),
    n_pairs_abs_rho_ge_0.5 = sum(module.pair.table$abs_rho >= 0.5, na.rm = TRUE),
    n_pairs_FDR_lt_0.05 = sum(module.pair.table$FDR < 0.05, na.rm = TRUE),
    interpretation = ifelse(
        max(module.pair.table$abs_rho, na.rm = TRUE) < 0.3,
        "weak_module_level_correlations_relatively_distinct_axes",
        ifelse(
            max(module.pair.table$abs_rho, na.rm = TRUE) < 0.5,
            "moderate_module_level_correlations_partial_redundancy",
            "substantial_module_level_correlations_not_independent_axes"
        )
    ),
    stringsAsFactors = FALSE
)

utils::write.table(
    module.eigenscore.cor,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".module_eigenscore_correlation_matrix.tsv")),
    sep = "\t", col.names = NA, quote = FALSE
)
utils::write.table(
    module.eigenscore.p,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".module_eigenscore_pvalue_matrix.tsv")),
    sep = "\t", col.names = NA, quote = FALSE
)
utils::write.table(
    module.pair.table,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".module_eigenscore_pairwise_correlations.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)
utils::write.table(
    module.relationship.summary,
    file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".module_relationship_summary.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

module.annotation.df <- module.summary[, c("consensus_module", "module_size"), drop = FALSE]
module.annotation.df <- module.annotation.df[match(colnames(module.eigenscore.matrix), module.annotation.df$consensus_module), , drop = FALSE]
rownames(module.annotation.df) <- module.annotation.df$consensus_module
module.annotation.df$consensus_module <- NULL

if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(
        module.eigenscore.cor,
        color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
        breaks = seq(-1, 1, length.out = 102),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_row = module.annotation.df,
        annotation_col = module.annotation.df,
        border_color = NA,
        display_numbers = round(module.eigenscore.cor, 2),
        number_color = "black",
        fontsize_number = 8,
        main = "Spearman correlations among non-genomic consensus module eigenscores",
        filename = file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".module_eigenscore_correlation_heatmap.pdf")),
        width = 7,
        height = 6
    )
} else {
    warning("pheatmap is not installed; skipped module eigenscore correlation heatmap.")
}

message("Non-genomic consensus module eigenscore relationship summary:")
print(module.relationship.summary)

if (requireNamespace("cluster", quietly = TRUE)) {

    silhouette.feature.order <- colnames(non.genomic.consensus.matrix)
    silhouette.cor <- suppressWarnings(stats::cor(
        non.genomic.consensus.matrix[, silhouette.feature.order, drop = FALSE],
        method = "spearman",
        use = "pairwise.complete.obs"
    ))
    silhouette.cor[!is.finite(silhouette.cor)] <- 0
    silhouette.cor <- (silhouette.cor + t(silhouette.cor)) / 2
    diag(silhouette.cor) <- 1
    silhouette.cor[silhouette.cor > 1] <- 1
    silhouette.cor[silhouette.cor < -1] <- -1
    silhouette.dist <- stats::as.dist(1 - silhouette.cor)

    get_consensus_class_for_k <- function(k) {
        current.class <- non.genomic.consensus.result[[k]]$consensusClass
        if (is.null(names(current.class)) ||
            length(current.class) != length(silhouette.feature.order) ||
            !all(names(current.class) %in% silhouette.feature.order)) {
            names(current.class) <- silhouette.feature.order
        }
        current.class[silhouette.feature.order]
    }

    calc_silhouette_for_k <- function(k) {
        current.class <- get_consensus_class_for_k(k)
        current.sil <- cluster::silhouette(as.integer(current.class), silhouette.dist)
        current.width <- current.sil[, "sil_width"]
        data.frame(
            k = k,
            mean_silhouette_width = mean(current.width, na.rm = TRUE),
            median_silhouette_width = stats::median(current.width, na.rm = TRUE),
            min_silhouette_width = min(current.width, na.rm = TRUE),
            n_negative_silhouette = sum(current.width < 0, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }

    silhouette.k.selection <- do.call(
        rbind,
        lapply(non.genomic.consensus.k.eval, calc_silhouette_for_k)
    )
    silhouette.k.selection$selected_by_silhouette <- FALSE
    best.silhouette <- max(silhouette.k.selection$mean_silhouette_width, na.rm = TRUE)
    selected.k.silhouette <- min(silhouette.k.selection$k[
        silhouette.k.selection$mean_silhouette_width == best.silhouette
    ])
    silhouette.k.selection$selected_by_silhouette <- silhouette.k.selection$k == selected.k.silhouette
    silhouette.k.selection$selected_by_PAC <- silhouette.k.selection$k == selected.k

    utils::write.table(
        silhouette.k.selection,
        file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".silhouette_k_selection.tsv")),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    p.silhouette <- ggplot2::ggplot(
        silhouette.k.selection,
        ggplot2::aes(x = k, y = mean_silhouette_width)
    ) +
        ggplot2::geom_line(color = "grey40", linewidth = 0.5) +
        ggplot2::geom_point(ggplot2::aes(fill = selected_by_silhouette), shape = 21, size = 3, color = "black") +
        ggplot2::geom_vline(xintercept = selected.k, linetype = "dashed", color = "#c82621", linewidth = 0.4) +
        ggplot2::scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#c82621"), guide = "none") +
        ggplot2::theme_bw() +
        ggplot2::labs(
            x = "Number of consensus modules (k)",
            y = "Mean silhouette width",
            title = paste0("Silhouette-selected k = ", selected.k.silhouette)
        )

    ggplot2::ggsave(
        file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".silhouette_k_selection.pdf")),
        p.silhouette,
        width = 6,
        height = 4
    )

    silhouette.selected.class <- get_consensus_class_for_k(selected.k.silhouette)
    silhouette.selected <- cluster::silhouette(as.integer(silhouette.selected.class), silhouette.dist)
    silhouette.selected.features <- silhouette.feature.order
    silhouette.feature.table <- data.frame(
        feature = silhouette.selected.features,
        original_category = unname(feature.category[silhouette.selected.features]),
        silhouette_selected_module = paste0("M", as.integer(silhouette.selected.class[silhouette.selected.features])),
        silhouette_width = silhouette.selected[, "sil_width"],
        selected_k_by_silhouette = selected.k.silhouette,
        selected_k_by_PAC = selected.k,
        stringsAsFactors = FALSE
    )

    utils::write.table(
        silhouette.feature.table,
        file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".silhouette_feature_widths.selected_k.tsv")),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    silhouette.summary <- data.frame(
        selected_k_by_PAC = selected.k,
        selected_k_by_silhouette = selected.k.silhouette,
        PAC_and_silhouette_agree = selected.k == selected.k.silhouette,
        max_mean_silhouette_width = best.silhouette,
        stringsAsFactors = FALSE
    )
    utils::write.table(
        silhouette.summary,
        file.path(non.genomic.consensus.output.dir, paste0(non.genomic.consensus.file.prefix, ".silhouette_selection_summary.tsv")),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    message("Alternative k selection by average silhouette width:")
    print(silhouette.summary)
} else {
    warning("cluster package is not installed; skipped Section 6 silhouette-based k selection.")
}
