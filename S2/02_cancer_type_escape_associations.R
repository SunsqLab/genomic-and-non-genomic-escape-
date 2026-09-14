# S2: Co-occurrence of immune escape features
# In RStudio, set the working directory to this file's directory, edit the paths,
# then Source this file or run it from top to bottom. The functions are sourced below.
# From another directory: source("path/to/02_cancer_type_escape_associations.R", chdir = TRUE).
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
library(igraph)
library(corrplot)
library(dplyr)

Dir.output="results"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive = TRUE) }
Dir.output <- normalizePath(Dir.output, mustWork = TRUE)

source(file.path(functions.dir, "DatasetLabels.R"))
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA'))
source(file.path(functions.dir, "CombineData.XTeam.R"))
profile.escape=CombineData.XTeam(ID.xteams, file.query=file.path(escape.query, "ProfilingImmuneEvading.rds"))
B2M.loh=CombineData.XTeam(ID.xteams, file.query=b2m.query)

source(file.path(functions.dir, "GetInfor.PatientCenter.R"))

patient.center=readRDS(file.path(PanCancer.path, "PatientCenter", "PatientCenter.rds"))

# Primary samples; preserve the source cancer-specific B2M replacement rules.
profile.escape=lapply(profile.escape, function(x){
            tmp=lapply(x, function(y){
                        y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
                    })
            tmp=Reduce(function(x, y) dplyr::full_join(x, y, by = "SampleID"), tmp)
            tmp=tmp    %>%
                dplyr::mutate(all.APM.mut=ifelse(all.APM.mut>0, TRUE, FALSE), HLA.mut=ifelse(HLA.mut>0, TRUE, FALSE))
            tmp=Filter(function(x) !all(is.na(x)), tmp)
        })      %>%     setNames(ID.xteams)

profile.escape <- mapply(function(df, b2m.mat){
    df$B2M.LOH <- NA
    if(!is.null(b2m.mat) && "B2M" %in% rownames(b2m.mat)){
        idx <- match(df$SampleID, colnames(b2m.mat))
        matched <- !is.na(idx)
        df$B2M.LOH[matched] <- b2m.mat["B2M", idx[matched]] == 1
    }
    df
}, profile.escape, B2M.loh[names(profile.escape)], SIMPLIFY = FALSE)

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

# Genomic pairs: direction-selected one-sided Fisher tests; BH within each cancer.
source(file.path(functions.dir, "EscapePairFisher.R"))
temp.genome.features=genome.features[!genome.features%in%c("HLA.mut", "B2M.mut")]
stat.data=EscapePairFisher(data=profile.escape, escape.feature=temp.genome.features)
fisher.count.cols <- c("Neither", "A_not_B", "B_not_A", "Both")
for(col in setdiff(fisher.count.cols, colnames(stat.data))){
    stat.data[[col]] <- NA
}
supp.genome.genome <- stat.data %>%
    dplyr::transmute(
        featureA = geneA,
        featureB = geneB,
        Neither = Neither,
        A_not_B = A_not_B,
        B_not_A = B_not_A,
        Both = Both,
        `Odds Ratio` = odds_ratio,
        FDR = p_adj,
        Cancer = cancer
    )
write.csv(supp.genome.genome, file=file.path(Dir.output, "S3.G-G features association.csv"), row.names=FALSE)

data=stat.data     %>%
    dplyr::group_by(geneA, geneB)       %>%
    dplyr::reframe( positive=sum(p_adj<0.05 & odds_ratio>1, na.rm=T), negative=sum(p_adj<0.05 & odds_ratio<1, na.rm=T),
                    pos_pvalue_cancers = paste(unique(na.omit(cancer[p_adj < 0.05 & odds_ratio > 1])),collapse = ","),
                    neg_pvalue_cancers = paste(unique(na.omit(cancer[p_adj < 0.05 & odds_ratio < 1])),collapse = ","))      %>%
        as.data.frame()

data=reshape2::melt(data, id.vars = c("geneA", "geneB","pos_pvalue_cancers","neg_pvalue_cancers"), variable.name = "significance", value.name = "pvalue_freq")
data <- data %>%
    dplyr::mutate(
        pvalue_cancers = ifelse(significance == "positive", pos_pvalue_cancers, neg_pvalue_cancers),
        geneA.type = names(all.features)[match(geneA, all.features)],
        geneB.type = names(all.features)[match(geneB, all.features)]
    ) %>%
    dplyr::select(-pos_pvalue_cancers, -neg_pvalue_cancers)

network.data=data[which(data$pvalue_freq>1), ]

nodes.data=data.frame(nodes=unique(c(network.data$geneA, network.data$geneB)))
nodes.data$fea.type=names(all.features)[match(nodes.data$nodes, all.features)]

net=graph_from_data_frame(network.data, directed=FALSE, vertices=nodes.data)
write_graph(net, file = file.path(Dir.output, "1.genomic_escape_associations.by_cancer.graphml"), format = "graphml")

# Panel-wide complete cases; preserve the original BH family.
source(file.path(functions.dir, "EscapePairSpearman.R"))
stat.data=EscapePairSpearman(data=profile.escape, escape.feature=non.genomic)
data=stat.data     %>%
        dplyr::group_by(feaA, feaB, label)       %>%
        dplyr::reframe(mean_cor=mean(cor, na.rm=TRUE),

            positive=sum(p_adj<0.05 & cor>0, na.rm=T), negative=sum(p_adj<0.05 & cor<0, na.rm=T),
            pos_pvalue_cancers = paste(unique(na.omit(cancer[p_adj < 0.05 & cor > 0])),collapse = ","),
            neg_pvalue_cancers = paste(unique(na.omit(cancer[p_adj < 0.05 & cor < 0])),collapse = ","))     %>%
        as.data.frame()
data=reshape2::melt(data, id.vars = c("feaA", "feaB", "mean_cor", "label","pos_pvalue_cancers","neg_pvalue_cancers"), variable.name = "significance", value.name = "pvalue_freq")

data <- data %>%
    dplyr::group_by(feaA, feaB) %>%
    dplyr::slice_max(pvalue_freq, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        pvalue_cancers = ifelse(significance == "positive", pos_pvalue_cancers, neg_pvalue_cancers)
    ) %>%
    dplyr::select(-pos_pvalue_cancers, -neg_pvalue_cancers)

supp.nongenome.nongenome <- data %>%
    dplyr::transmute(
        featureA = as.character(feaA),
        featureB = as.character(feaB),
        mean_cor = mean_cor,
        label = label,
        significance = as.character(significance),
        pvalue_freq = pvalue_freq,
        pvalue_cancers = pvalue_cancers
    )
write.csv(supp.nongenome.nongenome, file=file.path(Dir.output, "S5.NG-NG features association.csv"), row.names=FALSE)

data.cor <- data %>%
  dplyr::mutate(
    cor = ifelse(
      significance == "positive",  pvalue_freq,
      ifelse(significance == "negative", -pvalue_freq, NA)
    )
  ) %>%
  dplyr::select(feaA, feaB, cor) %>%
  tidyr::pivot_wider(names_from = feaB, values_from = cor) %>%
  tibble::column_to_rownames("feaA") %>%
  as.matrix()

var.order <- non.genomic
data.cor <- data.cor[var.order, var.order]

num.mat <- round(abs(data.cor))
num.mat[is.na(data.cor)] <- NA

tl.col <- names(all.features)[match(var.order, all.features)]
col.RWB <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(200)

pdf(file.path(Dir.output, "2.non_genomic_escape_associations.by_cancer.pdf"),width = 16, height = 12)
corrplot(data.cor,method = "square",col = col.RWB,
  outline = "grey",order = "original",diag = TRUE,tl.cex = 1.2,tl.col = tl.col,tl.srt = 45,
  addgrid.col = "grey",type = "lower",
  na.label = " ", is.corr = FALSE
)

n <- nrow(data.cor)
for (i in seq_len(n)) {
  for (j in seq_len(n)) {

    if (i > j && !is.na(num.mat[i, j]) && num.mat[i, j] > 0) {
      text(x = j,y = n - i+1,labels = num.mat[i, j],cex = 0.7,col = "black")}
  }}
dev.off()

# Pairwise complete cases; at least 5 samples per group; BH within each cancer.
source(file.path(functions.dir, "EscapePairWilcoxon.R"))
stat.data=EscapePairWilcoxon(data=profile.escape, disc.feature=unlist(genome.features), cont.feature=non.genomic)

data=stat.data     %>%
        dplyr::group_by(feaA, feaB)       %>%
        dplyr::reframe(pos=sum(p_adj<0.05 & fold_change>1, na.rm=T), neg=sum(p_adj<0.05 & fold_change<1, na.rm=T),
                        pos_pvalue_cancers = paste(unique(na.omit(cancer[p_adj < 0.05 & fold_change > 1])),collapse = ","),
                        neg_pvalue_cancers = paste(unique(na.omit(cancer[p_adj < 0.05 & fold_change < 1])),collapse = ","))     %>%
        as.data.frame()
data=reshape2::melt(data, id.vars = c("feaA", "feaB","pos_pvalue_cancers","neg_pvalue_cancers"), variable.name = "significance", value.name = "pvalue_freq")

data <- data %>%

    dplyr::group_by(feaA, feaB) %>%
    dplyr::slice_max(pvalue_freq, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        pvalue_cancers = ifelse(significance == "pos", pos_pvalue_cancers, neg_pvalue_cancers),
        feaA = factor(feaA, levels = genome.features),
        feaB = factor(feaB, levels = rev(non.genomic)),
        color.map = ifelse(significance == "pos",pvalue_freq ,-pvalue_freq ),
        feaA.class = factor(rep(names(features), lengths(features))[match(feaA, unlist(features))], levels = c("APMalt", "ActMalt", "Checkalt")),
        feaB.class = factor(rep(names(features), lengths(features))[match(feaB, unlist(features))], levels = c("Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo")),
    )%>%
    dplyr::select(-pos_pvalue_cancers, -neg_pvalue_cancers)

supp.genome.nongenome <- data %>%
    dplyr::transmute(
        featureA = as.character(feaA),
        featureB = as.character(feaB),
        significance = as.character(significance),
        pvalue_freq = pvalue_freq,
        pvalue_cancers = pvalue_cancers,
        featureA.class = as.character(feaA.class),
        featureB.class = as.character(feaB.class)
    )
write.csv(supp.genome.nongenome, file=file.path(Dir.output, "S4.G-NG features association.csv"), row.names=FALSE)
feature_col_map <- setNames(names(all.features), all.features)
p <- ggplot(data, aes(x = feaA, y = feaB, fill = color.map)) +
    geom_tile(color = "grey85", size = 0.3) +
    geom_text(aes(label = pvalue_freq), size = 3, color = "black") +
    scale_fill_gradient2(
        low = "#4B71B3", mid = "white", high = "#D66B6B",
        midpoint = 0, na.value = "white", name = "No. significant cancers",
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
ggsave(file.path(Dir.output, "3.genomic_non_genomic_associations.by_cancer.pdf"), p, width = 4, height = 9)
