# S5: Clonal-subclonal mutation analysis
# This script evaluates selective pressure on immune escape alterations using TCGA
# mutation or copy-number data. It performs the indicated stratified analyses and
# saves gene- or feature-level statistical summaries and figures.

DataCenter.file <- "/path/to/DataCenter.rds"
PanCancer.path <- "/path/to/PanCancer_TCGA.dataset"
escape.query <- "Results/BioImmune/ImmuneEscape"
b2m.query <- "Results/BioGenomics/B2M_biallelic_inactivation_matrix.rds"
functions.dir <- normalizePath("functions", mustWork = TRUE)
exclusion.type <- "cna"

library(dplyr)
library(ggplot2)

out.dir <- "results"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
out.dir <- normalizePath(out.dir, mustWork = TRUE)

source(file.path(functions.dir, "DatasetLabels.R"))
source(file.path(functions.dir, "CombineData.XTeam.R"))
source(file.path(functions.dir, "GetInfor.PatientCenter.R"))
ID.xteam <- "PanCancer_TCGA.dataset"

ID.xteams <- get(ID.xteam)

interested.metric <- c("SampleType", "CancerType", "TMB", "TMB.subclone", "TMB.clone", "Neoantigen.burden.clone", "Neoantigen.burden.subclone", "Liuwei_2025", "Purity", "Ploidy", "MSISubtype2")

patient.center <- readRDS(file.path(PanCancer.path, "PatientCenter", "PatientCenter.rds"))
patient.center <- GetInfor.PatientCenter(patient.center, colNames = interested.metric)
patient.center <- patient.center[!is.na(patient.center$SampleType) & grepl("Primary", patient.center$SampleType), ]
rownames(patient.center) <- NULL

patient.center$CancerType <- ifelse(patient.center$CancerType %in% c("COAD", "READ"), "CRC", patient.center$CancerType)
patient.center$MSI.Subtype <- ifelse(patient.center$MSISubtype2 == "MSI-H", "MSI", ifelse(patient.center$MSISubtype2 %in% c("MSS", "MSI-L"), "MSS", NA))

all.features.list <- list(
    c("HLA.LOH", "B2M.LOH", "HLA.mut", "B2M.mut", "all.APM.mut"),
    c("CD274.deepAmp"),
    c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
    c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),

    "mean.neo.exprs"
)
all.features <- unlist(all.features.list)
names(all.features) <- rep(c("APMalt", "Checkalt", "ActMalt", "Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo"), lengths(all.features.list))

B2M.loh <- CombineData.XTeam(ID.xteams, file.query = b2m.query)
B2M.loh.df <- dplyr::bind_rows(lapply(B2M.loh, function(x) {
    data.frame(
        SampleID = colnames(x),
        B2M.LOH = x["B2M", ] == 1,
        stringsAsFactors = FALSE
    )
}))

profile.escape <- CombineData.XTeam(ID.xteams, file.query = file.path(escape.query, "ProfilingImmuneEvading.Binarization.rds"))

profile.escape <- dplyr::bind_rows(lapply(names(profile.escape), function(name.1) {
    purrr::reduce(profile.escape[[name.1]], full_join, by = "SampleID")
}))%>%
    dplyr::select(-dplyr::any_of("B2M.LOH")) %>%
    dplyr::left_join(B2M.loh.df, by = "SampleID")

patient.escape.profiles.df <- dplyr::inner_join(patient.center, profile.escape, by = "SampleID") %>%
    dplyr::bind_rows(dplyr::mutate(., CancerType = "PanCancer")) %>%
    data.frame()

tcga.mut <- list()
tcga.mut$PanCancer_TCGA <- readRDS(file.path(PanCancer.path, "OMICSData", "Mutations.data.rds"))

tcga.mut <- purrr::map(tcga.mut, function(cancer) {
    cancer <- cancer %>%
        dplyr::filter(SampleID %in% patient.escape.profiles.df$SampleID)
})

Cancer.order <- c("PanCancer", "SKCM", "LUSC", "LUAD", "BLCA", "DLBC", "CRC", "STAD", "ESCA", "HNSC", "CESC", "LIHC", "UCEC", "OV", "KIRP", "KIRC", "GBM", "UCS", "SARC", "BRCA", "PAAD", "CHOL", "ACC", "LGG", "MESO", "PRAD", "KICH", "TGCT", "THYM", "UVM", "THCA", "PCPG")

source(file.path(functions.dir, "CloneTmb.Diff.R"))

data <- patient.escape.profiles.df[patient.escape.profiles.df$Liuwei_2025 == "PASS", ]
Clonal.tmb.diff <- CloneTmb.Diff(data, all.features, "CancerType", "TMB.clone")

Clonal.tmb.diff$CancerType <- factor(Clonal.tmb.diff$CancerType, levels = Cancer.order)
SubClonal.tmb.diff <- CloneTmb.Diff(data, all.features, "CancerType", "TMB.subclone")
SubClonal.tmb.diff$CancerType <- factor(SubClonal.tmb.diff$CancerType, levels = Cancer.order)

source(file.path(functions.dir, "PlotRectHeatmapWithSignif.R"))
Clonal.tmb.diff.genomic <- Clonal.tmb.diff[Clonal.tmb.diff$category %in% c("APMalt", "ActMalt", "Checkalt"), ]
SubClonal.tmb.diff.genomic <- SubClonal.tmb.diff[SubClonal.tmb.diff$category %in% c("APMalt", "ActMalt", "Checkalt"), ]
p1 <- PlotRectHeatmapWithSignif(Clonal.tmb.diff.genomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TMB.clone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)
p2 <- PlotRectHeatmapWithSignif(SubClonal.tmb.diff.genomic, x = "CancerType", y = "EscapeMetric", statistic.p = "p.value", tag = "TMB.subclone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)

plot.list <- patchwork::wrap_plots(p1, p2, ncol = 1)
ggsave(file.path(out.dir, "TMB_clonality.genomic_escape.pdf"), plot.list, width = 15, height = 10)

Clonal.tmb.diff.nongenomic <- Clonal.tmb.diff[!Clonal.tmb.diff$category %in% c("APMalt", "ActMalt", "Checkalt"), ]
SubClonal.tmb.diff.nongenomic <- SubClonal.tmb.diff[!SubClonal.tmb.diff$category %in% c("APMalt", "ActMalt", "Checkalt"), ]
p1 <- PlotRectHeatmapWithSignif(Clonal.tmb.diff.nongenomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TMB.clone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)
p2 <- PlotRectHeatmapWithSignif(SubClonal.tmb.diff.nongenomic, x = "CancerType", y = "EscapeMetric", statistic.p = "p.value", tag = "TMB.subclone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)

plot.list <- patchwork::wrap_plots(p1, p2, ncol = 1)
ggsave(file.path(out.dir, "TMB_clonality.non_genomic_escape.pdf"), plot.list, width = 15, height = 20)

data <- patient.escape.profiles.df[patient.escape.profiles.df$Liuwei_2025 == "PASS", ]

source(file.path(functions.dir, "LinearRegression.FC.R"))

escape.features <- all.features

origin.cancer <- data %>%
    dplyr::filter(CancerType != "PanCancer") %>%
    dplyr::distinct(SampleID, CancerType.origin = CancerType)

data.model.long <- data %>%
    dplyr::left_join(origin.cancer, by = "SampleID") %>%
    dplyr::mutate(
        CancerType.origin = factor(CancerType.origin),
        MSI.Subtype = factor(MSI.Subtype, levels = c("MSS", "MSI")),
        log10TMB = log10(TMB + 1),
        log10TMB.clone = log10(TMB.clone + 1),
        log10TMB.subclone = log10(TMB.subclone + 1)
    ) %>%
    tidyr::pivot_longer(
        cols = tidyr::all_of(unname(escape.features)),
        names_to = "EscapeMetric",
        values_to = "EscapeStatus"
    )

outcome.map <- data.frame(
    TMB.type = c("TMB.clone", "TMB.subclone"),
    Outcome = c("log10TMB.clone", "log10TMB.subclone"),
    stringsAsFactors = FALSE
)

RunTmbModelSet <- function(model.data,
                           cancer.types,
                           adjust.vars,
                           analysis.name,
                           min.n.per.group = 5) {

    fit_one <- function(cancer, metric, outcome, tmb.type) {

        model.cols <- c(outcome, "EscapeStatus", adjust.vars)
        current.data <- model.data %>%
            dplyr::filter(CancerType == cancer, EscapeMetric == metric) %>%
            dplyr::select(dplyr::all_of(c(model.cols, "SampleID"))) %>%
            dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::all_of(model.cols))))

        group.count <- table(current.data$EscapeStatus)
        n.escape <- if ("TRUE" %in% names(group.count)) unname(group.count["TRUE"]) else 0
        n.non.escape <- if ("FALSE" %in% names(group.count)) unname(group.count["FALSE"]) else 0

        if (n.escape < min.n.per.group || n.non.escape < min.n.per.group) {
            return(list(result = NULL, full.result = NULL))
        }

        fit.result <- tryCatch(
            LinearRegression.FC(
                data = current.data,
                outcome = outcome,
                exposure = "EscapeStatus",
                adjust.vars = adjust.vars,
                reference.levels = c(EscapeStatus = "FALSE"),
                return.vars = "all",
                outcome.log2 = FALSE
            ),
            error = function(e) {
                structure(list(), error.message = conditionMessage(e))
            }
        )

        if (is.list(fit.result) && length(fit.result) == 0) {
            return(list(result = NULL, full.result = NULL))
        }

        full.result <- fit.result %>%
            dplyr::mutate(
                Analysis = analysis.name,
                Regression.Variable = Variable,
                CancerType = cancer,
                EscapeMetric = metric,
                TMB.type = tmb.type,
                Effect.scale = "linear coefficient on log10-transformed outcome",
                p.value.raw = p.value,
                N.escape = n.escape,
                N.non_escape = n.non.escape,
                Adjust.vars = paste(adjust.vars, collapse = " + "),
                category = names(escape.features)[match(metric, unname(escape.features))]
            ) %>%
            data.frame()

        result <- full.result %>%
            dplyr::filter(Regression.Variable == "EscapeStatus") %>%
            dplyr::mutate(
                Variable = metric,

                fold.change = Estimate,
                FoldChange = 10^Estimate,
                FC.CI.min = 10^CI.min,
                FC.CI.max = 10^CI.max
            ) %>%
            data.frame()

        list(result = result, full.result = full.result)
    }

    model.grid <- expand.grid(
        CancerType = cancer.types,
        EscapeMetric = unname(escape.features),
        Outcome = outcome.map$Outcome,
        stringsAsFactors = FALSE
    ) %>%
        dplyr::left_join(outcome.map, by = "Outcome")

    model.list <- lapply(seq_len(nrow(model.grid)), function(i) {
        fit_one(
            cancer = model.grid$CancerType[i],
            metric = model.grid$EscapeMetric[i],
            outcome = model.grid$Outcome[i],
            tmb.type = model.grid$TMB.type[i]
        )
    })

    list(
        result = dplyr::bind_rows(lapply(model.list, `[[`, "result")),
        full.result = dplyr::bind_rows(lapply(model.list, `[[`, "full.result"))
    )
}

common.cancers <- setdiff(intersect(Cancer.order, unique(data.model.long$CancerType)), c("PanCancer", "CRC", "STAD", "UCEC"))
Adjusted.tmb.common <- RunTmbModelSet(
    model.data = data.model.long,
    cancer.types = common.cancers,
    adjust.vars = c("log10TMB", "Purity", "Ploidy"),
    analysis.name = "Cancer_specific"
)

Adjusted.tmb.PanCancer <- RunTmbModelSet(
    model.data = data.model.long,
    cancer.types = "PanCancer",
    adjust.vars = c("log10TMB", "Purity", "Ploidy", "CancerType.origin"),
    analysis.name = "PanCancer_adjust_CancerType"
)

msi.cancers <- intersect(c("CRC", "STAD", "UCEC"), unique(data.model.long$CancerType))
Adjusted.tmb.MSI <- RunTmbModelSet(
    model.data = data.model.long,
    cancer.types = msi.cancers,
    adjust.vars = c("log10TMB", "Purity", "Ploidy", "MSI.Subtype"),
    analysis.name = "Cancer_specific_adjust_MSI"
)

Adjusted.tmb.regression <- list(
    common = Adjusted.tmb.common,
    PanCancer = Adjusted.tmb.PanCancer,
    MSI = Adjusted.tmb.MSI
)

Adjusted.tmb.model.result <- dplyr::bind_rows(lapply(Adjusted.tmb.regression, `[[`, "result")) %>%
    dplyr::group_by(TMB.type, CancerType) %>%
    dplyr::mutate(
        p.value = p.adjust(p.value.raw, method = "fdr"),
        significance = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        CancerType = factor(CancerType, levels = Cancer.order),
        EscapeMetric = factor(EscapeMetric, levels = rev(unname(escape.features))),
        category = factor(category, levels = unique(names(escape.features)))
    ) %>%
    data.frame()

Adjusted.tmb.model.full.result <- dplyr::bind_rows(lapply(Adjusted.tmb.regression, `[[`, "full.result")) %>%
    dplyr::group_by(TMB.type, CancerType, Variable) %>%
    dplyr::mutate(
        p.value.raw = p.value,
        p.value.fdr = p.adjust(p.value.raw, method = "fdr")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        CancerType = factor(CancerType, levels = Cancer.order),
        EscapeMetric = factor(EscapeMetric, levels = rev(unname(escape.features)))
    ) %>%
    dplyr::transmute(CancerType,TMB.type,EscapeMetric,Model.variable = Variable,Model.term = Term,Estimate,CI.min,CI.max,
        p.value.raw,p.value.fdr,global.p.value,Reference,N,N.escape,N.non_escape,Adjust.vars) %>%
    data.frame()

Clonal.tmb.diff.adjusted <- Adjusted.tmb.model.result %>%
    dplyr::filter(TMB.type == "TMB.clone")
SubClonal.tmb.diff.adjusted <- Adjusted.tmb.model.result %>%
    dplyr::filter(TMB.type == "TMB.subclone")

write.csv(Adjusted.tmb.model.full.result, file = file.path(out.dir, "AdjustedTMB.escape.regression.model.result.all_covariates.csv"), row.names = FALSE)

genomic.category <- c("APMalt", "ActMalt", "Checkalt")
Clonal.tmb.diff.adjusted.genomic <- Clonal.tmb.diff.adjusted[Clonal.tmb.diff.adjusted$category %in% genomic.category, ]
SubClonal.tmb.diff.adjusted.genomic <- SubClonal.tmb.diff.adjusted[SubClonal.tmb.diff.adjusted$category %in% genomic.category, ]
Clonal.tmb.diff.adjusted.nongenomic <- Clonal.tmb.diff.adjusted[!Clonal.tmb.diff.adjusted$category %in% genomic.category, ]
SubClonal.tmb.diff.adjusted.nongenomic <- SubClonal.tmb.diff.adjusted[!SubClonal.tmb.diff.adjusted$category %in% genomic.category, ]

p1 <- PlotRectHeatmapWithSignif(Clonal.tmb.diff.adjusted.genomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TMB.clone_adjusted_beta", size_title = "|beta|", right_plot_type = "proportion", right_plot_xlim = 0.5)
p2 <- PlotRectHeatmapWithSignif(SubClonal.tmb.diff.adjusted.genomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TMB.subclone_adjusted_beta", size_title = "|beta|", right_plot_type = "proportion", right_plot_xlim = 0.5)
plot.list <- patchwork::wrap_plots(p1, p2, ncol = 1)
ggsave(file.path(out.dir, "TMB_clonality.adjusted.genomic_escape.pdf"), plot.list, width = 15, height = 10)

p1 <- PlotRectHeatmapWithSignif(Clonal.tmb.diff.adjusted.nongenomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TMB.clone_adjusted_beta", size_title = "|beta|", right_plot_type = "proportion", right_plot_xlim = 0.5)
p2 <- PlotRectHeatmapWithSignif(SubClonal.tmb.diff.adjusted.nongenomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TMB.subclone_adjusted_beta", size_title = "|beta|", right_plot_type = "proportion", right_plot_xlim = 0.5)
plot.list <- patchwork::wrap_plots(p1, p2, ncol = 1)
ggsave(file.path(out.dir, "TMB_clonality.adjusted.non_genomic_escape.pdf"), plot.list, width = 15, height = 20)

all.APM.gene <- c(

    "HLA-A", "HLA-B", "HLA-C",

    "B2M", "CALR", "TAP1", "TAP2", "TAPBP", "CIITA", "RFX5", "NLRC5",

    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "CANX", "CD4", "CD74", "CD8A", "CD8B", "CREB1", "CTSB", "CTSL", "CTSS", "ERAP1", "ERAP2", "FAS", "HLA-E", "HLA-F", "HLA-G", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B",
    "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "HSPBP1", "IFI30", "IFNG", "IRF1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "KLRD1", "LGMN", "MEX3B", "NFYA", "NFYB", "NFYC", "PDIA3", "PSMA7", "PSMB10", "PSMB11",
    "PSMB6", "PSMB8", "PSMB9", "PSME1", "PSME2", "PSME3", "PSMF1", "RFXANK", "RFXAP", "TNF"
)
IFNG.gene <- c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")
Escape.gene <- unique(c("CD58", all.APM.gene, "IDH1", IFNG.gene))
names(Escape.gene) <- c("CD58", rep("all.APM", times = length(all.APM.gene)), "IDH1", rep("IFNG.pathway", times = length(IFNG.gene)))

tcga.mut$PanCancer_TCGA <- tcga.mut$PanCancer_TCGA[tcga.mut$PanCancer_TCGA$SampleID %in% patient.escape.profiles.df$SampleID[patient.escape.profiles.df$Liuwei_2025 == "PASS"], ]

type.nonsilent <- c(
    "In_Frame_Ins", "In_Frame_Del", "Nonstop_Mutation", "Translation_Start_Site", "Frame_Shift_Ins", "Frame_Shift_Del",
    "Splice_Site", "Nonsense_Mutation", "Missense_Mutation", "Startloss_Mutation"
)

clonal.nonsyn.enrichment.mut.level <- tcga.mut$PanCancer_TCGA %>%

    dplyr::filter(geneSymbol %in% Escape.gene, !is.na(isClonal), Level5 %in% type.nonsilent) %>%

    dplyr::group_by(geneSymbol, isClonal) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%

    tidyr::complete(geneSymbol = Escape.gene, isClonal = c(TRUE, FALSE), fill = list(count = 0)) %>%
    tidyr::pivot_wider(names_from = isClonal, values_from = count) %>%
    dplyr::rename(`Clonal` = `TRUE`, `SubClonal` = `FALSE`) %>%

    dplyr::mutate(
        GeneType = names(Escape.gene)[match(geneSymbol, Escape.gene)]
    ) %>%
    data.frame()

clonal.nonsyn.enrichment.mut.level$GeneType <- factor(clonal.nonsyn.enrichment.mut.level$GeneType, levels = c("IDH1", "IFNG.pathway", "CD58", "all.APM"))

data.plot.2 <- split(clonal.nonsyn.enrichment.mut.level, clonal.nonsyn.enrichment.mut.level$GeneType)
wilcox.test(clonal.nonsyn.enrichment.mut.level$Clonal, clonal.nonsyn.enrichment.mut.level$SubClonal, paired = TRUE, alternative = "greater")

source(file.path(functions.dir, "mirror_bar.R"))
p2 <- mirror_bar(data.plot.2, x = "geneSymbol", y.up = "Clonal", y.down = "SubClonal", title = "Distribution of Clonal and Subclonal Immune Escape Mutations", break_range = c(50, 310), break_scales = 0.1, facet_split = FALSE)
ggsave(filename = file.path(out.dir, "Escape.gene.distribution(mutation.level).pdf"), plot = p2, width = 211.171, height = 88.208, units = "mm")

data <- tcga.mut$PanCancer_TCGA %>%
    dplyr::filter(SampleID %in% patient.escape.profiles.df$SampleID[patient.escape.profiles.df$all.APM.mut | patient.escape.profiles.df$IFNG.pathway.mut | patient.escape.profiles.df$CD58.mut | patient.escape.profiles.df$IDH1.mut]) %>%

    dplyr::mutate(
        group = ifelse(geneSymbol %in% Escape.gene, "target_gene", "background_gene")
    ) %>%

    dplyr::filter(!is.na(isClonal), Level5 %in% type.nonsilent)

source(file.path(functions.dir, "PlotContingency.R"))
p <- plot_contingency_from_df(data, x_col = "group", y_col = "isClonal", title = "Clonal vs Subclonal Mutation Distribution in Escape and Non- Escape Genes", p_display = "raw")

ggsave(file.path(out.dir, "EscapeGene.Clonal.Enrich(mutation.level).pdf"), p)

genomic.metric <- as.character(all.features[names(all.features) %in% c("APMalt", "Checkalt", "ActMalt")])
genomic.mut.metric <- c("HLA.mut", "B2M.mut", "all.APM.mut", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut")
genomic.cna.metric <- c("HLA.LOH", "B2M.LOH", "CD274.deepAmp", "IFNG.pathway.HD", "CD58.HD")

Transcript.metric <- as.character(all.features[names(all.features) %in% c("Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo")])

if (exclusion.type == "all") {
    data <- patient.escape.profiles.df %>%
        dplyr::filter(

            Liuwei_2025 == "PASS",

            if_all(all_of(genomic.metric), ~ !is.na(.) & . == FALSE)
        )
} else if (exclusion.type == "mut") {
    data <- patient.escape.profiles.df %>%
        dplyr::filter(

            Liuwei_2025 == "PASS",

            rowSums(pick(all_of(genomic.mut.metric)), na.rm = TRUE) == 0,

            if_all(all_of(genomic.mut.metric), ~ !is.na(.) & . == FALSE)
        )
} else if (exclusion.type == "cna") {
    data <- patient.escape.profiles.df %>%
        dplyr::filter(

            Liuwei_2025 == "PASS",

            rowSums(pick(all_of(genomic.cna.metric)), na.rm = TRUE) == 0,

            if_all(all_of(genomic.cna.metric), ~ !is.na(.) & . == FALSE)
        )
} else {
    stop("exclusion.type must be one of: all, mut, cna.")
}

Clonal.tmb.diff <- CloneTmb.Diff(data, all.features, "CancerType", "TMB.clone")

Clonal.tmb.diff$CancerType <- factor(Clonal.tmb.diff$CancerType, levels = Cancer.order)
Clonal.tmb.diff <- Clonal.tmb.diff[Clonal.tmb.diff$EscapeMetric %in% Transcript.metric, ]
SubClonal.tmb.diff <- CloneTmb.Diff(data, all.features, "CancerType", "TMB.subclone")
SubClonal.tmb.diff$CancerType <- factor(SubClonal.tmb.diff$CancerType, levels = Cancer.order)
SubClonal.tmb.diff <- SubClonal.tmb.diff[SubClonal.tmb.diff$EscapeMetric %in% Transcript.metric, ]

p1 <- PlotRectHeatmapWithSignif(Clonal.tmb.diff, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TMB.clone_median", right_plot_type = "proportion", right_plot_xlim = 0.4)
p2 <- PlotRectHeatmapWithSignif(SubClonal.tmb.diff, x = "CancerType", y = "EscapeMetric", statistic.p = "p.value", tag = "TMB.subclone_median", right_plot_type = "proportion", right_plot_xlim = 0.4)

plot.list <- patchwork::wrap_plots(p1, p2, ncol = 1)

ggsave(file.path(out.dir, paste0("TMB_clonality.exclude_", exclusion.type, ".pdf")), plot.list, width = 15, height = 20)

data <- patient.escape.profiles.df[patient.escape.profiles.df$Liuwei_2025 == "PASS", ]
Clonal.tnb.diff <- CloneTmb.Diff(data, all.features, "CancerType", "Neoantigen.burden.clone")

Clonal.tnb.diff$CancerType <- factor(Clonal.tnb.diff$CancerType, levels = Cancer.order)
SubClonal.tnb.diff <- CloneTmb.Diff(data, all.features, "CancerType", "Neoantigen.burden.subclone")
SubClonal.tnb.diff$CancerType <- factor(SubClonal.tnb.diff$CancerType, levels = Cancer.order)

Clonal.tnb.diff.genomic <- Clonal.tnb.diff[Clonal.tnb.diff$category %in% c("APMalt", "Checkalt", "ActMalt"), ]
SubClonal.tnb.diff.genomic <- SubClonal.tnb.diff[SubClonal.tnb.diff$category %in% c("APMalt", "Checkalt", "ActMalt"), ]
p1 <- PlotRectHeatmapWithSignif(Clonal.tnb.diff.genomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TNB.clone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)
p2 <- PlotRectHeatmapWithSignif(SubClonal.tnb.diff.genomic, x = "CancerType", y = "EscapeMetric", statistic.p = "p.value", tag = "TNB.subclone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)

plot.list <- patchwork::wrap_plots(p1, p2, ncol = 1)
ggsave(file.path(out.dir, "TNB_clonality.genomic_escape.pdf"), plot.list, width = 15, height = 10)

Clonal.tnb.diff.nongenomic <- Clonal.tnb.diff[!Clonal.tnb.diff$category %in% c("APMalt", "Checkalt", "ActMalt"), ]
SubClonal.tnb.diff.nongenomic <- SubClonal.tnb.diff[!SubClonal.tnb.diff$category %in% c("APMalt", "Checkalt", "ActMalt"), ]
p1 <- PlotRectHeatmapWithSignif(Clonal.tnb.diff.nongenomic, x = "CancerType", y = "EscapeMetric", statistic.type = "fold.change", statistic.p = "p.value", tag = "TNB.clone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)
p2 <- PlotRectHeatmapWithSignif(SubClonal.tnb.diff.nongenomic, x = "CancerType", y = "EscapeMetric", statistic.p = "p.value", tag = "TNB.subclone_median", right_plot_type = "proportion", right_plot_xlim = 0.5)

plot.list <- patchwork::wrap_plots(p1, p2, ncol = 1)

ggsave(file.path(out.dir, "TNB_clonality.non_genomic_escape.pdf"), plot.list, width = 15, height = 20)
