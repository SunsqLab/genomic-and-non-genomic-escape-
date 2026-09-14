# S6: ITH by escape status
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

ID.xteam <- "PanCancer_TCGA.dataset"
source("/pub5/xiaoyun/BioX/Bioc/0.BioData/PanCancer/DataCenter/[[[Navigation]]]/[MultiDatasetProcessing]/0.integrate_dataset_labels.R")
ID.xteams <- get(ID.xteam)
is.pass = TRUE
evasion.list.all = list(
    APM.alt = c("HLA.LOH", "B2M.BiallelicInactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),
    Checkalt = c("CD274.deepAmp"),
    ActMalt = c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
    Check.exp = c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    SupprCell = c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    SupprSig = c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    ActGenedown = c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown = c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo = c("mean.neo.exprs")
)
query.col = c("prop.subclonal.mut", "MATH")
output.dir = "/WorkSpace/heshengyuan/01ImmuneEditing/10ImmuneEscape/50.ImmuneEscape/[Q]BaselineEffectsOfImmuneEscape/ImmuneEscapeAndITH/Revision"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
output.files = unlist(lapply(query.col, function(col){
    c(file.path(output.dir, paste(col, c("dot_heatmap_[log2FoldChangesignificance].pdf", "dot_heatmap_[significance].pdf", "metric_association_heatmap_[FCsignificance].pdf", "tile_heatmap_[significance].pdf"), sep = "_")),
        file.path(output.dir, paste(col, "metric_association_heatmap_[FCsignificance].xlsx", sep = "_")))
}))

if (is.pass & all(file.exists(output.files))) {
    cat(blue("All result files already exist"))
} else {
    source("/pub5/xiaoyun/BioX/Bioc/0.BioData/PanCancer/DataCenter/[BioImmune]/[[[Questions]]]/[Q]ImmuneEscapeProfile/immune_escape_metrics_multi.R")
    for(col in query.col){
        result = ImmuneEscape_ColRelationX(
            ID.xteams,
            query.col = col,
            obs.cols = NULL,
            obs.cols.group = evasion.list.all,
            output.dir = output.dir
        )
    }

}
