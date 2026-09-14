# Definitions of the 50 immune escape features used in the final panel.
features <- list(
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
