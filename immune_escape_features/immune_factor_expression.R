# Extract individual immune-factor expression values and an ssGSEA cytokine-signature score.
ImmuneSuppreCytokine <- function(exprs.data) {
    library(dplyr)

    tmp.TGFB <- as.numeric(exprs.data["TGFB1", ])

# Score the predefined immunosuppressive cytokine signature by ssGSEA.
    Immunosuppressive_cytokine <- c("TGFB1", "IL1RN", "GSF3", "IL23A", "PTGS2", "IL10", "IL13", "CSF1", "MIF", "IL6", "INHBA", "IL4", "IDO1", "HDC", "LGALS3", "IL5", "CSF3", "THBS1", "KLRF1", "SH2D1B", "FYN", "LCK", "SH3BP2", "GRB2", "PTPRC", "SH2D1A", "LCP2", "ZAP70", "KLRC3", "KLRD1", "HLA-E", "TYROBP", "KIR2DS2", "KIR3DS1", "NCR2", "KLRC2", "HLA-C", "PCNA", "KMT2E")

library(GSVA)
    tmp.pheno <- gsva(expr = data.matrix(exprs.data), gset.idx.list = list(Immunosuppressive_cytokine), method = "ssgsea")

    result <- data.frame(SampleID = colnames(exprs.data), TGFB1 = tmp.TGFB, Immune_supress_cytokine = as.numeric(tmp.pheno))
    # Keep predefined suppressive and immune-activating genes in a fixed output order.
    Suppressor_gene <- c("SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2")
    Active_gene <- c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1")
    genes <- c(Suppressor_gene, Active_gene)

sup.factor <- exprs.data[match(genes,rownames(exprs.data)), , drop=TRUE]
    sup.factor <- data.frame(SampleID = colnames(exprs.data), t(sup.factor))
    colnames(sup.factor) <- c("SampleID", genes)

    result <- merge(result, sup.factor, by="SampleID")

return(result)
}
