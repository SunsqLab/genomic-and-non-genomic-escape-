# Summarize neoantigen-associated expression metrics from mutation-level data.
NeoantigenExprs <- function(mut.data
) {
    library(magrittr)

    if (!is.null(mut.data) & all(c("SampleID", "geneExp", "isNeoantigen", "isClonal") %in% colnames(mut.data))) {
        tmp.data <- mut.data[, c("SampleID", "geneExp", "isNeoantigen", "isClonal")]

        # Calculate mean expression and the proportion of expressed neoantigens for each sample.
        tmp <- tmp.data %>%
            dplyr::group_by(SampleID) %>%
            dplyr::reframe(
                mean.neo.exprs = mean(geneExp[which(isNeoantigen)], na.rm = T),
                mean.clonal.neo.exprs = mean(geneExp[which(isNeoantigen & isClonal)], na.rm = T),
                neo.exprs.ratio = sum(isNeoantigen & (geneExp > log2(1 + 1)), na.rm = T) / sum(isNeoantigen, na.rm = T),
                clonal.neo.exprs.ratio = sum(isNeoantigen & isClonal & (geneExp > log2(1 + 1)), na.rm = T) / sum(isNeoantigen & isClonal, na.rm = T),
                subclonal.neo.exprs.ratio = sum(isNeoantigen & (!isClonal) & (geneExp > log2(1 + 1)), na.rm = T) / sum(isNeoantigen & (!isClonal), na.rm = T)
            ) %>%
            as.data.frame()
        tmp[which(tmp == "NaN", arr.ind = T)] <- NA
    } else {
        tmp <- NULL
    }

    return(tmp)
}
