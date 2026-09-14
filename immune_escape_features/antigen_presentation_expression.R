# Extract HLA class I and CALR expression from a gene-expression matrix.
AntigenPresentDefectHLAexpr <- function(exprs.data
) {
    exp.data <- t(exprs.data[intersect(rownames(exprs.data), c("HLA-A", "HLA-B", "HLA-C")), , drop = FALSE])

    # Calculate the HLA score when HLA-A, HLA-B, and HLA-C are all available.
    if (ncol(exp.data) == 3) {
        HLA.score <- apply(exp.data, 1, function(x) {
            tmp <- exp(mean(log(x)))
        })
        result <- data.frame(rownames(exp.data), exp.data, HLA.score)
        colnames(result)[1] <- c("SampleID")
    } else if (ncol(exp.data) >= 1) {
        result <- data.frame(rownames(exp.data), exp.data)
        colnames(result)[1] <- c("SampleID")
    } else {
        result <- NULL
    }
    # Add CALR expression when it is available in the input matrix.
    calr.exp <- t(exprs.data[intersect(rownames(exprs.data), c("CALR")), , drop = FALSE])

    if (ncol(calr.exp) == 1) {
        calr.exp <- data.frame(rownames(calr.exp), calr.exp)
        colnames(calr.exp)[1] <- "SampleID"
        if (is.null(result)) {
            result <- calr.exp
        } else {
            result <- merge(result, calr.exp, by = "SampleID")
        }
    }

return(result)
}
