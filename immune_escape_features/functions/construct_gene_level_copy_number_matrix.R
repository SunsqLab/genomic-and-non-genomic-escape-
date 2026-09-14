# Build a binary gene-by-sample matrix for selected copy-number event classes.
ConstructCNVMatrix <- function(CN.data,
                               type.column = c("type", "CNProfile.type")[2],
                               CNV.types,
                               RefGenome = c("hg19", "hg38")) {
    library(reshape2)
    source(file.path("functions", "annotate_genes_in_copy_number.R"))

    # Assign genes only when their full span is covered by a copy-number segment.
    gene.cnv <- Annotate.genes.in.copynumber(CNV.infor = CN.data, ref.genome = RefGenome, gene.coverage.ratio = 100)
    gene.cnv$index <- paste(gene.cnv$chr, gene.cnv$CNV.start, gene.cnv$CNV.end, sep = "_")

    CN.data$index <- paste(CN.data$chromosome, CN.data$startPosition, CN.data$endPosition, sep = "_")
    library(data.table)
    setDT(CN.data)
    merged_data <- merge(
        CN.data,
        gene.cnv[, c("index", "gene_name")],
        by = "index",
        all.x = TRUE,
        allow.cartesian = TRUE
    )
    tmp.gene.cnv <- merged_data[!is.na(merged_data$gene_name), ]

    tmp <- dcast(tmp.gene.cnv, gene_name ~ SampleID, value.var = type.column)
    tmp <- as.data.frame(tmp)

    rownames(tmp) <- tmp$gene_name
    tmp <- tmp[, -1]

    result <- matrix(0, nrow = nrow(tmp), ncol = ncol(tmp), dimnames = list(rownames(tmp), colnames(tmp)))
    result[which(unlist(tmp) %in% CNV.types)] <- 1

    return(result)
}
