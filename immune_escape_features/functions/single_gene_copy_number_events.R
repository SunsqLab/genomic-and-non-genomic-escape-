# Extract selected copy-number event classes into a sample-by-gene table.
source(file.path("functions", "construct_gene_level_copy_number_matrix.R"))
source(file.path("functions", "extract_gene_coordinates.R"))

CNVGeneAmpDel <- function(CN.data,
    genes,
    RefGenome = "hg38",
    type.column = c("type", "CNProfile.type")[2],
    CNV.types = c("High-levelAmplication")
) {
    t.genes <- extract.coor.gene(geneSet = genes, genome.build = RefGenome)
    tmp.CN.data <- CN.data[which(CN.data$chromosome %in% as.character(t.genes$seqnames)), ]

    # Construct the selected gene-level copy-number event matrix.
    Allgene.deepAmp <- ConstructCNVMatrix(
        CN.data = tmp.CN.data,
        RefGenome = RefGenome,
        type.column = type.column,
        CNV.types = CNV.types
    )
    if (any(genes %in% rownames(Allgene.deepAmp))) {
        share.genes <- intersect(genes, rownames(Allgene.deepAmp))
        tmp <- Allgene.deepAmp[share.genes, ]
        tmp <- matrix(tmp, nrow = ncol(Allgene.deepAmp))

        genes.amp <- data.frame(colnames(Allgene.deepAmp), ifelse(tmp == 0, FALSE, TRUE))
        colnames(genes.amp) <- c("SampleID", paste0(share.genes, ".deepAmp"))
    } else {
        genes.amp <- NULL
    }
    return(genes.amp)
}
