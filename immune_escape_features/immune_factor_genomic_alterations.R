# Identify IFNG-pathway, IDH1, and CD58 genomic alterations from mutation and copy-number inputs.
library(dplyr)
source(file.path("functions", "construct_gene_level_copy_number_matrix.R"))
source(file.path("functions", "construct_gene_level_mutation_matrix.R"))
ImmuneCytokineGenomic <- function(mut.data,
                                  CN.data = NULL,
                                  CN.Matrix.data = NULL,
                                  RefGenome = c("hg19", "hg38")
) {
    # Detect homozygous deletions from raw segments or a discrete gene-level copy-number matrix.
    if (!is.null(CN.data)) {
        tmp.CN.data <- CN.data[which(CN.data$chromosome %in% c("chr1", "chr2", "chr4", "chr6", "chr9", "chr11", "chr21")), ]
        Allgene.HD <- ConstructCNVMatrix(CN.data = tmp.CN.data, RefGenome = RefGenome, CNV.types = "HD", type.column = "CNProfile.type")

        IFNG.gene <- intersect(c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1"), rownames(Allgene.HD))
        df.IFNG <- Allgene.HD[intersect(IFNG.gene, rownames(Allgene.HD)), , drop = FALSE]
        IFNG.HD <- data.frame(SampleID = colnames(Allgene.HD), IFNG.pathway.HD = colSums(df.IFNG) != 0)

        df.CD58 <- Allgene.HD[intersect("CD58", rownames(Allgene.HD)), , drop = FALSE]
        CD58.HD <- data.frame(SampleID = colnames(Allgene.HD), CD58.HD = colSums(df.CD58) != 0)

        factor.HD <- dplyr::full_join(IFNG.HD, CD58.HD, by = "SampleID")

    } else if (!is.null(CN.Matrix.data)) {
        IFNG.pathway.HD <- CN.Matrix.data[intersect(rownames(CN.Matrix.data), c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")), , drop = FALSE]
        CD58.HD <- CN.Matrix.data[intersect(rownames(CN.Matrix.data), "CD58"), , drop = FALSE]
        factor.HD <- data.frame(
            SampleID = colnames(CN.Matrix.data),
            IFNG.pathway.HD = colSums(IFNG.pathway.HD == -2) > 0,
            CD58.HD = colSums(CD58.HD == -2) > 0
        )
    } else {
        factor.HD <- NULL
    }
    # Use nonsilent mutations and retain pathway-level and individual-gene status.
    if (!is.null(mut.data)) {
        gene.mut <- ConstructMutationMatrix(mutation.data = mut.data, mut.types = "nonsilent")

        IFNG.gene <- intersect(c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1"), rownames(gene.mut))
        IFNG.pathway.mut <- data.frame(
            SampleID = colnames(gene.mut),
            IFNG.pathway.mut = ifelse(colSums(gene.mut[intersect(IFNG.gene, rownames(gene.mut)), , drop = FALSE]) == 0, FALSE, TRUE)
        )

        IFNG.gene <- c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")
        for (g in IFNG.gene) {
            if (g %in% rownames(gene.mut)) {
             IFNG.pathway.mut[[paste0(g, ".mut")]] <- as.logical(gene.mut[g, ])
            } else {
            IFNG.pathway.mut[[paste0(g, ".mut")]] <- FALSE
            }
        }

        if ("IDH1" %in% rownames(gene.mut)) {
            IDH1.mut <- data.frame(SampleID = colnames(gene.mut), IDH1.mut = gene.mut["IDH1", ] != 0)
        } else {
            IDH1.mut <- data.frame(SampleID = colnames(gene.mut), IDH1.mut = FALSE)
        }

        if ("CD58" %in% rownames(gene.mut)) {
            CD58.mut <- data.frame(SampleID = colnames(gene.mut), CD58.mut = gene.mut["CD58", ] != 0)
        } else {
            CD58.mut <- data.frame(SampleID = colnames(gene.mut), CD58.mut = FALSE)
        }

        factor.mut <- Reduce(
            function(x, y) merge(x, y, by = "SampleID", all = TRUE),
            list(IFNG.pathway.mut, IDH1.mut, CD58.mut)
        )
    } else {
        factor.mut <- NULL
    }

result <- Reduce(
        function(x, y) merge(x, y, by = "SampleID", all = TRUE, sort = FALSE),
        Filter(Negate(is.null), list(factor.HD, factor.mut))
    )
    return(result)
}
