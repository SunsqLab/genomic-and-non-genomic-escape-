# Annotate genes fully covered by copy-number segments.
Annotate.genes.in.copynumber <- function(CNV.infor, ref.genome = "hg38", gene.coverage.ratio = 100) {
    library(IRanges)
    library(GenomicRanges)
    source(file.path("functions", "overlap_ranges.R"))

    query <- GRanges(
        seqnames = CNV.infor[, "chromosome"],
        ranges = paste0(CNV.infor[, "startPosition"], "-", CNV.infor[, "endPosition"])
    )

    if (ref.genome == "hg19") {
        gene.infor <- unique(get(load(file.path("reference_data", "gene_coordinates_hg19.RData")))[, c("seqnames", "start", "end", "gene_name")])
        subject <- GRanges(
            seqnames = gene.infor[, "seqnames"],
            ranges = paste0(gene.infor[, "start"], "-", gene.infor[, "end"]),
            state = gene.infor[, "gene_name"]
        )
    }
    if (ref.genome == "hg38") {
        gene.infor <- unique(get(load(file.path("reference_data", "gene_coordinates_hg38.RData")))[, c("seqnames", "start", "end", "gene_name")])
        subject <- GRanges(
            seqnames = gene.infor[, "seqnames"],
            ranges = paste0(gene.infor[, "start"], "-", gene.infor[, "end"]),
            state = gene.infor[, "gene_name"]
        )
    }
    overlap.information <- olRanges(query, subject, output = "df")

    if (is.na(overlap.information[1, 1])) {
        anno.CNV.gene <- NULL
    }

    if (!is.na(overlap.information[1, 1])) {
        # Retain genes meeting the requested percentage of subject-gene coverage.
        overlap.CNV.gene <- overlap.information[which(overlap.information[, "OLpercS"] >= gene.coverage.ratio), ]

        overlap.CNV.gene[, dim(overlap.CNV.gene)[2] + 1] <- paste0(
            overlap.CNV.gene[, "space"], "-", overlap.CNV.gene[, "Sstart"], "-", overlap.CNV.gene[, "Send"]
        )
        colnames(overlap.CNV.gene)[dim(overlap.CNV.gene)[2]] <- "overlap.index"
        gene.infor[, dim(gene.infor)[2] + 1] <- paste0(
            gene.infor[, "seqnames"], "-", gene.infor[, "start"], "-", gene.infor[, "end"]
        )
        colnames(gene.infor)[dim(gene.infor)[2]] <- "gene.index"

        library(data.table)
        setDT(overlap.CNV.gene)
        setDT(gene.infor)
        anno.CNV.gene <- merge(
            overlap.CNV.gene,
            gene.infor,
            by.x = "overlap.index",
            by.y = "gene.index",
            sort = FALSE,
            allow.cartesian = TRUE
        )

        anno.CNV.gene <- anno.CNV.gene[, c("space", "Qstart", "Qend", "Sstart", "Send", "gene_name")]
        colnames(anno.CNV.gene) <- c("chr", "CNV.start", "CNV.end", "gene.start", "gene.end", "gene_name")
        anno.CNV.gene <- unique(anno.CNV.gene)
    }

    return(anno.CNV.gene)
}
