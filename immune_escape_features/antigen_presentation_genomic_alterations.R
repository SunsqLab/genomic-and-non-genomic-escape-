# Identify antigen-presentation defects from mutation and copy-number data.
library(dplyr)
source(file.path("functions", "construct_gene_level_copy_number_matrix.R"))
source(file.path("functions", "construct_gene_level_mutation_matrix.R"))
AntigenPresentDefectAlteration <- function(mut.data,
                                           CN.data = NULL,
                                           CN.Matrix.data = NULL,
                                           RefGenome = c("hg19", "hg38")
) {

    # Use CNProfile.type to identify LOH-related copy-number events.
    if (!is.null(CN.data)) {
        tmp.CN.data <- CN.data[which(CN.data$chromosome %in% c("chr6", "chr15")), ]
        Allgene.LOH <- ConstructCNVMatrix(CN.data = tmp.CN.data, RefGenome = RefGenome, type.column = "CNProfile.type", CNV.types = c("CNLOH", "HD", "LOH", "ModerateGain&LOH", "MinorGain&LOH", "High-levelAmplication&LOH"))

        # Determine HLA class I LOH across HLA-A, HLA-B, and HLA-C.
        HLAgene.LOH <- data.frame(SampleID = colnames(Allgene.LOH), HLA.LOH = ifelse(colSums(Allgene.LOH[intersect(c("HLA-A", "HLA-B", "HLA-C"), rownames(Allgene.LOH)), , drop = FALSE]) == 0, FALSE, TRUE))

        if ("B2M" %in% rownames(Allgene.LOH)) {
            B2M.LOH <- data.frame(SampleID = colnames(Allgene.LOH), B2M.LOH = ifelse(Allgene.LOH[c("B2M"), ] == 0, FALSE, TRUE))
        } else {
            B2M.LOH <- data.frame(SampleID = colnames(Allgene.LOH), B2M.LOH = FALSE)
        }

        tmp.CN.data <- CN.data[which(CN.data$chromosome %in% c("chr1", "chr6", "chr16", "chr19")), ]
        Allgene.HD <- ConstructCNVMatrix(CN.data = tmp.CN.data, RefGenome = RefGenome, CNV.types = "HD", type.column = "CNProfile.type")

        genes <- c("CALR", "TAP1", "TAP2", "TAPBP", "NLRC5", "CIITA", "RFX5")
        APM.pathway.HD <- data.frame(
            SampleID = colnames(Allgene.HD),
            sapply(setNames(genes, paste0(genes, ".HD")), function(g) {
                if (g %in% rownames(Allgene.HD)) Allgene.HD[g, ] != 0 else rep(FALSE, ncol(Allgene.HD))
            })
        )
    } else if (!is.null(CN.Matrix.data)) {
        tmp <- CN.Matrix.data[intersect(rownames(CN.Matrix.data), c("HLA-A", "HLA-B", "HLA-C")), , drop = FALSE]
        HLAgene.LOH <- data.frame(SampleID = colnames(CN.Matrix.data), HLA.LOH = colSums(tmp == -1 | tmp == -2) > 0)
        tmp <- CN.Matrix.data[intersect(rownames(CN.Matrix.data), "B2M"), , drop = FALSE]
        B2M.LOH <- data.frame(SampleID = colnames(CN.Matrix.data), B2M.LOH = colSums(tmp == -1 | tmp == -2) > 0)

        genes <- c("CALR", "TAP1", "TAP2", "TAPBP", "NLRC5", "CIITA", "RFX5")
        APM.pathway.HD <- data.frame(
            SampleID = colnames(CN.Matrix.data),
            sapply(setNames(genes, paste0(genes, ".HD")), function(g) {
                if (g %in% rownames(CN.Matrix.data)) CN.Matrix.data[g, ] == -2 else rep(FALSE, ncol(CN.Matrix.data))
            }),
            check.names = FALSE
        )
    } else {
        HLAgene.LOH <- NULL
        B2M.LOH <- NULL
        APM.pathway.HD <- NULL
    }

    # Use nonsilent mutations consistently for all antigen-presentation genes.
    if (!is.null(mut.data)) {
        gene.mut <- ConstructMutationMatrix(mutation.data = mut.data, mut.types = "nonsilent")

        all.APM.gene <- c(
            "HLA-A", "HLA-B", "HLA-C",
            "B2M", "CALR", "TAP1", "TAP2", "TAPBP", "CIITA", "RFX5", "NLRC5",
            "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
            "CANX", "CD4", "CD74", "CD8A", "CD8B", "CREB1", "CTSB", "CTSL", "CTSS", "ERAP1", "ERAP2", "FAS", "HLA-E", "HLA-F", "HLA-G", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "HSPBP1", "IFI30", "IFNG", "IRF1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "KLRD1", "LGMN", "MEX3B", "NFYA", "NFYB", "NFYC", "PDIA3", "PSMA7", "PSMB10", "PSMB11", "PSMB6", "PSMB8", "PSMB9", "PSME1", "PSME2", "PSME3", "PSMF1", "RFXANK", "RFXAP", "TNF"
        )

        inter.APM <- intersect(rownames(gene.mut), all.APM.gene)
        tmp <- matrix(ifelse(t(gene.mut[inter.APM, , drop = FALSE]) == 0, FALSE, TRUE), ncol = length(inter.APM))
        if (ncol(tmp) > 0) {
            APM.mut <- data.frame(SampleID = colnames(gene.mut), tmp)
        } else {
            APM.mut <- data.frame(SampleID = colnames(gene.mut))
        }
        colnames(APM.mut) <- c("SampleID", inter.APM)
        APM.mut[, setdiff(c(all.APM.gene), inter.APM)] <- FALSE

        colnames(APM.mut) <- c("SampleID", paste0(tail(colnames(APM.mut), -1), ".mut"))

        APM.mut$HLA.mut <- rowSums(APM.mut[, c("HLA-A.mut", "HLA-B.mut", "HLA-C.mut")], na.rm = TRUE)
        APM.mut$all.APM.mut <- rowSums(APM.mut[, paste0(all.APM.gene, ".mut")], na.rm = TRUE)
    } else {
        APM.mut <- NULL
    }

    tmp.result <- list(HLAgene.LOH, B2M.LOH, APM.pathway.HD, APM.mut)
    result <- Reduce(dplyr::full_join, tmp.result[lengths(tmp.result) != 0])
    return(result)
}
