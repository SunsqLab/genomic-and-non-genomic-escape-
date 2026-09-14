# Construct unbinarized immune escape feature tables for downstream binarization.
source("epigenetic_escape_features.R")
source("antigen_presentation_features.R")
source("immune_checkpoint_features.R")
source("immune_microenvironment_features.R")
source("neoantigen_expression_features.R")
EscapeFeatures <- function(exprs.data = NULL,
                           mut.data = NULL,
                           CN.data = NULL,
                           CN.Matrix.data = NULL,
                           RefGenome = c("hg19", "hg38")
) {
    exprs.data <- na.omit(exprs.data)

    APM.feature <- AntigenPresentDefect(exprs.data = exprs.data, mut.data = mut.data, CN.data = CN.data, CN.Matrix.data = CN.Matrix.data, RefGenome = RefGenome)

    # Separate expression-derived antigen-presentation features from genomic alteration features.
    t.index <- intersect(colnames(APM.feature), c("SampleID", "CALR", "HLA.A", "HLA.B", "HLA.C", "HLA.score"))

    if (length(t.index) != 1) {
        HLA.exp <- APM.feature[, t.index, drop = FALSE]
    } else {
        HLA.exp <- NULL
    }

    APM.alt <- APM.feature[, setdiff(colnames(APM.feature), c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR")), drop = FALSE]

    Checkpoint.feature <- ImmuneCheckpoint(exprs.data = exprs.data, CN.data = CN.data, CN.Matrix.data = CN.Matrix.data, RefGenome = RefGenome)

    ImmuneSuppre.feature <- ImmuneSuppreEnvir(exprs.data = exprs.data, mut.data = mut.data, CN.data = CN.data, CN.Matrix.data = CN.Matrix.data, RefGenome = RefGenome)

    Suppre.Cell <- ImmuneSuppre.feature$Suppre.Cell

    Immune.Factor <- ImmuneSuppre.feature$Immune.Factor

    Factor.Genomic <- ImmuneSuppre.feature$Factor.Genomic

    NeoantigenExp.feature <- NeoantigenExprs(mut.data = mut.data)

    Epigenetic.feature <- EpigeneticAlt(CN.data = CN.data, RefGenome = RefGenome)

result <- list(
        APM.alt = APM.alt, Checkpoint.feature = Checkpoint.feature, Suppre.Cell = Suppre.Cell, Immune.Factor = Immune.Factor, Factor.Genomic = Factor.Genomic,
        HLA.exp = HLA.exp, NeoantigenExp.feature = NeoantigenExp.feature, Epigenetic.feature = Epigenetic.feature
    )
    return(result)
}
