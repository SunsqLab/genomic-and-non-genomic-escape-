# Original gene-by-sample mutation-matrix function.
source(file.path(functions.dir, "MutationType.R"))
source(file.path(functions.dir, "FilterGenomicProfiler.R"))
ConstructMutationMatrix <- function(mutation.data,
                                    mut.types = "nonsilent",
                                    min.num.frequent = 1) {
    library(dplyr)
    library(tidyr)
    library(tibble)
    mut.types <- MutationType(type.mut = mut.types)

    mutation.data <- mutation.data[mutation.data$geneSymbol != "" & !is.na(mutation.data$geneSymbol), ]

    mut.mat <- mutation.data %>%
        dplyr::select(SampleID, geneSymbol, Level5) %>%
        dplyr::group_by(SampleID, geneSymbol) %>%

        dplyr::summarize(mutations = as.numeric(any(Level5 %in% mut.types)), .groups = "drop") %>%

        spread(SampleID, mutations) %>%
        column_to_rownames("geneSymbol") %>%
        as.matrix()
    mut.mat[is.na(mut.mat)] <- 0
    mut.mat <- mut.mat[rownames(mut.mat) != "", ]

    mut.mat <- FilterGenomicProfiler(genomic.profiler = mut.mat, min.num.frequent = min.num.frequent)

    return(mut.mat)
}
