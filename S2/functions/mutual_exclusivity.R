# Original pairwise contingency counts, raw odds ratios and one-sided Fisher tests.
odds_ratio_matrix_fun <- function(mat, alternative = c("less", "greater"), cut_off = 0.1 , p = 0.05, fdr = 0.05, verbose = T) {

    all_choose <- t(combn(rownames(mat), m = 2))
    if (verbose) message("=> There are ", nrow(all_choose), " candidate pairs")

    if (verbose) message("==> Calculating odds ratio ...")
    res_list <- apply(all_choose, 1, function(x) {
        tables <- table(mat[x[1],], mat[x[2],])
        odds_ratio <- (tables[2,2] * tables[1,1])/(tables[2,1] * tables[1,2])
        p_value <- fisher.test(matrix(c(tables[2,2], tables[1,2], tables[2,1], tables[1,1]), ncol = 2), alternative = alternative)$p.value
        detail_res <- c(tables[1,1], tables[2,1], tables[1,2], tables[2,2], odds_ratio, p_value)
        detail_res
    })

    res_list <- as.data.frame(t(res_list))
    colnames(res_list) <- c("Neither", "A_not_B", "B_not_A", "Both", "odds_ratio", "p_value")
    res_list$q_value <- p.adjust(res_list$p_value, method = "BH")
    res_list <- cbind(all_choose, res_list)
    colnames(res_list)[1:2] <- c("geneA", "geneB")
    if (alternative == "less") {
        if (verbose) message("==> Choosing mutually exclusive pairs ...")
        res_list <- dplyr::filter(res_list, p_value < p & q_value < fdr)
        res_list <- dplyr::filter(res_list, odds_ratio < cut_off)
        if (verbose) message("==> Get ", nrow(res_list)," mutually exclusive pairs ...")
    } else {
        if (verbose) message("==> Choosing co-occurring pairs ...")
        res_list <- dplyr::filter(res_list, p_value < p & q_value < fdr)
        res_list <- dplyr::filter(res_list, odds_ratio > cut_off)
        if (verbose) message("==> Get ", nrow(res_list)," co-occurring pairs ...")
    }
    if (verbose) message("=> Done")
    return(res_list)
}
