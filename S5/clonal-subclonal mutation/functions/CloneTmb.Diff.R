# Original log10 median-ratio effect, two-sided Wilcoxon tests and within-cancer BH adjustment.
CloneTmb.Diff <- function(data, feature, Cancer, Target) {
    library(dplyr)
    library(rlang)

    data.plot <- tidyr::pivot_longer(data, cols = tidyr::all_of(unname(feature)), names_to = "EscapeMetric", values_to = "EscapeStatus") %>%
        data.frame()

    res <- data.plot %>%

        dplyr::select(!!sym(Cancer), EscapeMetric, EscapeStatus, !!sym(Target)) %>%
        dplyr::filter(complete.cases(.)) %>%
        dplyr::group_by(!!sym(Cancer), EscapeMetric) %>%

        dplyr::filter(n_distinct(EscapeStatus) == 2 & all(table(EscapeStatus) >= 5)) %>%

        dplyr::summarise(
            fold.change = round(log10(median(.data[[Target]][EscapeStatus == TRUE]) + exp(-6)) - log10(median(.data[[Target]][EscapeStatus == FALSE]) + exp(-6)), 3),
            p.value = wilcox.test(.data[[Target]] ~ EscapeStatus, exact = FALSE)$p.value,
            .groups = "drop"
        ) %>%

        dplyr::group_by(!!sym(Cancer)) %>%
        dplyr::mutate(
            p.value = p.adjust(p.value, "BH"),
        ) %>%
        dplyr::ungroup() %>%

        dplyr::group_by(!!sym(Cancer)) %>%
        dplyr::mutate(
            significance = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "")),
            EscapeMetric = factor(EscapeMetric, levels = rev(feature)),

            category = if (is.null(names(feature))) NA else factor(names(feature)[match(as.character(EscapeMetric), unname(feature))], levels = unique(names(feature)))
        ) %>%
        dplyr::ungroup() %>%
        data.frame(.)

    return(res)
}
