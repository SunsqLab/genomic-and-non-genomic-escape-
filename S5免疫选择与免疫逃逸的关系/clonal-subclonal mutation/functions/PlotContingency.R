# Original mutation-event contingency plot; the main script explicitly requests p_display = 'raw'.
plot_contingency_from_df <- function(data,
                                     x_col,
                                     y_col,
                                     title = "Distribution Comparison",
                                     p_display = c("formatted", "raw", NULL),
                                     colors = c("#AABCDB", "#C0D6EA")) {

    require(ggplot2)
    require(dplyr)
    require(scales)

    if (!x_col %in% names(data)) stop(paste("Column", x_col, "not found"))
    if (!y_col %in% names(data)) stop(paste("Column", y_col, "not found"))

    if (length(unique(data[[x_col]])) != 2) stop(paste("Column", x_col, "must have exactly 2 unique values"))
    if (length(unique(data[[y_col]])) != 2) stop(paste("Column", y_col, "must have exactly 2 unique values"))

    contingency_table <- table(data[[y_col]], data[[x_col]])

    data.plot <- as.data.frame(contingency_table) %>%
        rename(y_category = Var1, x_category = Var2, Freq = Freq) %>%
        group_by(x_category) %>%
        mutate(Proportion = Freq / sum(Freq)) %>%
        ungroup()

    data.plot$Label <- paste0(
        format(data.plot$Freq, big.mark = ","),
        "\n(", round(data.plot$Proportion * 100), "%)"
    )

    p_label <- NULL

    if (!is.null(p_display)) {
        p_value <- fisher.test(contingency_table)$p.value

        if (p_display == "formatted") {

            p_label <- dplyr::case_when(
                p_value < 0.001 ~ "p < 0.001",
                p_value < 0.01 ~ "p < 0.01",
                p_value < 0.05 ~ "p < 0.05",
                TRUE ~ paste("p =", format(round(p_value, 3), nsmall = 3))
            )
        } else if (p_display == "raw") {
            p_label <- paste("p =", format(p_value, scientific = TRUE, digits = 3))
        }
    }

    p <- ggplot(data.plot, aes(x = x_category, y = Proportion, fill = y_category)) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4, fontface = "bold"
        ) +
        scale_fill_manual(values = colors) +
        scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
        labs(title = title, x = x_col, y = "Proportion", fill = y_col) +
        theme_classic() +

        annotate("text", x = 1.5, y = 1.05, label = p_label, size = 4)

    return(p)
}
