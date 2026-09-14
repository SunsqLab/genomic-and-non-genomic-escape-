# Original heatmap and side bars; only the optional effect-size label is new.
library(ggplot2)
library(patchwork)
library(rlang)
library(dplyr)

PlotRectHeatmapWithSignif <- function(plot.data = NULL,
                                      x,
                                      y,
                                      statistic.type = "fold.change",
                                      statistic.p = "p.value",
                                      Metric.category = "category",
                                      tag = "Clonal.TMB",
                                      right_plot_type = "count",
                                      right_plot_xlim = NULL,
                                      size_title = "|log10FC|") {

    no_facet <- is.na(Metric.category)

    data_plot <- plot.data %>%
        dplyr::mutate(
            significance = cut(
                as.numeric(.data[[statistic.p]]),
                breaks = c(0, 0.001, 0.01, 0.05, 1),
                include.lowest = TRUE,
                right = FALSE,
                labels = c("***", "**", "*", "")
            )
        )
    size_values <- abs(data_plot[[statistic.type]])
    size_values <- size_values[is.finite(size_values)]
    if (length(size_values) == 0) {
        size_breaks <- waiver()
    } else if (length(unique(size_values)) == 1) {
        size_breaks <- unique(size_values)
    } else {
        size_range <- range(size_values)
        size_breaks <- pretty(size_range, n = 4)
        size_breaks <- size_breaks[size_breaks >= size_range[1] & size_breaks <= size_range[2]]
        if (length(size_breaks) == 0) {
            size_breaks <- size_range
        }
    }

    p_main <- ggplot(data_plot, aes(x = !!sym(x), y = !!sym(y))) +

        geom_hline(
            yintercept = 1:(length(unique(data_plot[[y]])) - 1) + 0.5,
            color = "grey90", linewidth = 0.5
        ) +
        geom_vline(
            xintercept = 1:(length(unique(data_plot[[x]])) - 1) + 0.5,
            color = "grey90", linewidth = 0.5
        ) +
        geom_point(aes(
            size = abs(!!sym(statistic.type)),
            color = factor(!!sym(statistic.type) > 0),
            alpha = !is.na(significance) & significance != ""
        ), shape = 15) +
        scale_alpha_manual(values = c(0.2, 1), guide = "none") +
        scale_color_manual(
            values = c("#009ccc", "#fe0000"),
            labels = c("Negative", "Positive"),
            name = "Direction",
            guide = guide_legend(order = 1)
        ) +
        scale_size_continuous(
            range = c(1, 5),
            breaks = size_breaks,
            name = size_title,
            guide = guide_legend(
                direction = "horizontal",
                title.position = "top",
                order = 2,
                nrow = 1,
                byrow = TRUE,
                override.aes = list(color = "black")
            )
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 20)),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "grey"),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.text.y = element_text(size = 11),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.spacing.x = unit(1, "cm"),
            legend.margin = margin(t = 10),
            strip.text.y = element_text(size = 11),
            strip.background = element_rect(fill = "grey92", color = "grey50", linewidth = 0.5),
            panel.spacing.y = unit(0.5, "lines")
        ) +
        labs(x = NULL, y = NULL, title = paste("Comparison between", tag))

    if (!no_facet) {
        p_main <- p_main +
            facet_grid(
                rows = vars(!!sym(Metric.category)),
                scales = "free_y",
                space = "free_y",
                switch = "y"
            )
    }

    # Original denominator: displayed cohorts, including PanCancer if present.
    total_cancer_types <- dplyr::n_distinct(data_plot[[x]])

    if (no_facet) {
        freq_data <- data_plot %>%
            dplyr::group_by(!!sym(y))
    } else {
        freq_data <- data_plot %>%
            dplyr::group_by(!!sym(Metric.category), !!sym(y))
    }

    if (right_plot_type == "proportion") {
        freq_data <- freq_data %>%
            dplyr::summarise(
                high = sum(!!sym(statistic.type) > 0 & significance != "", na.rm = TRUE) / total_cancer_types,
                low = sum(!!sym(statistic.type) < 0 & significance != "", na.rm = TRUE) / total_cancer_types,
                .groups = "drop"
            )
    } else if (right_plot_type == "count") {
        freq_data <- freq_data %>%
            dplyr::summarise(
                high = sum(!!sym(statistic.type) > 0 & significance != "", na.rm = TRUE),
                low = sum(!!sym(statistic.type) < 0 & significance != "", na.rm = TRUE),
                .groups = "drop"
            )
    } else {
        stop("right_plot_type must be either 'proportion' or 'count'.")
    }

    freq_data <- freq_data %>%
        tidyr::pivot_longer(
            cols = c(high, low),
            names_to = "direction",
            values_to = "value"
        )

    p_right <- ggplot(freq_data, aes(x = value, y = !!sym(y), fill = direction)) +
        geom_col(position = "dodge", width = 0.7) +
        scale_fill_manual(
            values = c("high" = "#fe0000", "low" = "#009ccc"),
            labels = c("high" = "Sig-Up", "low" = "Sig-Down"),
            name = "Significance",
            guide = guide_legend(order = 3)
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            panel.grid.major.y = element_blank()
        )

    if (right_plot_type == "proportion") {
        xlim_upper <- if (!is.null(right_plot_xlim)) right_plot_xlim else 1.0
        p_right <- p_right +
            scale_x_continuous(
                labels = scales::percent_format(accuracy = 1),
                expand = c(0, 0),
                limits = c(0, xlim_upper)
            )
    } else {
        p_right <- p_right +
            scale_x_continuous(expand = c(0, 0))
    }

    if (!no_facet) {
        p_right <- p_right +
            facet_grid(
                rows = vars(!!sym(Metric.category)),
                scales = "free_y",
                space = "free_y"
            )
    }

    combined_plot <- p_main + p_right +
        plot_layout(widths = c(4, 1))

    return(combined_plot)
}
