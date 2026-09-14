# Original mirrored mutation counts with a two-sided paired Wilcoxon annotation.
mirror_bar <- function(data, x, y.up, y.down,
                       paired = TRUE,
                       title = NULL,
                       is.sort = TRUE,
                       show_class_blocks = FALSE,
                       break_range = NULL,
                       break_scales = 0.5,

                       facet_split = FALSE
) {
  library(dplyr)
  library(rlang)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(forcats)
  if (!is.null(break_range)) library(ggbreak)

  plot_df <- purrr::map_dfr(
    data,
    ~ dplyr::select(.x, !!sym(x), !!sym(y.up), !!sym(y.down)),
    .id = "Class"
  ) %>%
    dplyr::filter(!is.na(!!sym(x)), !is.na(!!sym(y.up)), !is.na(!!sym(y.down))) %>%
    dplyr::mutate(Total = !!sym(y.up) + !!sym(y.down))
    plot_df$Class <- factor(plot_df$Class, levels = names(data))

  if (is.sort) {
    plot_df <- plot_df %>%
      dplyr::group_by(Class) %>%
      dplyr::mutate(!!sym(x) := forcats::fct_reorder(!!sym(x), Total, .desc = TRUE)) %>%
      dplyr::arrange(Class, desc(Total)) %>%
      dplyr::ungroup()
  }

  u <- dplyr::pull(plot_df, !!sym(y.up))
  d <- dplyr::pull(plot_df, !!sym(y.down))
  pval <- stats::wilcox.test(u, d, paired = paired)$p.value

  down_nm <- as_label(ensym(y.down))
  plot_long <- plot_df %>%
    dplyr::select(!!sym(x), !!sym(y.up), !!sym(y.down), Class) %>%
    tidyr::pivot_longer(
      cols = c(!!sym(y.up), !!sym(y.down)),
      names_to = "group",
      values_to = "value"
    ) %>%
    dplyr::mutate(value = dplyr::if_else(group == down_nm, -value, value))


  p <- ggplot(plot_long, aes(x = !!sym(x), y = value, fill = group)) +
    geom_col(width = 0.8) +

    scale_fill_manual(values = c("#C0D6EA", "#AABCDB")) +
    scale_y_continuous(labels = function(z) abs(z), expand = expansion(mult = c(0.05, 0.1))) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
    labs(
      title = if (is.null(title)) "Mirror Bar Chart" else title,
      subtitle = paste0("Paired Wilcoxon p = ", format(pval, digits = 3, scientific = TRUE)),
      x = NULL, y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      axis.ticks.y.left = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length.y = unit(0.15, "cm"),

      axis.line.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.y.right = element_blank(),

      axis.line.y.left = element_line(linewidth = 0.3, color = "black"),
      axis.line.x = element_line(linewidth = 0.3, color = "black"),

      axis.ticks.x = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length.x = unit(0.15, "cm"),

      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y = element_text(size = 5),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),

      legend.position = "top",
      legend.justification = "center",
      legend.margin = margin(t = 0, b = 5),
      legend.key.size = unit(0.4, "cm"),

      strip.background = element_blank(),
      strip.text = element_text(size = 10, face = "bold", margin = margin(b = 5))
    )

  if (facet_split) {
    p <- p + facet_grid(~ Class, scales = "free_x", space = "free_x")
  }

  if (!is.null(break_range)) {
    p <- p + ggbreak::scale_y_break(break_range, space = 0.1, scales = break_scales)
  }

  return(p)
}
