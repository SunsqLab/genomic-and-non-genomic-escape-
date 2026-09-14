# S6: ITH correlation statistics
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

if (length(query.col) == 1) {
  query.col.name = query.col[1]
} else {
  query.col.name = paste0(paste(query.col[1], query.col[2], sep = "_", collapse = "_"), "_", length(query.col))
}
source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/MatrixProfileData/BasicProcessing/extract_significant_differences.R")
input.data <- StandardDiff(input.data)
# Keep only cohort IDs present in input.data for subsequent display in the specified order.
ID.xteams.available <- ID.xteams[ID.xteams %in% unique(input.data$ID.xteam)]
input.data <- input.data %>%
  dplyr::mutate(metric_pair = paste(query.col, "vs", obs.col)) %>%
  dplyr::arrange(factor(ID.xteam, levels = ID.xteams.available))
if (!is.null(obs.cols.group)) {
  if (!is.null(names(obs.cols.group)) && all(names(obs.cols.group) != "")) {
    all.obs.cols.group.data <- data.frame(
      obs.col = unlist(obs.cols.group),
      group = factor(rep(names(obs.cols.group), lengths(obs.cols.group)), levels = names(obs.cols.group)),
      row.names = unlist(obs.cols.group)
    )
    input.data <- dplyr::left_join(input.data, all.obs.cols.group.data, by = c("obs.col"))
    input.data$obs.col <- factor(input.data$obs.col, levels = rev(all.obs.cols.group.data$obs.col))
    Metric.category <- "group"
  }
} else {
  Metric.category <- NULL
}
source("/pub5/xiaoyun/BioY/heshengyuan/0.Function/DataCenterFunctions/log2_transform_and_clip.R")
if("odds.ratio" %in% names(input.data)){
  input.data <- TruncateLog2value(input.data, process.col = "odds.ratio", limit = 3, out_col = "log2OR")
}
source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_dot_heatmap.R")
p <- PlotCircleHeatmapWithSignif(
  plot.data       = input.data,
  x               = "ID.xteam",
  y               = "obs.col",
  statistic.type  = NULL,
  statistic.p     = "p.adj",
  Metric.category = NULL,
  title           = "Significance dot heatmap"
)
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(output.dir, paste0(query.col.name, "_dot_heatmap_", "[significance].pdf")), p, width = 12, height = 10)
source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_tile_heatmap.R")
p <- PlotRectHeatmapWithSignif2(
  plot.data = input.data,
  x = "query.col",
  y = "obs.col",
  fill.col = "IsSig",
  facet.category = "ID.xteam",
  title = "Significance dot heatmap"
)
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(output.dir, paste0(query.col.name, "_tile_heatmap_", "[significance].pdf")), p, width = 12, height = 10)
source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_heatmap_with_barplot2.R")
row.groups = lapply(query.col, function(x) {
  return(paste0(x, " vs ", obs.cols))
})
names(row.groups) <- query.col
if (any(!is.na(input.data$log2FoldChange))) {
  heatmap.list <- PlotEffectHeatmap(
    data = input.data,
    row.col = "metric_pair",
    col.col = "ID.xteam",
    coef.col = "log2FoldChange",
    p.col = "p.adj",
    row.groups = row.groups,
    legend.title = "Log2 Fold Change",
    cluster.cols = FALSE, column.names.rot = 90, show.sig.barplot = TRUE,
    output.dir = output.dir,
    output.file.name = paste0(query.col.name, "_metric_association_heatmap_", "[FCsignificance]")
  )
}
# If correlation has values, draw a correlation heatmap.
if (any(!is.na(input.data$correlation))) {
  heatmap.list <- PlotEffectHeatmap(
    data = input.data,
    row.col = "metric_pair",
    col.col = "ID.xteam",
    coef.col = "correlation",
    p.col = "p.adj",
    row.groups = row.groups,
    legend.title = "Correlation",
    cluster.cols = FALSE, column.names.rot = 90, show.sig.barplot = TRUE,
    output.dir = output.dir,
    output.file.name = paste0(query.col.name, "_metric_association_heatmap_", "[correlationsignificance]")
  )
}
if (any(!is.na(input.data$log2FoldChange))) {
  p <- PlotCircleHeatmapWithSignif(
    plot.data       = input.data,
    x               = "ID.xteam",
    y               = "obs.col",
    statistic.type  = "log2FoldChange",
    statistic.p     = "p.adj",
    Metric.category = Metric.category,
    title           = "Significance dot heatmap"
  )
  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output.dir, paste0(query.col.name, "_dot_heatmap_", "[log2FoldChangesignificance].pdf")), p, width = 12, height = 10)
}
if (any(!is.na(input.data$odds.ratio))) {
  source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_hotspot_plot.R")
  input.data$OR <- factor(
    ifelse(input.data$odds.ratio > 1, "OR>1", "OR<1"),
    levels = c("OR>1", "OR<1")
  )
  p1 <- PlotPointHeatmapWithSignif(
    plot.data = input.data,
    p.type = "p.adj",
    x = "obs.col",
    y = "ID.xteam",
    color = "OR"
  )
  plot.filename <- file.path(output.dir, paste0(query.col.name, "_metric_association_heatmap_", "[ORsignificance].pdf"))
  ggsave(plot.filename, p1, width = 16, height = 10)
}
if (any(!is.na(input.data$odds.ratio)) & !is.null(obs.cols.group) & !is.null(Metric.category)) {
  source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_hotspot_plot.R")
  input.data$OR <- factor(
    ifelse(input.data$odds.ratio > 1, "OR>1", "OR<1"),
    levels = c("OR>1", "OR<1")
  )
  p1 <- PlotPointHeatmapWithSignif(
    plot.data = input.data,
    p.type = "p.adj",
    x = "obs.col",
    y = "ID.xteam",
    group = "group",
    color = "OR",
    facet = "cols"
  )
  plot.filename <- file.path(output.dir, paste0(query.col.name, "_metric_association_heatmap_", "[ORsignificance].pdf"))
  ggsave(plot.filename, p1, width = 16, height = 10)
}
if (any(!is.na(input.data$log2FoldChange)) & !is.null(obs.cols.group) & !is.null(Metric.category)) {
  source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_hotspot_plot.R")
  input.data$log2FC <- factor(
    ifelse(input.data$log2FoldChange > 0, "log2FC>0", "log2FC<0"),
    levels = c("log2FC>0", "log2FC<0")
  )
  p <- PlotPointHeatmapWithSignif(
    plot.data = input.data,
    p.type = "p.adj",
    x = "ID.xteam",
    y = "obs.col",
    color = "log2FC",
    size.value = "3class",
    group = Metric.category,
    plot.function = "dot_plot2"
  )
  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output.dir, paste0(query.col.name, "_dot_plot_", "[log2FoldChangesignificance].pdf")), p, width = 12, height = 10)
}
if (any(!is.na(input.data$log2FoldChange)) & !is.null(obs.cols.group) & !is.null(Metric.category)) {
  source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_square_heatmap.R")
  p_rect <- PlotRectHeatmapWithSignif(
    plot.data = input.data,
    x = "ID.xteam",
    y = "obs.col",
    statistic.type = "log2FoldChange",
    statistic.p = "p.adj",
    Metric.category = Metric.category,
    tag = paste0(query.col.name, "_Relationship_RectangularHeatmap"),
    right_plot_type = "proportion"
  )
  plot.filename <- file.path(output.dir, paste0(query.col.name, "_metric_association_heatmap_", "[log2FoldChangesignificancesquare_heatmap].pdf"))
  ggsave(plot.filename, p_rect, width = 12, height = 10)
}
if (any(!is.na(input.data$log2FoldChange)) & !is.null(obs.cols.group) & !is.null(Metric.category)) {
  source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_dot_heatmap.R")
  p_circle <- PlotCircleHeatmapWithSignif(
    plot.data = input.data,
    x = "ID.xteam",
    y = "obs.col",
    statistic.type = "log2FoldChange",
    statistic.p = "p.adj",
    Metric.category = Metric.category,
    title = paste0(query.col.name, "_Relationship_CircularHeatmap")
  )
  plot.filename <- file.path(output.dir, paste0(query.col.name, "_metric_association_heatmap_", "[log2FoldChangesignificancecircular_heatmap].pdf"))
  ggsave(plot.filename, p_circle, width = 12, height = 10)
}

if (any(!is.na(input.data$odds.ratio))   & exists("obs.cols.group") & exists("Metric.category")) {
  source("/pub5/xiaoyun/BioY/heshengyuan/1.ImmuneEditing/09ImmuneEditingInSelectedGenomicRegions/00Function/significance_dot_heatmap.R")
  input.data$log2OR[is.infinite(input.data$log2OR) | is.nan(input.data$log2OR)] <- NA
  p_circle <- PlotCircleHeatmapWithSignif(
    plot.data = input.data,
    x = "ID.xteam",
    y = "obs.col",
    statistic.type = "log2OR",
    statistic.p = "p.adj",
    Metric.category = Metric.category,
    title = paste0(query.col.name, "_Relationship_CircularHeatmap_OR"),
    add.signif.box = TRUE
  )
  plot.filename <- file.path(output.dir, paste0(query.col.name, "_metric_association_heatmap_", "[ORsignificancecircular_heatmap].pdf"))
  ggsave(plot.filename, p_circle, width = 12, height = 10)
}
