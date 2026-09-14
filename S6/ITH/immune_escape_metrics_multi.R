# S6: ITH immune escape metrics multi
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

ImmuneEscape_ColRelationX <- function(
    ID.xteams,
    query.col = "MATH",
    obs.cols = NULL,
    obs.cols.group = list(
      APM.alt = c("HLA.LOH", "B2M.LOH", "HLA.mut", "B2M.mut", "all.APM.mut"),
      Checkalt = c("CD274.deepAmp"),
      ActMalt = c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
      Check.exp = c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
      SupprCell = c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
      SupprSig = c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
      ActGenedown = c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
      APMdown = c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
      Neo = c("mean.neo.exprs")
    ),
    output.dir = NULL
    ) {
  require(dplyr)
  if (is.null(obs.cols) & !is.null(obs.cols.group)) {
    obs.cols = unique(unlist(obs.cols.group))
  }

# Integrate PatientCenter data from multiple datasets and preprocess it.
  source("/pub5/xiaoyun/BioX/Bioc/0.BioData/PanCancer/DataCenter/[[Management]]/Tools/[[[IntegratedData]]]/CombineDataMulti.XTeam.R")
  source("/pub5/xiaoyun/BioX/Bioc/0.BioData/PanCancer/DataCenter/[BioImmune]/[[[Questions]]]/[Q]ImmuneEscapeProfile/immune_escape_metrics_single.R")
  input.data = CombineDataMulti.XTeam(ID.xteams,
    file.query = c("PatientCenter/PatientCenter.rds", "Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds"),
    is.rbind = TRUE,
    process.fun = ImmuneEscape_ColRelation,
    query.col = query.col, obs.cols = obs.cols,
    obs.cols.group = obs.cols.group,
    output.dir = NULL
  )

  source("/pub5/xiaoyun/BioY/heshengyuan/1.ImmuneEditing/10ImmuneEscape/01ImmuneEscapeAndITH/[MultiDataIntegration]ImmuneEscapeAndContinuousMetrics/correlation_statistics.R", local = TRUE)

  return(input.data)
}
