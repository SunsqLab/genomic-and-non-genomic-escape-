# S6: ITH immune escape metrics single
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

ImmuneEscape_ColRelation <- function(
    patient.center,
    ImmuneEscapeProfile.binary,
    query.col = "MATH",
    obs.cols = NULL,
    evasion.list.all = list(
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
    adjust.method = "BH",
    output.dir = NULL,
    num.as.continue_query.col=FALSE,
    num.as.continue_obs.cols=FALSE
    )
  {
  library(dplyr)

  source("/pub5/xiaoyun/BioX/Bioc/0.BioData/PanCancer/DataCenter/[BioImmune]/[[[Questions]]]/[Q]ImmuneEscapeProfile/[2.Tools]/merge_binary_escape_profile_with_patient_center.R")
  t.patient.center <- AppendBinaryFeatures(patient.center, ImmuneEscapeProfile.binary, evasion.list.all)
  if (is.null(obs.cols)) {
    obs.cols = unlist(evasion.list.all)
  }

  source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/ClinicalPhenotypes/AnalysisAndPlots/CorrelationTests/patient_center_correlation.R")
  if (length(query.col) == 1) {
    result <- PCAllTest(t.patient.center, query.col = query.col, obs.cols = obs.cols, output.dir = output.dir, adjust.method = adjust.method, num.as.continue_query.col=num.as.continue_query.col, num.as.continue_obs.cols=num.as.continue_obs.cols)
  } else {

    result = lapply(query.col, function(x) {
      result <- PCAllTest(t.patient.center, query.col = x, obs.cols = obs.cols, output.dir = output.dir, adjust.method = adjust.method, num.as.continue_query.col=num.as.continue_query.col, num.as.continue_obs.cols=num.as.continue_obs.cols)
      return(result)
    })
    result = do.call(rbind, result)

    result$p.adj = p.adjust(result$p.value, method = adjust.method)
  }
  return(result)
}
