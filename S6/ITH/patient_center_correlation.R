# S6: ITH patient center correlation
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

PCAllTest <- function(
      patient.center,
      query.col = "BergerParkerDominance",
      obs.cols = c(
        "AgeDiagnosis", "Gender", "AnatomicSite_class",
        "Histologicaltype",
        "TNM.stage_4class",
        "MSISubtype_class", "SmokingStatus_class",
        "ResidualTumor", "VenousInvasion", "LymphaticInvasion", "PerineuralInvasion",
        "SynchronousCRC", "isColonPolyps"
      ),
      output.dir = NULL,
      return.plot = FALSE
){
    library(dplyr)

# Sample filtering and preparation.
    source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/ClinicalPhenotypes/0.Filter/FilterPatientCenter.R")

# Factor conversion.
    source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/ClinicalPhenotypes/Tools/DiscretizeRecodeFactorize/auto_factorize_patient_center.R")
    tmp <- PC.AutoFactorize(patient.center)
    t.patient.center <- tmp$patient.center

    source("/pub5/xiaoyun/BioX/R/BasicOperation/DataFrame/detect_column_types.R")
    t.check = Check_DisCon_Vector(PCCols(t.patient.center, query.col))
    if (t.check) {
      source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/ClinicalPhenotypes/AnalysisAndPlots/CorrelationTests/patient_center_group_differences.R")
# Difference tests (the underlying function selects anova/kruskal for continuous variables and chisq/fisher for categorical variables).
      result <- Diff.PC(
        patient.center = t.patient.center,
        group.col = query.col,
        obs.col = obs.cols,
        out.dir = output.dir,
        return.plot = return.plot
      )
    } else {
      source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/ClinicalPhenotypes/AnalysisAndPlots/CorrelationTests/patient_center_continuous_associations.R")
      result <- ComprehensiveCorr.PC(
        patient.center = t.patient.center,
        target.var = query.col,
        obs.col = obs.cols,
        out.dir = output.dir,
        return.plot = return.plot
      )
    }

    return(result)
  }
