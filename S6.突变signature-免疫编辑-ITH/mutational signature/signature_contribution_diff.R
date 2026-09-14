# S6: Signature contribution differences
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/ClinicalPhenotypes/AnalysisAndPlots/CorrelationTests/patient_center_continuous_group_differences.R")

MutSigContributionDiff <- function(patient.center, group.data, group.cols=NULL, signature.data, fdr.method="BH") {
# Extract contribution data from the signature object as continuous variables for association analysis.

    sig.data.used = signature.data[["COSMIC.SData"]][["Contribution"]]
    sig.data.used[["SampleID"]] = rownames(sig.data.used)

# Call Diff.PC2 to analyze signature-contribution differences across binary-feature groups.

    analysis.results = Diff.PC2(
        patient.center = patient.center,
        group.data = group.data,
        group.cols = group.cols,
        feature.data = sig.data.used,
        feature.cols = NULL,
        fdr.method = fdr.method
    )
    analysis.results = as.data.frame(analysis.results)
    return(analysis.results)
}

if(FALSE){

    test1 = pbapply::pblapply(setNames(ID.xteams, ID.xteams), cl = 4, function(ID.xteam){
        tt1 = readRDS(file.path("/IData/DataCenter/TCGA", ID.xteam, "PatientCenter/PatientCenter.rds"))
        tt2 = readRDS(file.path("/IData/DataCenter/TCGA", ID.xteam, "Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds"))
        tt3 = readRDS(file.path("/IData/DataCenter/TCGA", ID.xteam, "Results/BioGenomics/05.MutationSignatures/mutation_signatures.rds"))

        test = MutSigContributionDiff(patient.center = tt1, group.data = tt2, group.cols=t.evasion.mechanism, signature.data = tt3, fdr.method="BH")

        return(test)
    })

}
