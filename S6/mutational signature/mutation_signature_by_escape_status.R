# S6: Mutation signature by escape status
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

source("/pub5/xiaoyun/BioX/Bioc/0.BioData/PanCancer/DataCenter/[[Management]]/Tools/[[[IntegratedData]]]/CombineDataMulti.XTeam.R")
t.evasion.mechanism = c("all.APM.mut", "IFNG.pathway.mut")
file.query.vec = c('PatientCenter/PatientCenter.rds',
    'Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds',
    'Results/BioGenomics/05.MutationSignatures/mutation_signatures.rds')
output.dir = "/WorkSpace/heshengyuan/01ImmuneEditing/10ImmuneEscape/50.ImmuneEscape/[Q]BaselineEffectsOfImmuneEscape/02ImmuneEscapeAndSignatureContribution/Revision"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
output.file = file.path(file.path(output.dir, "Heatmaps"), paste("[PanCancer][Heatmap]signature_contribution_FC_and_p_by_escape_group.pdf", sep = ""))

if (is.pass & file.exists(output.file)) {
    cat(blue(paste0("Results already exist; please check them")))
} else {
# Use CombineDataMulti.XTeam to process multiple datasets.
    source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[UnifiedDataType]]/GenomicData/MutationSignatures/AssociationAnalysis/signature_contribution_diff.R")
    all.data <- CombineDataMulti.XTeam(
        ID.xteams = ID.xteams,
        file.query.vec = file.query.vec,
        is.rbind = TRUE,
        process.fun = MutSigContributionDiff,
        group.cols = t.evasion.mechanism,
        fdr.method = "BH"
    )
    analysis.data <- all.data
    colnames(analysis.data) <- c("ID.xteam", "obs.col", "comparison_note", "value", "p.value", "p.adj", "feature")

# Draw the heatmap.
    source("/pub5/xiaoyun/BioX/PLOT/Heatmaps/SignificanceHeatmaps/significance_dot_heatmap.R")
    tt.output.dir = file.path(output.dir, "Heatmaps")
    if(!dir.exists(tt.output.dir)) { dir.create(tt.output.dir, showWarnings = FALSE, recursive = TRUE)  }
    t.output.file = file.path(tt.output.dir, paste("[PanCancer][Heatmap]signature_contribution_FC_and_p_by_escape_group.pdf", sep = ""))

    pdf(t.output.file, width = 12, height = 10)
    for(ev in t.evasion.mechanism){
        message("Processing escape mechanism: ", ev)
        t.df1 = dplyr::filter(analysis.data, feature == ev)
        p = PlotCircleHeatmapWithSignif(plot.data = t.df1, x = "ID.xteam", y = "obs.col", statistic.type = "value", statistic.p = "p.adj",
            Metric.category = NULL, add.signif.box = TRUE, title = ev)
        print(p)
    }
    dev.off()
}
