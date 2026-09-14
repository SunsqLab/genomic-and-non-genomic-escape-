# S7: Genomic escape subgroup heatmap
# This script defines immune escape subgroups and evaluates their molecular or
# clinical characteristics across TCGA cohorts. It performs the indicated subgroup
# comparisons or survival models and saves the resulting figures.

library(magrittr); library(cluster); library(ggplot2); library(ComplexHeatmap); library(circlize); library(openxlsx)
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]Multiple_immune_escape_subgroup_patterns/Question1"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R")
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c("FPPP_TCGA", "LAML_TCGA")))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.escape=CombineData.XTeam(ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.rds")
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds")
profile.escape.binarization=lapply(ID.xteams, function(x){
        data=profile.escape.binarization[[x]]
        tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), data)
        tmp$HLA.mut=profile.escape[[x]]$APM.alt$HLA.mut[match(tmp$SampleID, profile.escape[[x]]$APM.alt$SampleID)]
        tmp$all.APM.mut=profile.escape[[x]]$APM.alt$all.APM.mut[match(tmp$SampleID, profile.escape[[x]]$APM.alt$SampleID)]
        tmp$B2M.inactive=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
        tmp$B2M.inactive=ifelse(tmp$B2M.inactive==1, TRUE, FALSE)
        return(tmp)
})      %>%     setNames(ID.xteams)
profile.escape.binarization$PanCancer_TCGA=do.call(rbind, profile.escape.binarization)

source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS("/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds")
profile.clini=GetInfor.PatientCenter(patient.center, colNames=c("TMB", "CancerType"))

profile.Immunecluster=read.xlsx(file.path(Dir.output, "Immune landscape of cancer/NIHMS958212-supplement-2.xlsx"))
profile.Immunecluster=profile.Immunecluster[, c("TCGA.Participant.Barcode", "Immune.Subtype")]  %>%
    setNames(c("PatientID", "ImmuneCluster"))

profile.exprs=CombineData.XTeam(ID.xteams, file.query='OMICSData/Exprs.data.rds')
profile.exprs$CRC_TCGA=do.call(cbind, CombineData.XTeam(c('COAD_TCGA', 'READ_TCGA'), file.query='OMICSData/Exprs.data.rds'))
source('/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/ImmuneEscape/ImmuneEscapeMechanisms/ImmuneCellInfiltration/ImmuneCellInfiltrationScores.R')
profile.Danaher.cell=lapply(profile.exprs, function(x){
    ImmuneInfiltraScore(exprs.data=x)
})      %>%     setNames(ID.xteams)
profile.Danaher.cell$PanCancer_TCGA=do.call(rbind, profile.Danaher.cell)

outDir=file.path(Dir.output, "All_escape_features_binary_distance")
profile.cluster=readRDS(file.path(outDir, "Es.features.sum.combinedTrans.20260510.cluster.rds"))

genome.features.list=list(
    APMalt=c("HLA.LOH", "B2M.inactive", "HLA.mut", "B2M.mut", "all.APM.mut"),
    Checkalt=c("CD274.deepAmp"),
    ActMalt=c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut")
)
genome.features=unlist(genome.features.list)
names(genome.features)=rep(c("#c82621", "#F6C141", "#fa8b69ff"), lengths(genome.features.list))

cell.types=c("Cytotoxic.cell", "CD8.Tcell", "NK.cell")

heatmap.list=lapply(names(profile.escape.binarization), function(x){
    combined.escape=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), list(profile.escape.binarization[[x]], profile.clini, profile.cluster$Genomic, profile.Danaher.cell[[x]]))
    tmp.data=Filter(function(x)  !all(is.na(x)), combined.escape)
    data_pre=tmp.data[, c(intersect(genome.features, colnames(tmp.data)), "SampleID", "TMB", "CancerType", "cluster", cell.types)]   %>%
        dplyr::filter(if_all(all_of(c("cluster", "TMB", cell.types)), ~ !is.na(.)))
    data_pre$ImmuneCluster=profile.Immunecluster$ImmuneCluster[match(substr(data_pre$SampleID, 1, 12), profile.Immunecluster$PatientID)]
    data_pre$ImmuneCluster[is.na(data_pre$ImmuneCluster)]="Other"
    data_pre=data_pre[, colSums(is.na(data_pre)) == 0]

    anno.data=data_pre

    heatmap.data=data_pre[, intersect(colnames(data_pre), c(genome.features, cell.types))]      %>%
        dplyr::mutate(all.APM.mut=ifelse(all.APM.mut==0, FALSE, TRUE), HLA.mut=ifelse(HLA.mut==0, FALSE, TRUE))      %>%
        dplyr::mutate(across(where(is.logical), ~ ifelse(is.na(.), NA, ifelse(., "Yes", ""))))

    cancer_colors=setNames(c("#1f77b4", "#00b4d8", "#90e0ef", "#caf0f8", "#fb6f92", "#ff8fab", "#ffc2d1", "#ffe5ec",
            "#38a3a5", "#57cc99", "#80ed99", "#c7f9cc", "#5e548e", "#9f86c0", "#be95c4", "#e0b1cb",
            "#2b9348", "#80b918", "#d4d700", "#eeef20", "#bb3e03", "#ca6702", "#ee9b00", "#e9d8a6",
            "#51ccd1", "#8be8d7", "#a0f1da", "#b4fadc", "#9d6b53", "#cd9777", "#deab90", "#edc4b3")[1:length(unique(anno.data$CancerType))], unique(anno.data$CancerType))
    immue.cluster.color=c("C1"="#e63946", "C2"="#ec9a9a", "C3"="#f1faee", "C4"="#a8dadc", "C5"="#457b9d", "C6"="#1d3557", "Other"="lightgrey")
    col_anno=HeatmapAnnotation(
        cluster=anno_block(gp=gpar(fill=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"), col="black"), height=unit(5, "mm")),
        CancerType=anno.data$CancerType   %>%     as.factor(),
        ImmuneCluster=anno.data$ImmuneCluster          %>%     as.factor(),

        col=list(CancerType=cancer_colors, ImmuneCluster=immue.cluster.color)
    )

    alter_fun=list(
        background=function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                    gp=gpar(fill="#EEEEEE", col=NA))
        },
        Yes=function(x, y, w, h) {
            grid.rect(x, y, w*1, h*1, gp=gpar(fill="#376f9cff", col=NA))
        }
    )

    tmp.matrix=heatmap.data[, intersect(unlist(genome.features.list), colnames(heatmap.data))]
    heatmap1=oncoPrint(
        t(tmp.matrix),
        alter_fun=alter_fun,
        col=c("Yes"="#376f9cff"),
        column_split=data_pre$cluster,
        top_annotation=col_anno,
        column_gap=unit(2, "mm"),
        cluster_column_slices=TRUE,
        right_annotation=NULL,

        show_column_names=FALSE,
        heatmap_legend_param=list(title="isEscape"),
        row_title=paste0(x, " (n=", nrow(tmp.matrix), ") "),
        row_names_gp=gpar(col=names(genome.features)[match(colnames(tmp.matrix), genome.features)]),
        na_col="grey"
    )

    tmp.matrix=scale(heatmap.data[, cell.types])
    col_range=colorRamp2(c(min(tmp.matrix, na.rm=T), min(tmp.matrix, na.rm=T)/2, 0, max(tmp.matrix, na.rm=T)/2, max(tmp.matrix, na.rm=T)), c("# Prepare annotations and generate the heatmap.
    heatmap2=Heatmap(
            t(tmp.matrix), name="Cell.infiltra",
            show_column_names=FALSE,
            cluster_columns=FALSE,
            cluster_rows=FALSE,
            show_column_dend=FALSE,
            col=col_range)

    combined_heatmap=Reduce(`%v%`, list(heatmap1, heatmap2))
})

pdf(file.path(outDir, "1.genomic_escape_subgroup_pattern_based_on_binarized_escape_burden_combinedtranscriptomic.pdf"), 25, 6)
print(heatmap.list)
dev.off()
