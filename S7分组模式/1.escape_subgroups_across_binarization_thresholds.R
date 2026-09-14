# S7: Escape subgroups across binarization thresholds
# This script defines immune escape subgroups and evaluates their molecular or
# clinical characteristics across TCGA cohorts. It performs the indicated subgroup
# comparisons or survival models and saves the resulting figures.

library(magrittr); library(cluster); library(writexl); library(ggplot2)
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]Multiple_immune_escape_subgroup_patterns/Question1"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R")
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c("FPPP_TCGA", "LAML_TCGA")))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.escape=CombineData.XTeam(ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.rds")

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
patient.center=readRDS("/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds")
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds")
profile.escape.binarization=lapply(ID.xteams, function(x){
    data=profile.escape.binarization[[x]]
# Retain primary tumor samples only.
    tmp=lapply(data, function(y){
                y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
            })
    tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
    tmp$all.APM.mut=profile.escape[[x]]$APM.alt$all.APM.mut[match(tmp$SampleID, profile.escape[[x]]$APM.alt$SampleID)]
    tmp$biallelic.B2M.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
    tmp$biallelic.B2M.inactivation=ifelse(tmp$biallelic.B2M.inactivation==1, TRUE, FALSE)
    return(tmp)
})     %>%      setNames(ID.xteams)

profile.trans.combine.escape=profile.escape.binarization 	%>%
    do.call(what=rbind)     %>%
    dplyr::group_by(SampleID)    %>%
    dplyr::mutate(
        AntigenPresentGene.downregulation=any(HLA.A, HLA.B, HLA.C, HLA.score, CALR, na.rm=TRUE),
        Checkpoint.overexprs=any(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm=TRUE),
        ImmunoSuppressiveCell.overexprs=any(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm=TRUE),
        ImmunoSuppressiveSig.overexprs=any(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm=TRUE),
        ImmuneActivationGene.downregulation=any(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm=TRUE))         %>%
    as.data.frame()

# Apply the specified binarization threshold to define immune escape groups.
profile.escape.median.binary=CombineData.XTeam(ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.median.rds")
profile.escape.median.binary=lapply(ID.xteams, function(x){
    data=profile.escape.median.binary[[x]]
    tmp=lapply(data, function(y){
                y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
            })
    tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
    tmp$all.APM.mut=profile.escape[[x]]$APM.alt$all.APM.mut[match(tmp$SampleID, profile.escape[[x]]$APM.alt$SampleID)]
    tmp$biallelic.B2M.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
    tmp$biallelic.B2M.inactivation=ifelse(tmp$biallelic.B2M.inactivation==1, TRUE, FALSE)
    return(tmp)
})     %>%      setNames(ID.xteams)

profile.trans.median.escape=profile.escape.median.binary 	%>%
    do.call(what=rbind)     %>%
    dplyr::group_by(SampleID)    %>%
    dplyr::mutate(
        AntigenPresentGene.downregulation=any(HLA.A, HLA.B, HLA.C, HLA.score, CALR, na.rm=TRUE),
        Checkpoint.overexprs=any(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm=TRUE),
        ImmunoSuppressiveCell.overexprs=any(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm=TRUE),
        ImmunoSuppressiveSig.overexprs=any(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm=TRUE),
        ImmuneActivationGene.downregulation=any(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm=TRUE))         %>%
    as.data.frame()

profile.escape.2sd.binary=CombineData.XTeam(ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.Mean2SD.rds")
profile.escape.2sd.binary=lapply(ID.xteams, function(x){
    data=profile.escape.2sd.binary[[x]]
    tmp=lapply(data, function(y){
                y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
            })
    tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
    tmp$all.APM.mut=profile.escape[[x]]$APM.alt$all.APM.mut[match(tmp$SampleID, profile.escape[[x]]$APM.alt$SampleID)]
    tmp$biallelic.B2M.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
    tmp$biallelic.B2M.inactivation=ifelse(tmp$biallelic.B2M.inactivation==1, TRUE, FALSE)
    return(tmp)
})     %>%      setNames(ID.xteams)

profile.trans.2sd.escape=profile.escape.2sd.binary 	%>%
    do.call(what=rbind)     %>%
    dplyr::group_by(SampleID)    %>%
    dplyr::mutate(
        AntigenPresentGene.downregulation=any(HLA.A, HLA.B, HLA.C, HLA.score, CALR, na.rm=TRUE),
        Checkpoint.overexprs=any(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm=TRUE),
        ImmunoSuppressiveCell.overexprs=any(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm=TRUE),
        ImmunoSuppressiveSig.overexprs=any(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm=TRUE),
        ImmuneActivationGene.downregulation=any(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm=TRUE))         %>%
    as.data.frame()

profile.exprs=CombineData.XTeam(ID.xteams, file.query="OMICSData/Exprs.data.rds")
profile.mut=CombineData.XTeam(ID.xteams, file.query="OMICSData/Mutations.data.rds")
profile.CN=CombineData.XTeam(ID.xteams, file.query="OMICSData/CN.rds")
common.sam=lapply(ID.xteams, function(x){
        Reduce(intersect, list(colnames(profile.exprs[[x]]), unique(profile.mut[[x]]$SampleID), unique(profile.CN[[x]]$SampleID)))
})          %>%     unlist()

combine.features=list(Checkexp="Checkpoint.overexprs", Supprcell="ImmunoSuppressiveCell.overexprs", Supprsig="ImmunoSuppressiveSig.overexprs",
    Actdown="ImmuneActivationGene.downregulation", APMdown="AntigenPresentGene.downregulation")
genome.features=list(APMalt=c("HLA.LOH", "biallelic.B2M.inactivation", "all.APM.mut"), ActMalt=c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"), Checkalt=c("CD274.deepAmp"))
all.features=c(combine.features, genome.features)

APMCheck.genomic.features=c(c("HLA.LOH", "biallelic.B2M.inactivation", "all.APM.mut"), "CD274.deepAmp")
ActCheck.genomic.features=c(c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"), "CD274.deepAmp")
CheckSupprAct.nonG.features=c("Checkpoint.overexprs", "ImmunoSuppressiveCell.overexprs", "ImmunoSuppressiveSig.overexprs", "ImmuneActivationGene.downregulation")
CheckSupprAPM.nonG.features=c("Checkpoint.overexprs", "ImmunoSuppressiveCell.overexprs", "ImmunoSuppressiveSig.overexprs", "AntigenPresentGene.downregulation")
SupprActAPM.nonG.features=c("ImmunoSuppressiveCell.overexprs", "ImmunoSuppressiveSig.overexprs", "ImmuneActivationGene.downregulation", "AntigenPresentGene.downregulation")
CheckActAPM.nonG.features=c("ImmuneActivationGene.downregulation", "AntigenPresentGene.downregulation")
CheckSupprsigActAPM.nonG.features=c("ImmunoSuppressiveSig.overexprs", "ImmuneActivationGene.downregulation", "AntigenPresentGene.downregulation")
CheckSupprcellActAPM.nonG.features=c("ImmunoSuppressiveCell.overexprs", "ImmunoSuppressiveSig.overexprs", "ImmuneActivationGene.downregulation", "AntigenPresentGene.downregulation")

leave.one.category.out.list=lapply(1:length(all.features), function(sort.feature) {
    drop.features=all.features[[sort.feature]]
    keep.features=setdiff(unlist(all.features), drop.features)
}) %>%      setNames(paste0("remove_", names(all.features)))

tmp.data=Filter(function(x)  !all(is.na(x)), profile.trans.combine.escape)
distinct.combine.features=c(list(all.features, genome.features, APMCheck.genomic.features, ActCheck.genomic.features,
    combine.features, CheckSupprAct.nonG.features, CheckSupprAPM.nonG.features, SupprActAPM.nonG.features,
    CheckActAPM.nonG.features, CheckSupprsigActAPM.nonG.features), leave.one.category.out.list)
names(distinct.combine.features)=c("all.features", "Genomic", "APMaltCheckalt", "ActMaltCheckalt", "nonGenomic", "CheckexpSupprActdown", "CheckexpSupprAPMdown",
    "SupprActAPMdown", "CheckexpActAPMdown", "CheckexpSupprsigActAPMdown", names(leave.one.category.out.list))

cluster.data=lapply(distinct.combine.features, function(plus.features){
    plus.features=unlist(plus.features)
    tmp.data[, intersect(colnames(tmp.data), c(plus.features, "SampleID"))]    %>%
        dplyr::filter(SampleID %in% common.sam)          %>%
        dplyr::mutate(
            Escape.burden=rowSums(dplyr::select(., -SampleID))
        )       %>%
        dplyr::mutate(
            es.feature.group=dplyr::case_when(
                0<Escape.burden & Escape.burden<3 ~ "1-2",
                Escape.burden>3 ~ "4+",
                TRUE ~ as.character(Escape.burden)),
            cluster=paste0("E", es.feature.group)
        )      %>%
        dplyr::select(SampleID, cluster, Escape.burden)
})          %>%     setNames(names(distinct.combine.features))

tmp.median.data=Filter(function(x)  !all(is.na(x)), profile.trans.median.escape)

cluster.median.data=lapply(distinct.combine.features, function(plus.features){
    plus.features=unlist(plus.features)
    tmp.median.data[, intersect(colnames(tmp.median.data), c(plus.features, "SampleID"))]    %>%
        dplyr::filter(SampleID %in% common.sam)          %>%
        dplyr::mutate(
            Escape.burden=rowSums(dplyr::select(., -SampleID))
        )       %>%
        dplyr::mutate(
            es.feature.group=dplyr::case_when(
                0<Escape.burden & Escape.burden<3 ~ "1-2",
                Escape.burden>3 ~ "4+",
                TRUE ~ as.character(Escape.burden)),
            cluster=paste0("E", es.feature.group)
        )      %>%
        dplyr::select(SampleID, cluster, Escape.burden)
})          %>%     setNames(names(distinct.combine.features))

tmp.2sd.data=Filter(function(x)  !all(is.na(x)), profile.trans.2sd.escape)

cluster.2sd.data=lapply(distinct.combine.features, function(plus.features){
    plus.features=unlist(plus.features)
    tmp.2sd.data[, intersect(colnames(tmp.2sd.data), c(plus.features, "SampleID"))]    %>%
        dplyr::filter(SampleID %in% common.sam)          %>%
        dplyr::mutate(
            Escape.burden=rowSums(dplyr::select(., -SampleID))
        )       %>%
        dplyr::mutate(
            es.feature.group=dplyr::case_when(
                0<Escape.burden & Escape.burden<3 ~ "1-2",
                Escape.burden>3 ~ "4+",
                TRUE ~ as.character(Escape.burden)),
            cluster=paste0("E", es.feature.group)
        )      %>%
        dplyr::select(SampleID, cluster, Escape.burden)
})          %>%     setNames(names(distinct.combine.features))

outDir=file.path(Dir.output, "All_escape_features_binary_distance")
saveRDS(cluster.data, file=file.path(outDir, "Es.features.sum.combinedTrans.20260510.cluster.rds"))

saveRDS(cluster.median.data, file=file.path(outDir, "Es.features.median.Trans.20260610.cluster.rds"))
saveRDS(cluster.2sd.data, file=file.path(outDir, "Es.features.2sd.Trans.20260610.cluster.rds"))
