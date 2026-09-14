# S8: ICB escape subgroups by genomic and transcriptomic burden
# This script evaluates immune escape features or burden-based subgroups in ICB
# cohorts. It tests associations with treatment response or survival and saves the
# resulting statistical summaries and figures.

library(magrittr); library(cluster); library(writexl); library(ggplot2)
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]AssociationBetweenEscapeGroupsAndICBOutcome/Question2"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
ICB.ID.xteams=c(
    "Liu_NatureMedicine_2019", "VanAllen_Science_2015", "Gide_CancerCell_2019", "Hugo_Cell_2016", "Freeman_CellReportsMedicine_2022", "Abbott_ClinicalCancerResearch_2021",
    "Ravi_NatureGenetics_2023", "Alban_NatureMedicine_2024")
profile.exprs=CombineData.XTeam(ICB.ID.xteams, file.query="OMICSData/Exprs.data.rds")
ICB.ID.xteams=ICB.ID.xteams[lengths(profile.exprs)!=0]
common.sam.list=lapply(ICB.ID.xteams, function(x){
    tmp=colnames(profile.exprs[[x]])
})      %>%     setNames(ICB.ID.xteams)

escape.data=CombineData.XTeam(ICB.ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ReferenceBasedBatchCorrectedEscapeProfile/ReferenceBasedBatchCorrectedEscapeProfile.rds")
Anagnostou.binary=readRDS("/IData2/DataCenter/LungCancer/Anagnostou_NatureCancer_2020/Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds")
Anagnostou.escape=readRDS("/IData2/DataCenter/LungCancer/Anagnostou_NatureCancer_2020/Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.rds")
Van_Hugo_binary=CombineData.XTeam(c("VanAllen_Science_2015", "Hugo_Cell_2016"), "Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
profile.escape.binary=lapply(ICB.ID.xteams, function(x){
    tmp.escape=escape.data[[x]]$orig.escape.data
    tmp=escape.data[[x]]$escape.binarization
    if(x=="Anagnostou_NatureCancer_2020"){
        tmp.escape=Anagnostou.escape
        tmp=Anagnostou.binary
    }
    tmp=tmp[lengths(tmp)!=0]
    tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
    if(x%in%c("VanAllen_Science_2015", "Hugo_Cell_2016")){
        pos=match(as.character(tmp$SampleID), colnames(Van_Hugo_binary[[x]]))
        tmp$B2M.bi.inactive=Van_Hugo_binary[[x]][1, pos]
    }
    tmp$all.APM.mut=tmp.escape$APM.alt$all.APM.mut[match(tmp$SampleID, tmp.escape$APM.alt$SampleID)]
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

# Apply the stated eligibility and data-quality restrictions before analysis.
profile.trans.combine.escape=lapply(ICB.ID.xteams, function(x){
    tmp=profile.escape.binary[[x]]        	%>%
        dplyr::group_by(SampleID)       %>%
        dplyr::mutate(
            AntigenPresentGene.downregulation=any(HLA.A, HLA.B, HLA.C, HLA.score, CALR, na.rm=TRUE),
            Checkpoint.overexprs=any(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm=TRUE),
            ImmunoSuppressiveCell.overexprs=any(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm=TRUE),
            ImmunoSuppressiveSig.overexprs=any(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm=TRUE),
            ImmuneActivationGene.downregulation=any(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm=TRUE))         %>%
        as.data.frame()
})      %>%     setNames(ICB.ID.xteams)

combine.features=c("Checkpoint.overexprs", "ImmunoSuppressiveCellSig.overexprs", "ImmuneActivationGene.downregulation", "AntigenPresentGene.downregulation")
genome.features=c(c("HLA.LOH", "B2M.bi.inactive", "all.APM.mut"), c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"), "CD274.deepAmp")
APMCheck.genomic.features=c(c("HLA.LOH", "B2M.bi.inactive", "all.APM.mut"), "CD274.deepAmp")
ActCheck.genomic.features=c(c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"), "CD274.deepAmp")
CheckSupprAct.nonG.features=c("Checkpoint.overexprs", "ImmunoSuppressiveCellSig.overexprs", "ImmuneActivationGene.downregulation")
CheckSupprAPM.nonG.features=c("Checkpoint.overexprs", "ImmunoSuppressiveCellSig.overexprs", "AntigenPresentGene.downregulation")
SupprActAPM.nonG.features=c("ImmunoSuppressiveCellSig.overexprs", "ImmuneActivationGene.downregulation", "AntigenPresentGene.downregulation")

tmp.data=Filter(function(x)  !all(is.na(x)), profile.trans.combine.escape)

cluster.data.list=lapply(names(tmp.data), function(x){
    cluster.data=lapply(list(c(combine.features, genome.features), combine.features, genome.features, APMCheck.genomic.features, ActCheck.genomic.features, CheckSupprAct.nonG.features,
         CheckSupprAPM.nonG.features, SupprActAPM.nonG.features), function(plus.features){

        if(length(intersect(colnames(tmp.data[[x]]), c(plus.features)))>0){
            tmp=tmp.data[[x]][, intersect(colnames(tmp.data[[x]]), c(plus.features, "SampleID"))]    %>%
                dplyr::filter(SampleID %in% common.sam.list[[x]])          %>%
                dplyr::mutate(
                    es.feature.num=rowSums(dplyr::select(., -SampleID))
                )       %>%
                dplyr::mutate(
                    es.feature.group=dplyr::case_when(
                        0<es.feature.num & es.feature.num<3 ~ "1-2",
                        es.feature.num>3 ~ "4+",
                        TRUE ~ as.character(es.feature.num)),
                    cluster=paste0("E", es.feature.group)
                )      %>%
                dplyr::select(SampleID, cluster, es.feature.num)
        }else{
            tmp=NULL
        }
        return(tmp)
    })          %>%     setNames(c("all.features", "nonGenomic", "Genomic", "APMCheck.G", "ActCheck.G", "CheckSupprAct.nonG", "CheckSupprAPM.nonG", "SupprActAPM.nonG"))
})      %>%     setNames(names(tmp.data))

saveRDS(cluster.data.list, file=file.path(Dir.output, "ICB.Es.features.combinedTrans.rds"))
