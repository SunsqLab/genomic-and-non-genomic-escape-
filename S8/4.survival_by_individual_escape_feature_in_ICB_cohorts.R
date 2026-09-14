# S8: Survival by individual escape feature in ICB cohorts
# This script evaluates immune escape features or burden-based subgroups in ICB
# cohorts. It tests associations with treatment response or survival and saves the
# resulting statistical summaries and figures.

library(magrittr); library(ggplot2); library(readxl); library(survival); library(survminer); library(cowplot)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]AssociationBetweenEscapeGroupsAndICBOutcome/Question2'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
ICB.ID.xteams=c(
    "Liu_NatureMedicine_2019", "VanAllen_Science_2015", "Gide_CancerCell_2019", "Hugo_Cell_2016", "Freeman_CellReportsMedicine_2022", "Abbott_ClinicalCancerResearch_2021",
    "Ravi_NatureGenetics_2023", "Alban_NatureMedicine_2024")
escape.data=CombineData.XTeam(ICB.ID.xteams, file.query="Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ReferenceBasedBatchCorrectedEscapeProfile/ReferenceBasedBatchCorrectedEscapeProfile.rds")
Anagnostou.binary=readRDS("/IData2/DataCenter/LungCancer/Anagnostou_NatureCancer_2020/Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds")
Van_Hugo_binary=CombineData.XTeam(c("VanAllen_Science_2015", "Hugo_Cell_2016"), "Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
profile.escape.binary=lapply(ICB.ID.xteams, function(x){
    tmp=escape.data[[x]]$escape.binarization
    if(x=="Anagnostou_NatureCancer_2020"){
        tmp=Anagnostou.binary
    }
    if(x=="Braun_NatureMedicine_2020"){
        tmp=tmp[c("APM.alt", "Factor.Genomic")]
    }
    tmp=tmp[lengths(tmp)!=0]
    tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
    if(x%in%c("VanAllen_Science_2015", "Hugo_Cell_2016")){
        pos=match(as.character(tmp$SampleID), colnames(Van_Hugo_binary[[x]]))
        tmp$B2M.bi.inactive=Van_Hugo_binary[[x]][1, pos]
        tmp$B2M.bi.inactive=ifelse(tmp$B2M.bi.inactive==1, TRUE, FALSE)
    }
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
patient.center=CombineData.XTeam(ICB.ID.xteams, file.query="PatientCenter/PatientCenter.rds")
profile.ICB.clinic=lapply(ICB.ID.xteams, function(x){
    if(x!="Gide_CancerCell_2019"){
        tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("TreatmentResponse", "SampleTime", "TMB", "OS", "OS.time", "TreatmentHistory", "CancerType"))   %>%
# Retain eligible samples and remove records that do not meet the analysis criteria.
            dplyr::filter(TreatmentHistory!="[NA(8)]")
    }else{
        tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("TreatmentResponse", "SampleTime", "OS", "OS.time", "TreatmentHistory", "CancerType"))
        tmp$TMB=NA
    }
    if(x=="Liu_NatureMedicine_2019"){
        tmp$SampleTime="Pre-Treatment"
    }
    tmp=tmp[grepl("PRE", toupper(tmp$SampleTime)), ]
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

color.list=list("Skin cancer"=c("grey", "#ce3375"), "Lung cancer"=c("grey", "#60c8b3"))

all.features.list=list(
    APMalt=c("HLA.LOH", "B2M.bi.inactive", "HLA.mut", "B2M.mut", "all.APM.mut"),
    Checkalt=c("CD274.deepAmp"),
    ActMalt=c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
    Checkexp=c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    Supprcell=c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    Supprsig=c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    Actdown=c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown=c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo=c("mean.neo.exprs")
)
all.features=unlist(all.features.list)

pdf(file.path(Dir.output, paste0("2.1.individual_escape_featureeffect_onICBOS.pdf")), width=35, height=25)
plot.data=lapply(ICB.ID.xteams, function(x){
    if(!is.null(profile.ICB.clinic[[x]])){
        data=dplyr::full_join(profile.escape.binary[[x]], profile.ICB.clinic[[x]], by="SampleID")          %>%
            dplyr::filter(!is.na(SampleTime))       %>%
            dplyr::mutate(
                OS=ifelse(OS.time > 365*5, 0, OS),
                OS.time=ifelse(OS.time > 365*5, 365*5, OS.time)
            )

        plots=lapply(intersect(all.features, colnames(data)), function(feature){
            if(length(intersect(all.features, colnames(data)))>0 && length(table(data[, feature]))==2 && min(table(data[, feature]))>5){
                formula=as.formula(paste("Surv(OS.time, OS) ~", feature))
                fit=survfit(formula, data=data)
                fit$call$formula=formula

                p=ggsurvplot(fit, data=data, pval=TRUE, conf.int=TRUE, risk.table=TRUE,
                    tables.y.text=FALSE, legend.title="", palette=color.list[[data$CancerType[1]]], title=x)
            }else{
                NULL
            }
        })
        print(arrange_ggsurvplots(plots[lengths(plots)!=0], ncol=10,  nrow=5, print=F))

    }
})
dev.off()
