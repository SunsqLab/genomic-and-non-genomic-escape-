# S8: Logistic association of individual escape features with ICB response
# This script evaluates immune escape features or burden-based subgroups in ICB
# cohorts. It tests associations with treatment response or survival and saves the
# resulting statistical summaries and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(scatterpie); library(broom)
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]AssociationBetweenEscapeGroupsAndICBOutcome/Question2"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
ICB.ID.xteams=c(
    "Liu_NatureMedicine_2019", "VanAllen_Science_2015", "Gide_CancerCell_2019", "Hugo_Cell_2016", "Freeman_CellReportsMedicine_2022", "Abbott_ClinicalCancerResearch_2021",
    "Ravi_NatureGenetics_2023")
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
        tmp$B2M.biallelic.inactivation=Van_Hugo_binary[[x]][1, pos]
        tmp$B2M.biallelic.inactivation=ifelse(tmp$B2M.biallelic.inactivation==1, TRUE, FALSE)
    }
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
patient.center=CombineData.XTeam(ICB.ID.xteams, file.query="PatientCenter/PatientCenter.rds")
profile.ICB.clinic=lapply(ICB.ID.xteams, function(x){
    if(x=="Gide_CancerCell_2019"){
        tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("TreatmentResponse", "SampleTime", "TreatmentHistory", "CancerType"))
        tmp$TMB=NA
    }else{
        tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("TreatmentResponse", "SampleTime", "TMB", "TreatmentHistory", "CancerType"))   %>%
# Retain eligible samples and remove records that do not meet the analysis criteria.
            dplyr::filter(TreatmentHistory!="[NA(8)]")
    }
    if(x=="Liu_NatureMedicine_2019"){
        tmp$SampleTime="Pre-Treatment"
    }
    tmp=tmp[grepl("PRE", toupper(tmp$SampleTime)), ]
    tmp$log10.TMB=log10(tmp$TMB)
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

all.features.list=list(
    APMalt=c("HLA.LOH", "B2M.biallelic.inactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),
    Checkalt=c("CD274.deepAmp"),
    ActMalt=c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
    Checkexp=c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    Supprcell=c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    Supprsig=c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    Actdown=c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown=c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo=c("mean.neo.exprs")
)
all.features=unname(unlist(all.features.list))
names(all.features)=rep(names(all.features.list), lengths(all.features.list))

OR.RECIST=list(response=c("Complete response", "Partial response", "Very good partial response", "Response"), nonresponse=c("Stable disease", "Progressive disease", "Non-response"))

plot.data=lapply(names(profile.escape.binary), function(x){
    data=dplyr::full_join(profile.escape.binary[[x]], profile.ICB.clinic[[x]], by="SampleID")          %>%
        dplyr::filter(!is.na(SampleTime)) %>%
        dplyr::mutate(TreatmentResponse=case_when(
                    TreatmentResponse %in% OR.RECIST$response ~ 1,
                    TreatmentResponse %in% OR.RECIST$nonresponse ~ 0,
                    TRUE ~ NA),
                )    %>%
        dplyr::select(where(~ !all(is.na(.x))))
})      %>%     setNames(names(profile.escape.binary))

library(logistf)
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/ImmuneBiomarkerStat/Res.biomarker.Logist.R")
logistic.result=lapply(names(profile.escape.binary), function(x){
    tmp.data=plot.data[[x]]

    result=lapply(intersect(all.features, names(tmp.data)), function(feature){
        data=tmp.data[, c(feature, "TreatmentResponse", "TreatmentHistory", "CancerType")]        %>%     na.omit()
        data$sam.size=ifelse(nrow(data)>100, "large.size", "small.size")
        if(length(table(data[, c(feature, "TreatmentResponse")]))>3 & nrow(data)>29){
            data[, feature]=factor(data[, feature], levels=c("FALSE", "TRUE"))

            if(min(table(data[, c(feature, "TreatmentResponse")]))==0){
                fit=logistf(as.formula(paste("TreatmentResponse ~", feature)), data)
                data.frame(p.value=fit$prob[2], OR=fit$coefficients[2], OR.confint.lower=exp(fit$ci.lower[2]), OR.confint.upper=exp(fit$ci.upper[2]), feature, Data=x, CancerType=data$CancerType[1], sam.size=data$sam.size[1])
            }else{
                re=Res.biomarker.Logist(data, biomarker_col=feature, response_col="TreatmentResponse")
                data.frame(re, feature, Data=x, CancerType=data$CancerType[1], sam.size=data$sam.size[1])
            }
        }else{
            data.frame(p.value=NA, OR=NA, OR.confint.lower=NA, OR.confint.upper=NA, feature, Data=x, CancerType=data$CancerType[1], sam.size=data$sam.size[1])
        }
    })      %>%     do.call(what=rbind)
})      %>%     do.call(what=rbind)
logistic.result$OR.p_value=ifelse(logistic.result$OR>1, log10(logistic.result$p.value)*(-1), log10(logistic.result$p.value))
logistic.result$feature=factor(logistic.result$feature, levels=na.omit(unique(logistic.result$feature)[match(all.features, unique(logistic.result$feature))]))

source("/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethod/Script/ImmuneEscape/1.0.AssociationBetweenEscapePatternsAndICBOutcome/PlotGroupedLollipopChart.R")

data=logistic.result  %>%
    dplyr::mutate(
        CancerType=factor(CancerType, levels=c("Skin cancer", "Lung cancer", "Kidney cancer")),
        Data=factor(Data, levels=unique(Data[order(CancerType)])),
        feature.type=names(all.features)[match(feature, all.features)],
        feature.type=factor(feature.type, levels=names(all.features.list))
    )
show.Data=levels(data$Data)
data.shape=setNames(rep_len(c(1, 7, 8, 9, 15, 17, 21, 22, 23, 24), length(show.Data)), show.Data)
data.color=c("Skin cancer"="#ce3375", "Lung cancer"="#60c8b3")

p=Plot_BarPoint(data, x="feature", y="OR.p_value", groups="Data", fill.col="CancerType", color.col="CancerType",
        fill=data.color, col=data.color, shape=data.shape, y_label="-log10 (p.value)")+
        ggh4x::facet_grid2(. ~ feature.type, scales="free", space ="free_x", axes="x")+
        guides(fill=guide_legend(nrow=4), color=guide_legend(nrow=2), shape=guide_legend(nrow=3), alpha="none")
ggsave(file.path(Dir.output, paste0("1.3.responseandnonresponsesamplesescape_featuredifference[][logistic].pdf")), p, width=20, height=5)
