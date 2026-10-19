# S8: ORR association of individual escape features with ICB response
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
        tmp$B2M.bi.inactive=Van_Hugo_binary[[x]][1, pos]
        tmp$B2M.bi.inactive=ifelse(tmp$B2M.bi.inactive==1, TRUE, FALSE)
    }
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
patient.center=CombineData.XTeam(ICB.ID.xteams, file.query="PatientCenter/PatientCenter.rds")
profile.ICB.clinic=lapply(ICB.ID.xteams, function(x){
    tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("TreatmentResponse", "SampleTime", "TreatmentHistory", "CancerType"))   %>%
# Retain eligible samples and remove records that do not meet the analysis criteria.
        dplyr::filter(TreatmentHistory!="[NA(8)]")
    if(x=="Liu_NatureMedicine_2019"){
        tmp$SampleTime="Pre-Treatment"
    }
    tmp=tmp[grepl("PRE", toupper(tmp$SampleTime)), ]
    return(tmp)
})  %>%     setNames(ICB.ID.xteams)

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
all.features=unname(unlist(all.features.list))
names(all.features)=rep(names(all.features.list), lengths(all.features.list))
all.features.color=rep(c("#c82621", "#F6C141", "#fa8b69ff", "#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"), lengths(all.features.list))    %>%
    setNames(unlist(all.features.list))

OR.RECIST=list(response=c("Complete response", "Partial response", "Very good partial response", "Response"), nonresponse=c("Stable disease", "Progressive disease", "Non-response"))

plot.data=lapply(names(profile.escape.binary), function(x){
    data=dplyr::full_join(profile.escape.binary[[x]], profile.ICB.clinic[[x]], by="SampleID")          %>%
        dplyr::filter(!is.na(SampleTime)) %>%
        dplyr::mutate(TreatmentResponse=case_when(
                    TreatmentResponse %in% OR.RECIST$response ~ "R",
                    TreatmentResponse %in% OR.RECIST$nonresponse ~ "NR",
                    TRUE ~ NA)) %>%
        dplyr::select(where(~ !all(is.na(.x))))
})      %>%     setNames(names(profile.escape.binary))

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/BioStat/AssociationBetweenTwoCategoricalVariables.FisherExactTest.R")
ORR.result=lapply(names(profile.escape.binary), function(x){
    tmp.data=plot.data[[x]]

    tmp=lapply(intersect(all.features, names(tmp.data)), function(feature){
        data=tmp.data[, c(feature, "TreatmentResponse", "TreatmentHistory", "CancerType")]        %>%             na.omit()
        data$sam.size=ifelse(nrow(data)>100, "large.size", "small.size")

        if(nrow(data)>29 & length(table(data[c(feature, "TreatmentResponse")]))==4){
            re=BinaryVarFisher(data, variable1=feature, variable2="TreatmentResponse")
            data.frame(re, feature=feature, Data=x, CancerType=data$CancerType[1], sam.size=data$sam.size[1])
        }else{
            NULL
        }
    })      %>%     do.call(what=rbind)
})        %>%     do.call(what=rbind)     %>%       na.omit()
ORR.result$feature=factor(ORR.result$feature, levels=intersect(all.features, unique(ORR.result$feature)))
ORR.result$direction=ifelse(ORR.result$odds_ratio<1, "odds.ratio<1", ifelse(ORR.result$odds_ratio==1, "odds.ratio=1", "odds.ratio>1"))

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/DotPlot.R")
data=ORR.result %>%
    dplyr::mutate(significance=case_when(p_value<0.01 ~ "p<0.01", p_value<0.05 ~ "p<0.05", TRUE ~ "no.sig"),
        CancerType=factor(CancerType, levels=c("Skin cancer", "Lung cancer")),
        Data=factor(Data, levels=unique(Data[rev(order(CancerType))]))
        )

p=dot_plot(t.data.frame=data, x="feature", y="Data", size="significance", size.value=c("p<0.01"=4, "p<0.05"=3, "no.sig"=0.5),
        color="direction", color.value=c("odds.ratio<1"="#4393C3", "odds.ratio=1"="#D1E5F0", "odds.ratio>1"="#D6604D"),
        facet="rows") +
        theme(axis.text.x=element_text(size=10, color=all.features.color[levels(data$feature)], hjust=1, angle=45), axis.text.y=element_text(size=10))
ggsave(file.path(Dir.output, paste0("1.1.escape_featureandandICBobjective_response_rateassociation.pdf")), p, width=13, height=4)

source("/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethod/Script/ImmuneEscape/1.0.AssociationBetweenEscapePatternsAndICBOutcome/PlotObjectiveResponseRate.BarAndPieCharts.R")
data=ORR.result       %>%
    dplyr::mutate(Data=factor(Data, levels=unique(Data)))
data$TRUE.ORR=data$TRUE.R / (data$TRUE.R + data$TRUE.NR)
data$FALSE.ORR=data$FALSE.R / (data$FALSE.R + data$FALSE.NR)

plots=Plot_ORR_Bar(data=data, x="feature", y="Data",
        left.col=c("TRUE.R", "TRUE.NR"), right.col=c("FALSE.R", "FALSE.NR"),
        left.ORR="TRUE.ORR", right.ORR="FALSE.ORR", p.row.num=3, facet.color=all.features.color)

pdf(file.path(Dir.output, paste0("1.2.escape_featuretwo_groupssamplesICBresponse_rate.pdf")), width=30, height=10)
result=lapply(plots, function(p)    {   print(p)    })
dev.off()
