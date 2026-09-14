# S4: Effect of TMB on infiltration–escape associations by cancer
# This script evaluates immune infiltration, TMB, and immune escape relationships
# across TCGA cancer types. It prepares cohort-level profiles, applies the specified
# association models, and saves the resulting tables and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(scatterpie); library(Cairo); library(openxlsx)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_infiltration_affect_escape_feature_prevalence/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Binarized immune escape profiles stored as a list.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds')
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
# Retain primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData/DataCenter/TCGA/PanCancer_TCGA/PatientCenter/PatientCenter.rds')
profile.escape.binarization=lapply(ID.xteams, function(x){
            data=profile.escape.binarization[[x]]
			tmp=lapply(data, function(y){
						y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
					})
            tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
            tmp$biallelic.B2M.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
            tmp$biallelic.B2M.inactivation=ifelse(tmp$biallelic.B2M.inactivation==1, TRUE, FALSE)
            return(tmp)
		})      %>%     setNames(ID.xteams)

# Danaher immune cell infiltration scores stored as a list.
profile.exprs=CombineData.XTeam(ID.xteams, file.query='OMICSData/Exprs.data.rds')
profile.exprs$CRC_TCGA=do.call(cbind, CombineData.XTeam(c('COAD_TCGA', 'READ_TCGA'), file.query='OMICSData/Exprs.data.rds'))
source('/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/ImmuneEscape/ImmuneEscapeMechanisms/ImmuneCellInfiltration/ImmuneCellInfiltrationScores.R')
profile.Danaher.cell=lapply(profile.exprs, function(x){
    tmp=ImmuneInfiltraScore(exprs.data=x)
})      %>%     setNames(ID.xteams)

# TMB profiles stored as data frames.
profile.TMB=GetInfor.PatientCenter(patient.center, colNames=c("SampleID", "TMB"))

cell.types=c('Cytotoxic.cell', 'CD8.Tcell', 'NK.cell', 'DC', 'CD4.Tcell', 'Bcell', 'Neutrophils', 'Mast.cell')

all.features.list1=list(
    APMalt=c("HLA.LOH", "biallelic.B2M.inactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),
    Checkalt=c("CD274.deepAmp"),
    ActMalt=c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut")
)
all.features1=unname(unlist(all.features.list1))
names(all.features1)=rep(c("#c82621", "#F6C141", "#fa8b69ff"), lengths(all.features.list1))

all.features.list2=list(
    Checkexp=c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    Supprcell=c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    Supprsig=c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    Actdown=c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown=c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo=c("mean.neo.exprs")
)
all.features2=unname(unlist(all.features.list2))
names(all.features2)=rep(c("#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"), lengths(all.features.list2))

# Compare immune escape feature prevalence between high- and low-infiltration groups using the Wilcoxon test.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/ImmuneBiomarkerStat/Res.biomarker.Wilcox.R")
cell='Cytotoxic.cell'
fisher.result=lapply(ID.xteams, function(x){
    data=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(profile.escape.binarization[[x]], profile.Danaher.cell[[x]]), profile.TMB)

    tmp.result=lapply(c(all.features1, all.features2), function(feature){
        df=na.omit(data[, c(cell, 'TMB', feature)])
        df$TMB=ifelse(df$TMB<10, "lowTMB", "highTMB")
        lapply(c("highTMB", "lowTMB"), function(tmp.TMB){
            tmp.df=df[df$TMB==tmp.TMB, c(cell, feature)]
            tmp.df[, feature]=factor(tmp.df[, feature], levels=c("TRUE", "FALSE"))
            tmp=Res.biomarker.Wilcox(tmp.df, binary_col=feature, continuous_col=cell)
            if(!is.null(tmp)){
                data.frame(CancerType=gsub('_TCGA', '', x), cell=cell, escape.feature=feature, TMBstatus=tmp.TMB, tmp)
            }else{
                data.frame(CancerType=gsub('_TCGA', '', x), cell=cell, escape.feature=feature, TMBstatus=tmp.TMB, p_value=NA, fold_change=NA, TRUE.median.score=NA, FALSE.median.score=NA)
            }
        })      %>%         do.call(what=rbind)
    })      %>%         do.call(what=rbind)
    tmp.result$FDR=p.adjust(tmp.result$p_value, method="BH")
    return(tmp.result)
})      %>%         do.call(what=rbind)
fisher.result$significance=ifelse(fisher.result$FDR<0.01, 'FDR<0.01', ifelse(fisher.result$FDR<0.05, 'FDR<0.05', 'no.sig'))
fisher.result$direction=ifelse(fisher.result$fold_change>1, "positive", ifelse(fisher.result$fold_change<1, "negative", "other"))

fisher.result=fisher.result[fisher.result$CancerType%in%c("CRC", "STAD", "UCEC"), ]
fisher.result$CancerType=factor(fisher.result$CancerType, levels=c("CRC", "STAD", "UCEC"))
shade_df=data.frame(
    CancerType=c("CRC", "STAD", "UCEC"),
    fill=c("#ffe5d9", "#ffcad4", "#f4acb7")
)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/CircleDotMatrixPlot.R")
pdf(file.path(Dir.output, "2.3.TMB_effect_on_genomic_and_non_genomic_escape_infiltration_associations.circle_plot.pdf"), width=7, height=7)
plots=lapply(list(all.features1, all.features2), function(feature){
    tmp.data=fisher.result[fisher.result$escape.feature%in%feature, ]
    tmp.data$escape.feature=factor(tmp.data$escape.feature, levels=feature)

    split.data=split(tmp.data, tmp.data$TMBstatus)
    plots.p=lapply(split.data, function(data){
        p=Dot2Circle(data, x.col="escape.feature", y.col="CancerType",
            color.col="direction", size.col="significance",
            color.values=c("positive"="#e26d5c", "negative"="#6096ba", "other"="lightgrey"),
            size.values=c("FDR<0.01"=3.5, "FDR<0.05"=3, "no.sig"=1), plot.title=data$TMBstatus[1], shade_df)+
            theme(axis.text.x=element_text(color=names(feature)[match(levels(data$escape.feature), feature)]))
        print(p)
    })
})
dev.off()
