# S4: Multivariable logistic regression of escape features by cancer
# This script evaluates immune infiltration, TMB, and immune escape relationships
# across TCGA cancer types. It prepares cohort-level profiles, applies the specified
# association models, and saves the resulting tables and figures.

library(magrittr); library(ggplot2); library(cowplot); library(logistf); library(openxlsx)
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

# Danaher immune cell infiltration scores and CYT scores stored as a list.
profile.exprs=CombineData.XTeam(ID.xteams, file.query='OMICSData/Exprs.data.rds')
profile.exprs$CRC_TCGA=do.call(cbind, CombineData.XTeam(c('COAD_TCGA', 'READ_TCGA'), file.query='OMICSData/Exprs.data.rds'))
source('/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/ImmuneEscape/ImmuneEscapeMechanisms/ImmuneCellInfiltration/ImmuneCellInfiltrationScores.R')
profile.Danaher.cell=lapply(profile.exprs, function(x){
    tmp=ImmuneInfiltraScore(exprs.data=x)
})      %>%     setNames(ID.xteams)

# TMB and other clinical profiles stored as data frames.
profile.TMB=GetInfor.PatientCenter(patient.center, colNames=c("SampleID", "TMB", "Purity", "Ploidy"))

# Copy-number burden profiles stored as a list.
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
profile.CN=lapply(ID.xteams, function(x){
    tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("frac.cnv", "CancerType"))
    tmp$CancerType[tmp$CancerType%in%c("COAD", "READ")]="CRC"
    return(tmp)
})      %>%     setNames(ID.xteams)

# MSI subtype profiles stored as a list.
ID.xteams.MSI=c('CRC_TCGA', 'STAD_TCGA', 'UCEC_TCGA')
profile.MSI=lapply(ID.xteams.MSI, function(x){
    GetInfor.PatientCenter(patient.center[[x]], SampleID=profile.escape.binarization[[x]]$SampleID, colNames="MSISubtype2")
})      %>%     setNames(ID.xteams.MSI)

all.features.list=list(
    APMalt=c("HLA.LOH", "biallelic.B2M.inactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),
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

cell.types='Cytotoxic.cell'

# Perform multivariable logistic regression with LOH status as the outcome and TMB, immune infiltration, and CNV burden as predictors.
plots.data=lapply(ID.xteams, function(x){
    data=profile.escape.binarization[[x]]

    if(x%in%c("CRC_TCGA", "STAD_TCGA", "UCEC_TCGA")){
        data=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(data, profile.Danaher.cell[[x]], profile.CN[[x]], profile.TMB, profile.MSI[[x]]))      %>%
            dplyr::mutate(
                MSI.Subtype=ifelse(MSISubtype2=='MSI-H', 'MSI', ifelse(MSISubtype2 %in% c('MSS', 'MSI-L'), 'MSS', NA)))       %>%
            dplyr::filter(!is.na(MSI.Subtype))  %>%
            dplyr::mutate(MSI.Subtype=factor(MSI.Subtype, levels=c("MSS", "MSI")))
    }else{
        data=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(data, profile.Danaher.cell[[x]], profile.CN[[x]], profile.TMB))
    }
    data=data[!is.na(data$frac.cnv), ]

    tmp=lapply(intersect(all.features, names(data)), function(feature){
        if(sum(!is.na(data[, feature]))>0 && sum(!is.na(data$TMB))>0 && !min(table(data[, feature]))<10){

            data=data[!is.na(data$TMB) & !is.na(data[, cell.types]), ]
            data=Filter(function(col) !all(is.na(col)), data)

            data$log10.TMB=log10(data$TMB+1)
# Although logistf is designed to address complete separation, it may still fail to converge when variables are severely imbalanced or information is insufficient.
            tmp.cell=paste(cell.types, collapse="+")

            if(x%in%c("CRC_TCGA", "STAD_TCGA", "UCEC_TCGA")){
                formula=as.formula(paste0(feature, " ~ ", tmp.cell, " + log10.TMB + frac.cnv + Purity + MSI.Subtype"))
            }else{
                formula=as.formula(paste0(feature, " ~ ", tmp.cell, " + log10.TMB + frac.cnv + Purity"))
            }
            model=logistf::logistf(formula=formula, data=data, control=logistf.control(maxit=1000, maxstep=0.5))

            model.result=data.frame(coef=signif(model$coefficients[-1], 3), CI_2.5=model$ci.lower[-1], CI_97.5=model$ci.upper[-1], p_value=signif(model$prob[-1], 3))
            model.result$p_value=ifelse(model.result$p_value < 10^(-6), 10^(-6), model.result$p_value)
            result=data.frame(CancerType=gsub('_TCGA', '', x), variable=rownames(model.result), feature, model.result)
        }else{
            result=NULL
        }
        return(result)
    })
    tmp=do.call(rbind, tmp[lengths(tmp)!=0])
    tmp$FDR=p.adjust(tmp$p_value, method="BH")
    return(tmp)
})      %>%     do.call(what=rbind)
plots.data$variable[plots.data$variable=='log10.TMB']='TMB'
plots.data$variable[plots.data$variable=='frac.cnv']='CNA'
plots.data$variable[plots.data$variable=='MSI.SubtypeMSI']='MSI'

file_path=file.path(Dir.output, "Associations_between_immune_infiltration_and_escape_features.xlsx")
wb=loadWorkbook(file_path)
sheet.name='Table S3'
if (sheet.name %in% names(wb)) {
  removeWorksheet(wb, sheet.name)
}
addWorksheet(wb, sheet.name)
writeData(wb, sheet=sheet.name, x=paste0(sheet.name, ". Multivariate logistic regression analysis"), startRow=1, startCol=1)
writeData(wb, sheet=sheet.name, plots.data, startRow=2, headerStyle=createStyle(textDecoration="bold"))
saveWorkbook(wb, file=file_path, overwrite=TRUE)

colors=c("#2D3561", "#C05C7E", "#F3826F", "#FFB961", "#BCEAD5", "#124E96", "#0D8ABC", "#64A97B", "#a6761d", "#605EA1", "#967E76",
    "#b5e2fa", "#419197", "#7f4f24", "#D8B5DE", "#CEDF9F", "#678C40", "#FFAAAA", "#FF7777", "#AF1740", "#CDC1FF", "#3B1C32", "#6A1E55",
    "#003161", "#7ED4AD", "#8B5DFF", "#9694FF", "#1b9e77", "#1A1A19", "#9ceaef", "#CC2B52", "lightgrey")
cancer.color=setNames(colors, c(gsub("_TCGA", "", ID.xteams), "no.sig"))

pdf(file.path(Dir.output, "3.escape_feature_predictors_continuous_TMB_immune_infiltration_CNV_logistic_regression.pdf"), 5, 5)
plots=lapply(unique(plots.data$feature), function(feature){
    set.seed(8766)
    data=plots.data[plots.data$feature==feature & plots.data$variable!='(Intercept)', ]      %>%     na.omit()
    data$label=paste0(data$CancerType, "_", data$variable)
    data$label[!data$FDR<0.05]=""
    data$CancerType[data$FDR>0.05 | is.na(data$FDR)]="no.sig"
    data$CancerType=factor(data$CancerType, levels=c(unique(setdiff(data$CancerType, "no.sig")), "no.sig"))

    if(nrow(data)!=0){
        p=ggplot(data, aes(x=coef, y=log10(FDR)*(-1), label=label))+
            geom_hline(yintercept=(-1)*log10(0.05), linetype="dashed", col="gray80")+
            geom_vline(xintercept=0, linetype="dashed", col="gray80")+
            geom_point(aes(color=CancerType))+
            scale_color_manual(values=cancer.color)+
            labs(y='-log10 (FDR)', x='Regression coefficients', title=feature)+
            ggrepel::geom_text_repel(size=2.5, max.overlaps=Inf)+
            theme_bw()+
            theme(legend.position="bottom", axis.text = element_text(size=10), axis.title = element_text(size=10))
        p1=p+theme(legend.position='none')
        print(p)
        print(p1)
    }else{
        return(NULL)
    }
})
dev.off()
