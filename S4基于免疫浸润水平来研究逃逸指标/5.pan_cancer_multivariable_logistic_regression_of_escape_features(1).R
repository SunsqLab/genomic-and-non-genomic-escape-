# S4: Pan-cancer multivariable logistic regression of escape features
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
})      %>%     do.call(what=rbind)

# Danaher immune cell infiltration scores and CYT scores stored as a list.
profile.exprs=CombineData.XTeam(ID.xteams, file.query='OMICSData/Exprs.data.rds')
profile.exprs$CRC_TCGA=do.call(cbind, CombineData.XTeam(c('COAD_TCGA', 'READ_TCGA'), file.query='OMICSData/Exprs.data.rds'))
source('/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/ImmuneEscape/ImmuneEscapeMechanisms/ImmuneCellInfiltration/ImmuneCellInfiltrationScores.R')
profile.Danaher.cell=lapply(profile.exprs, function(x){
    tmp=ImmuneInfiltraScore(exprs.data=x)
})      %>%    do.call(what=rbind)

# TMB and other clinical profiles stored as data frames.
profile.TMB=GetInfor.PatientCenter(patient.center, colNames=c("SampleID", "TMB", "Purity", "Ploidy"))

# Copy-number burden profiles stored as a list.
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
profile.CN=lapply(ID.xteams, function(x){
    tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("frac.cnv", "CancerType"))
    tmp$CancerType[tmp$CancerType%in%c("COAD", "READ")]="CRC"
    return(tmp)
})      %>%     do.call(what=rbind)

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

# Perform multivariable logistic regression with LOH status as the outcome and TMB, immune infiltration, and CNV burden as predictors.
data=list(profile.escape.binarization, profile.Danaher.cell, profile.CN, profile.TMB)       %>%
    purrr::reduce(full_join, by="SampleID")   %>%
    dplyr::filter(!is.na(frac.cnv) & !is.na(TMB) & !is.na(Cytotoxic.cell)) %>%
    dplyr::mutate(log10.TMB=log10(TMB+1))       %>%
    select(where(~ !all(is.na(.x))))

pdf(file.path(Dir.output, paste0("4.1.escape_features_logistf_forest_with_cancer_type_as_fixed_covariate.pdf")), width=5, height=3)
plots=lapply(intersect(all.features, names(data)), function(feature){
    print(feature)
    if(sum(!is.na(data[, feature]))>0 && sum(!is.na(data$TMB))>0 && !min(table(data[, feature]))<10){

        formula=as.formula(paste0(feature, " ~ Cytotoxic.cell + log10.TMB + frac.cnv + Purity + CancerType"))
        model=logistf::logistf(formula=formula, data=data, control=logistf.control(maxit=1000, maxstep=0.5))

        logit_result=data.frame(term=names(model$coefficients), p.value=as.numeric(model$prob))      %>%
                dplyr::mutate(
                    OR=exp(as.numeric(model$coefficients)), CI.lower=exp(as.numeric(model$ci.lower)), CI.upper=exp(as.numeric(model$ci.upper)),
                    label=case_when(term == "Cytotoxic.cell" ~ "Cytotoxic cell", term == "log10.TMB" ~ "TMB", term == "frac.cnv" ~ "CNA", TRUE ~ term),
                    p.label=ifelse(p.value < 0.001, "p < 0.001", paste0("p=", signif(p.value, 2)))
                )        %>%
                dplyr::filter(term != "(Intercept)" & !grepl("^CancerType", term))

        p=ggplot(logit_result, aes(x=OR, y=reorder(label, OR))) +
                geom_vline(xintercept=1, linetype="dashed", linewidth=0.4) +
                geom_segment(aes(x=CI.lower, xend=CI.upper, y=label, yend=label), linewidth=0.6) +
                geom_point(size=2.5) +
                geom_text(aes(x=CI.upper * 1.2, label=p.label), hjust=0, size=3) +
                scale_x_log10() +
                labs(x="Odds ratio", y=NULL, title=feature) +
                theme_classic()
        print(p)
    }
})
dev.off()

# Use a mixed-effects model with cancer type as a random effect.
pdf(file.path(Dir.output, paste0("4.2.escape_features_mixed_effects_logistic_forest_with_cancer_type_as_random_effect.pdf")), width=5, height=3)
plots=lapply(intersect(all.features, names(data)), function(feature){
    print(feature)
    if(sum(!is.na(data[, feature]))>0 && sum(!is.na(data$TMB))>0 && !min(table(data[, feature]))<10){

        formula_mixed=as.formula(paste0(feature, " ~ Cytotoxic.cell + log10.TMB + frac.cnv + Purity + (1|CancerType)"))
        model_mixed=lme4::glmer(formula=formula_mixed, data=data, family=binomial(link="logit"), control=lme4::glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
        mixed_result=as.data.frame(summary(model_mixed)$coefficients)       %>%
            tibble::rownames_to_column("term") %>%
            dplyr::rename(beta=Estimate, SE=`Std. Error`, z.value=`z value`, p.value=`Pr(>|z|)`) %>%
            dplyr::mutate(OR=exp(beta), CI.lower=exp(beta - 1.96 * SE), CI.upper=exp(beta + 1.96 * SE),
                label=dplyr::case_when(term == "Cytotoxic.cell" ~ "Cytotoxic cell", term == "log10.TMB" ~ "TMB", term == "frac.cnv" ~ "CNA", TRUE ~ term),
                p.label=ifelse(p.value < 0.001, "p < 0.001", paste0("p=", signif(p.value, 2)))) %>%
            dplyr::filter(term != "(Intercept)")

        p=ggplot(mixed_result, aes(x=OR, y=reorder(label, OR))) +
            geom_vline(xintercept=1, linetype="dashed", linewidth=0.4) +
            geom_errorbarh(aes(xmin=CI.lower, xmax=CI.upper), height=0.2, linewidth=0.5) +
            geom_point(size=2.5) +
            geom_text(aes(x=CI.upper * 1.25, label=p.label), hjust=0, size=3) +
            scale_x_log10() +
            labs(x="Odds ratio", y=NULL, title=paste0("Mixed-effects logistic regression: ", feature)) +
            theme_classic()
        print(p)
    }
})
dev.off()
