# S7: Pan-cancer multivariable Cox analysis adjusted for cancer type
# This script defines immune escape subgroups and evaluates their molecular or
# clinical characteristics across TCGA cohorts. It performs the indicated subgroup
# comparisons or survival models and saves the resulting figures.

library(magrittr); library(ggplot2); library(readxl); library(writexl); library(survival); library(openxlsx); library(survminer)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Multiple_immune_escape_subgroup_patterns/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

outDir=file.path(Dir.output, 'All_escape_features_binary_distance')
profile.msd.cluster=readRDS(file.path(outDir, 'Es.features.sum.combinedTrans.20260510.cluster.rds'))
profile.median.cluster=readRDS(file.path(outDir, 'Es.features.median.Trans.20260610.cluster.rds'))
profile.2sd.cluster=readRDS(file.path(outDir, 'Es.features.2sd.Trans.20260610.cluster.rds'))

source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.surv=GetInfor.PatientCenter(patient.center, colNames=c("OS", "OS.time", "PFI", "PFI.time", "CancerType", "Purity", "Ploidy"))

profile.TMB=GetInfor.PatientCenter(patient.center, colNames=c("SampleID", "TMB"))

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/KM.curve.and.logrank.test-2.0.R")
if(!exists("profile.msd.cluster")){
    profile.msd.cluster=readRDS(file.path(outDir, 'Es.features.sum.combinedTrans.20260510.cluster.rds'))
}
profile.cluster=list(median=profile.median.cluster, m2sd=profile.2sd.cluster, msd=profile.msd.cluster)
tmp.data.list=lapply(profile.cluster, function(cluster.data){
    tmp.data=lapply(names(cluster.data), function(plus.features){
        escape.cluster.data=cluster.data[[plus.features]]
        if(is.list(escape.cluster.data) && !is.data.frame(escape.cluster.data)){
            escape.cluster.data=dplyr::bind_rows(escape.cluster.data)
        }
        tmp=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(escape.cluster.data, profile.surv, profile.TMB))   %>%
            na.omit()       %>%
            dplyr::mutate(
                OS=ifelse(OS.time > 365*5, 0, OS),
                OS.time=ifelse(OS.time > 365*5, 365*5, OS.time),
                PFI=ifelse(PFI.time > 365*5, 0, PFI),
                PFI.time=ifelse(PFI.time > 365*5, 365*5, PFI.time)
            )

        split.data=list(
            PanCancer_all=tmp,
            PanCancer_TMB_low_ltCancerMedian=dplyr::group_by(tmp, CancerType) %>%
                dplyr::filter(TMB < median(TMB, na.rm=TRUE)) %>%
                dplyr::ungroup(),
            PanCancer_TMB_high_geCancerMedian=dplyr::group_by(tmp, CancerType) %>%
                dplyr::filter(TMB >= median(TMB, na.rm=TRUE)) %>%
                dplyr::ungroup()
        )

        return(split.data)
    })          %>%     setNames(names(cluster.data))
})      %>%     setNames(names(profile.cluster))

Multivariable.Cox.result.list=lapply(tmp.data.list, function(tmp.data){
    Multivariable.Cox.result=lapply(names(tmp.data), function(plus.features){
        split.data=tmp.data[[plus.features]]

        cancer.result=lapply(names(split.data), function(x){
            data=split.data[[x]]
            df=na.omit(data[, c('cluster', 'OS', 'OS.time', 'PFI', 'PFI.time', 'CancerType')])
            df$cluster=relevel(factor(df$cluster, levels=c("E0", "E1-2", "E3", "E4+")), ref="E0")

                os_form=coxph(as.formula(paste0("Surv(OS.time, OS) ~ cluster + strata(CancerType)")), data=df)
                pfi_form=coxph(as.formula(paste0("Surv(PFI.time, PFI) ~ cluster + strata(CancerType)")), data=df)
                return(list(os_form, pfi_form))
        })      %>%     setNames(names(split.data))
    })      %>%     setNames(names(tmp.data))
})      %>%     setNames(names(tmp.data.list))

source("/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/Script/ImmuneEscape/0.9.ImmuneEscapeSubgroups/Function.SurvivalAnalysisCoxRegression.R")
Multivariable.Cox.plot.list=lapply(names(tmp.data.list), function(cutoff){
    print(cutoff)
    tmp.data=tmp.data.list[[cutoff]]
    Multivariable.Cox.plot=lapply(names(tmp.data), function(plus.features){
        print(plus.features)

        split.data=tmp.data[[plus.features]]
        tmp.dir=file.path(outDir, plus.features, cutoff)
        if(!dir.exists(tmp.dir)){ dir.create(tmp.dir, recursive=TRUE) }

        pdf(file.path(tmp.dir, '5.2.pan_cancer_multivariable_cox_adjusted_for_cancer_type.5yr.pdf'), 7, 4)
        cancer.result=lapply(names(split.data), function(x){
            print(x)
            data=split.data[[x]]
            surv.result=Multivariable.Cox.result.list[[cutoff]][[plus.features]][[x]]
            df=na.omit(data[, c('cluster', 'OS', 'OS.time', 'PFI', 'PFI.time', 'CancerType')])
            df$cluster=relevel(factor(df$cluster, levels=c("E0", "E1-2", "E3", "E4+")), ref="E0")
            if(all(table(df$OS)>20)){

                os_form=surv.result[[1]]
                pfi_form=surv.result[[2]]

                p=plot_multivariable_cox_forest(os_form, main=paste0(x, "_Cox models stratified by cancer type (OS); escape cluster ref: E0; adjusted by CancerType strata"))
            }else{
                p=NULL
            }
            print(p)
        })
        dev.off()
    })
})

Multivariable.Cox.excel.df=lapply(names(tmp.data.list), function(cutoff){
    tmp.data=tmp.data.list[[cutoff]]
    Multivariable.Cox.excel=lapply(names(tmp.data), function(plus.features){

        surv.result=Multivariable.Cox.result.list[[cutoff]][[plus.features]]
        cancer.result=lapply(names(surv.result), function(x){
            os_form=summary(surv.result[[x]][[1]])
            pfi_form=summary(surv.result[[x]][[2]])

            os_result=data.frame(
                CancerType=x, Escape.group=plus.features, cutoff=cutoff,
                Surv.group="OS",
                Variable=rownames(os_form$coefficients), HR=os_form$coefficients[, "exp(coef)"],
                lower95=os_form$conf.int[, "lower .95"], upper95=os_form$conf.int[, "upper .95"],
                C.index=os_form$concordance[1], C.index.se=os_form$concordance[2],
                p.value=os_form$coefficients[, "Pr(>|z|)"])

            pfi_result=data.frame(
                CancerType=x, Escape.group=plus.features, cutoff=cutoff,
                Surv.group="PFI",
                Variable=rownames(pfi_form$coefficients), HR=pfi_form$coefficients[, "exp(coef)"],
                lower95=pfi_form$conf.int[, "lower .95"], upper95=pfi_form$conf.int[, "upper .95"],
                C.index=pfi_form$concordance[1], C.index.se=pfi_form$concordance[2],
                p.value=pfi_form$coefficients[, "Pr(>|z|)"])
            result=rbind(os_result, pfi_result)
        })          %>%     do.call(what=rbind)
    })      %>%     do.call(what=rbind)
})      %>%      do.call(what=rbind)

PanCancer.Cindex.plot.df=Multivariable.Cox.excel.df %>%
    dplyr::filter(CancerType=="PanCancer_all") %>%
    dplyr::select(Escape.group, cutoff, Surv.group, C.index) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Escape.group=factor(Escape.group, levels=unique(Escape.group)),
        cutoff=factor(cutoff, levels=unique(cutoff)))

p=ggplot(PanCancer.Cindex.plot.df, aes(x=Escape.group, y=C.index, color=cutoff, group=cutoff)) +
    geom_line(linewidth=0.8) +
    geom_point(size=2) +
    scale_color_manual(values=c("msd"="#EC0404", "median"="#42B541", "m2sd"="#0499B4")) +
    facet_wrap(~Surv.group, nrow=1) +
    ylim(c(0.4, 0.6))+
    labs(title="PanCancer multivariable Cox C-index adjusted by CancerType", x="", y="C-index", color="binarization.cutoff") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid.minor=element_blank())
ggsave(filename=file.path(outDir, "6.2.pan_cancer_multivariable_cox_results_adjusted_for_cancer_type.pdf"), p, width=10, height=4)
