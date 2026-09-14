# S7: Survival across escape subgroups
# This script defines immune escape subgroups and evaluates their molecular or
# clinical characteristics across TCGA cohorts. It performs the indicated subgroup
# comparisons or survival models and saves the resulting figures.

library(magrittr); library(ggplot2); library(readxl); library(writexl); library(survival); library(openxlsx)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Multiple_immune_escape_subgroup_patterns/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

outDir=file.path(Dir.output, 'All_escape_features_binary_distance')
profile.cluster=readRDS(file.path(outDir, 'Es.features.sum.combinedTrans.20260510.cluster.rds'))

source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.surv=GetInfor.PatientCenter(patient.center, colNames=c("OS", "OS.time", "PFI", "PFI.time", "CancerType"))

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/KM.curve.and.logrank.test-2.0.R")
tmp.data=lapply(names(profile.cluster), function(plus.features){
    tmp=dplyr::full_join(profile.cluster[[plus.features]], profile.surv, by='SampleID')   %>%
        na.omit()       %>%
        dplyr::mutate(
            OS=ifelse(OS.time > 365*5, 0, OS),
            OS.time=ifelse(OS.time > 365*5, 365*5, OS.time),
            PFI=ifelse(PFI.time > 365*5, 0, PFI),
            PFI.time=ifelse(PFI.time > 365*5, 365*5, PFI.time)
        )

    split.data=split(tmp, tmp$CancerType)
    split.data[c('MESO', 'PCPG', 'TGCT', 'THCA')]=NULL
    split.data$PanCancer=tmp      %>%
        dplyr::mutate(CancerType="PanCancer")

    return(split.data)
})          %>%     setNames(names(profile.cluster))

plots.p.list=lapply(names(tmp.data), function(plus.features){
    print(plus.features)
    split.data=tmp.data[[plus.features]]
    plots.p=lapply(names(split.data), function(x){
        print(x)
        if(length(unique(split.data[[x]]$cluster))>1){
            p1=plot.surv(split.data[[x]],
                                group=split.data[[x]]$cluster,
                                median.time=F,
                                upper.time=365*5,
                                main=paste0(x, '_5yr_OS'),
                                endpoint="OS",
                                surv.median.line="hv",
                                risk.table=TRUE, xlab="Time (days)",
                                color=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7FFF", "#00A087FF", "#4DBBD5FF", "#3C5488FF", "#80b1d3", "#fb8072", "#bebada", "#8dd3c7", "#fccde5", "#bc80bd", "#d9d9d9", "#b3de69", "#fdb462", "#ffffb3")
            )
            p2=plot.surv(split.data[[x]],
                                group=split.data[[x]]$cluster,
                                median.time=F,
                                upper.time=365*5,
                                main=paste0(x, '_5yr_PFI'),
                                endpoint="PFI",
                                surv.median.line="hv",
                                risk.table=TRUE, xlab="Time (days)",
                                color=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7FFF", "#00A087FF", "#4DBBD5FF", "#3C5488FF", "#80b1d3", "#fb8072", "#bebada", "#8dd3c7", "#fccde5", "#bc80bd", "#d9d9d9", "#b3de69", "#fdb462", "#ffffb3")
            )
            p=list(p1, p2)
        }
    })
})      %>%     setNames(names(tmp.data))

tmp=lapply(names(tmp.data), function(plus.features){
    tmp.dir=file.path(outDir, plus.features)
    if(!dir.exists(tmp.dir)){ dir.create(tmp.dir, recursive=TRUE) }

    pdf(file.path(tmp.dir, '3.survival_across_escape_subgroups.5yr.pdf'), 8, 8)
    print(plots.p.list[[plus.features]])
    dev.off()
})

Univariate.Cox.result=lapply(names(tmp.data), function(plus.features){
    print(plus.features)
    split.data=tmp.data[[plus.features]]
    x="PanCancer"

    os_obj=Surv(time=split.data[[x]]$OS.time, event=split.data[[x]]$OS)
    pfi_obj=Surv(time=split.data[[x]]$PFI.time, event=split.data[[x]]$PFI)

    os_cox_model=coxph(os_obj ~ Escape.burden, data=split.data[[x]])
    pfi_cox_model=coxph(pfi_obj ~ Escape.burden, data=split.data[[x]])

    os_cox=summary(os_cox_model)
    pfi_cox=summary(pfi_cox_model)

    os_result=data.frame(
        Escape.group=plus.features, Surv.group="OS", Sample.num=nrow(split.data[[x]]),
        Variable = rownames(os_cox$coefficients), HR = os_cox$coefficients[, "exp(coef)"],
        lower95 = os_cox$conf.int[, "lower .95"], upper95 = os_cox$conf.int[, "upper .95"],
        p.value = os_cox$coefficients[, "Pr(>|z|)"])

    pfi_result=data.frame(
        Escape.group=plus.features, Surv.group="PFI", Sample.num=nrow(split.data[[x]]),
        Variable = rownames(pfi_cox$coefficients), HR = pfi_cox$coefficients[, "exp(coef)"],
        lower95 = pfi_cox$conf.int[, "lower .95"], upper95 = pfi_cox$conf.int[, "upper .95"],
        p.value = pfi_cox$coefficients[, "Pr(>|z|)"])
    surv.result=rbind(os_result, pfi_result)

    return(surv.result)
})      %>%     setNames(names(tmp.data))
result=do.call(rbind, Univariate.Cox.result)
