# S7: Survival across escape subgroups by tumor purity
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
profile.surv=GetInfor.PatientCenter(patient.center, colNames=c("OS", "OS.time", "PFI", "PFI.time", "CancerType", "Purity"))     %>%
    dplyr::group_by(CancerType)      %>%
    dplyr::mutate(Purity.group=ifelse(Purity>median(Purity, na.rm=T), "high Purity", "low Purity"))     %>%
    as.data.frame()

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

plots.p.list=lapply(c("Genomic", "nonGenomic"), function(plus.features){
    print(plus.features)
    split.data=tmp.data[[plus.features]]
    plots.p=lapply(names(split.data), function(x){
        print(x)
        data=split(split.data[[x]], split.data[[x]]$Purity.group)

        result=lapply(data, function(df){
            if(length(unique(df$cluster))>1){
                p1=plot.surv(df,
                                    group=df$cluster,
                                    median.time=F,
                                    upper.time=365*5,
                                    main=paste0(x, '_5yr_OS_', df$Purity.group[1]),
                                    endpoint="OS",
                                    surv.median.line="hv",
                                    risk.table=TRUE, xlab="Time (days)",
                                    color=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7FFF", "#00A087FF", "#4DBBD5FF", "#3C5488FF", "#80b1d3", "#fb8072", "#bebada", "#8dd3c7", "#fccde5", "#bc80bd", "#d9d9d9", "#b3de69", "#fdb462", "#ffffb3")
                )
                p2=plot.surv(df,
                                    group=df$cluster,
                                    median.time=F,
                                    upper.time=365*5,
                                    main=paste0(x, '_5yr_PFI_', df$Purity.group[1]),
                                    endpoint="PFI",
                                    surv.median.line="hv",
                                    risk.table=TRUE, xlab="Time (days)",
                                    color=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7FFF", "#00A087FF", "#4DBBD5FF", "#3C5488FF", "#80b1d3", "#fb8072", "#bebada", "#8dd3c7", "#fccde5", "#bc80bd", "#d9d9d9", "#b3de69", "#fdb462", "#ffffb3")
                )
                p=list(p1, p2)
            }
        })
    })
})      %>%     setNames(c("Genomic", "nonGenomic"))

tmp=lapply(c("Genomic", "nonGenomic"), function(plus.features){
    tmp.dir=file.path(outDir, plus.features)
    if(!dir.exists(tmp.dir)){ dir.create(tmp.dir, recursive=TRUE) }

    pdf(file.path(tmp.dir, '7.survival_across_escape_subgroups_by_tumor_purity.5yr.pdf'), 8, 8)
    print(plots.p.list[[plus.features]])
    dev.off()
})
