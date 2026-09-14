# S8: Survival by genomic escape burden using multivariable Cox analysis
# This script evaluates immune escape features or burden-based subgroups in ICB
# cohorts. It tests associations with treatment response or survival and saves the
# resulting statistical summaries and figures.

library(magrittr); library(reshape2); library(readxl); library(circlize); library(cowplot); library(survival); library(survminer)
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]AssociationBetweenEscapeGroupsAndICBOutcome/Question2"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
ICB.ID.xteams=c(
    "Liu_NatureMedicine_2019", "VanAllen_Science_2015", "Hugo_Cell_2016", "Freeman_CellReportsMedicine_2022", "Abbott_ClinicalCancerResearch_2021",
    "Ravi_NatureGenetics_2023", "Alban_NatureMedicine_2024")
profile.mut=CombineData.XTeam(ICB.ID.xteams, file.query="OMICSData/Mutations.data.rds")
profile.CN=CombineData.XTeam(ICB.ID.xteams, file.query="OMICSData/CN.rds")
profile.CNMatrix=CombineData.XTeam(ICB.ID.xteams, file.query="OMICSData/CN.MatrixProfiler.rds")
common.sam.list=lapply(ICB.ID.xteams, function(x){
    if(!is.null(profile.CN[[x]])){
        tmp=Reduce(intersect, list(unique(profile.mut[[x]]$SampleID), unique(profile.CN[[x]]$SampleID)))
    }else if(!is.null(profile.CNMatrix[[x]])){
        tmp=Reduce(intersect, list(unique(profile.mut[[x]]$SampleID), colnames(profile.CNMatrix[[x]])))
    }else{
        tmp=NULL
    }
    return(tmp)
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
        tmp$biallelic.B2M.inactivation=Van_Hugo_binary[[x]][1, pos]
    }
    tmp$all.APM.mut=tmp.escape$APM.alt$all.APM.mut[match(tmp$SampleID, tmp.escape$APM.alt$SampleID)]
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
# Retain eligible samples and remove records that do not meet the analysis criteria.
patient.center=CombineData.XTeam(ICB.ID.xteams, file.query="PatientCenter/PatientCenter.rds")
profile.ICB.clinic=lapply(ICB.ID.xteams, function(x){
    if(all(c("Purity", "Ploidy")%in%colnames(patient.center[[x]]$SampleInfo))){
        tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("TreatmentResponse", "SampleTime", "TreatmentHistory", "CancerType", "OS", "OS.time", "TMB"))   %>%
        dplyr::filter(TreatmentHistory!="[NA(8)]")
        if(x=="Liu_NatureMedicine_2019"){
            tmp$SampleTime="Pre-Treatment"
        }
        tmp=tmp[grepl("PRE", toupper(tmp$SampleTime)), ]
        return(tmp)
    }else{
        return(NULL)
    }
})      %>%     setNames(ICB.ID.xteams)

genome.features=c(c("HLA.LOH", "biallelic.B2M.inactivation", "all.APM.mut"), c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"), "CD274.deepAmp")

escape.num.list=lapply(ICB.ID.xteams, function(x){
    tmp.data=Filter(function(x) !all(is.na(x)), profile.escape.binary[[x]])     %>%
        dplyr::filter(SampleID%in%common.sam.list[[x]])

    cluster.data=tmp.data[, intersect(colnames(tmp.data), c(genome.features, "SampleID"))]    %>%
        dplyr::mutate(Escape.burden=rowSums(dplyr::select(., -SampleID)))       %>%
        dplyr::mutate(
            es.feature.group=dplyr::case_when(
                0<Escape.burden & Escape.burden<3 ~ "1-2",
                Escape.burden>3 ~ "4+",
                TRUE ~ as.character(Escape.burden)),
            es.feature.group=paste0("E", es.feature.group),
            es.feature.group=factor(es.feature.group, levels=c("E0", "E1-2", "E3", "E4+")))      %>%
        dplyr::select(SampleID, Escape.burden, es.feature.group)
})      %>%     setNames(ICB.ID.xteams)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/KM.curve.and.logrank.test-2.0.R")
pdf(file.path(Dir.output, "2.5.all_genomicescape_featureburden_groupsamplessurvival_difference.pdf"), 6, 7)
surv.data=lapply(names(escape.num.list), function(x){
    if(!is.null(profile.ICB.clinic[[x]])){
        tmp.data=dplyr::full_join(escape.num.list[[x]], profile.ICB.clinic[[x]], by="SampleID")   %>%
            dplyr::select(SampleTime, OS.time, OS, es.feature.group, PatientID)      %>%
            na.omit()      %>%
            dplyr::mutate(
                es.feature.group=factor(es.feature.group),
                OS=ifelse(OS.time > 365*5, 0, OS),
                OS.time=ifelse(OS.time > 365*5, 365*5, OS.time)
            )

        if(length(table(tmp.data$es.feature.group))>1 & min(table(tmp.data$es.feature.group))!=0){
            p=plot.surv(tmp.data,
                        group=tmp.data$es.feature.group,
                        median.time=F,
                        main=x,
                        endpoint="OS",
                        surv.median.line="hv",
                        risk.table=TRUE, xlab="Time (days)",
                        color=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7FFF", "#00A087FF", "#4DBBD5FF", "#3C5488FF", "#80b1d3", "#fb8072", "#bebada", "#8dd3c7", "#fccde5", "#bc80bd", "#d9d9d9", "#b3de69", "#fdb462", "#ffffb3")
            )
            print(p)
        }
    }
})
dev.off()

tmp.dir=file.path(Dir.output, "Genomic")
if(!dir.exists(tmp.dir)){ dir.create(tmp.dir, recursive=TRUE) }
pdf(file.path(tmp.dir, '3.escape_groupmultivariableCoxanalysis[based_onbinarizedescape_featurecumulative_count_combinedtranscriptomic][5yr].pdf'), 6, 4)
ICB.Multivariable.Cox.result=lapply(c("Liu_NatureMedicine_2019", "Ravi_NatureGenetics_2023"), function(x){
    if(!is.null(profile.ICB.clinic[[x]])){
        print(x)
        df=dplyr::full_join(escape.num.list[[x]], profile.ICB.clinic[[x]], by="SampleID")   %>%
            dplyr::select(SampleTime, OS.time, OS, Escape.burden, es.feature.group, PatientID, TMB)      %>%
            na.omit()      %>%
            dplyr::mutate(
                cluster=factor(es.feature.group),
                log10TMB=log10(TMB+1),
                OS=ifelse(OS.time > 365*5, 0, OS),
                OS.time=ifelse(OS.time > 365*5, 365*5, OS.time)
            )

        if(length(table(df$es.feature.group))>1 & min(table(df$es.feature.group))!=0){
            os_form2=coxph(as.formula("Surv(OS.time, OS) ~ cluster + log10TMB"), data=df)
            p2=ggforest(os_form2, data=df, main=paste0(x, "_Hazard ratio"), cpositions=c(0.02,-0.15, 0.25), fontsize=0.8, noDigits=2)
            print(p2)
        }
    }
})
dev.off()

tmp.dir=file.path(Dir.output, "Genomic")
if(!dir.exists(tmp.dir)){ dir.create(tmp.dir, recursive=TRUE) }
pdf(file.path(tmp.dir, '3.escape_groupmultivariableCoxanalysis[by_datasetTMBmedian][based_onbinarizedescape_featurecumulative_count_combinedtranscriptomic][5yr].pdf'), 6, 4)
ICB.Multivariable.Cox.result=lapply(names(escape.num.list), function(x){
    if(!is.null(profile.ICB.clinic[[x]])){
        print(x)
        df=dplyr::full_join(escape.num.list[[x]], profile.ICB.clinic[[x]], by="SampleID")   %>%
            dplyr::select(SampleTime, OS.time, OS, Escape.burden, es.feature.group, PatientID, TMB)      %>%
            na.omit()      %>%
            dplyr::mutate(
                cluster=factor(es.feature.group),
                log10TMB=log10(TMB+1),
                OS=ifelse(OS.time > 365*5, 0, OS),
                OS.time=ifelse(OS.time > 365*5, 365*5, OS.time)
            )

        cox.data.list=list(
            All=df
        )

        cox.result=lapply(names(cox.data.list), function(tmb.group){
            tmp.df=cox.data.list[[tmb.group]]     %>%
                dplyr::mutate(cluster=droplevels(cluster))
            print(paste0(x, " | ", tmb.group, " | n=", nrow(tmp.df), " | event=", sum(tmp.df$OS)))

            if(nrow(tmp.df)>0 &
               length(unique(tmp.df$cluster))>1 &
               length(unique(tmp.df$OS))>1 &
               sum(tmp.df$OS)>0){
                os_form2=tryCatch(
                    coxph(as.formula("Surv(OS.time, OS) ~ cluster + log10TMB"), data=tmp.df),
                    error=function(e){
                        print(paste0("Skip ", x, " | ", tmb.group, ": coxph failed - ", e$message))
                        return(NULL)
                    }
                )
                if(is.null(os_form2)){
                    return(NULL)
                }

                cox.sum=summary(os_form2)

                p2=tryCatch(
                    ggforest(
                        os_form2,
                        data=tmp.df,
                        main=paste0(x, "_", tmb.group, "_Hazard ratio"),
                        cpositions=c(0.02,-0.15, 0.25),
                        fontsize=0.8,
                        noDigits=2
                    ),
                    error=function(e){
                        print(paste0("Skip plot ", x, " | ", tmb.group, ": ggforest failed - ", e$message))
                        return(NULL)
                    }
                )
                if(!is.null(p2)){
                    print(p2)
                }
                return(os_form2)
            }else{
                print(paste0("Skip ", x, " | ", tmb.group, ": insufficient samples/events or only one escape group."))
                return(NULL)
            }
        })      %>%     setNames(names(cox.data.list))
        return(cox.result)
    }
})
dev.off()
