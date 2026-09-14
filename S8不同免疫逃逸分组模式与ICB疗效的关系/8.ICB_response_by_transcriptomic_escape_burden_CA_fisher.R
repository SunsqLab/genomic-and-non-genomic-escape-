# S8: ICB response by transcriptomic escape burden using CA and Fisher tests
# This script evaluates immune escape features or burden-based subgroups in ICB
# cohorts. It tests associations with treatment response or survival and saves the
# resulting statistical summaries and figures.

library(magrittr); library(reshape2); library(readxl); library(circlize); library(cowplot); library("survival"); library("survminer")
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]AssociationBetweenEscapeGroupsAndICBOutcome/Question2"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
ICB.ID.xteams=c(
    "Liu_NatureMedicine_2019", "VanAllen_Science_2015", "Gide_CancerCell_2019", "Hugo_Cell_2016", "Freeman_CellReportsMedicine_2022", "Abbott_ClinicalCancerResearch_2021",
    "Ravi_NatureGenetics_2023", "Alban_NatureMedicine_2024")
profile.exprs=CombineData.XTeam(ICB.ID.xteams, file.query="OMICSData/Exprs.data.rds")
ICB.ID.xteams=ICB.ID.xteams[lengths(profile.exprs)!=0]
common.sam.list=lapply(ICB.ID.xteams, function(x){
    tmp=colnames(profile.exprs[[x]])
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
        tmp$B2M.bi.inactive=Van_Hugo_binary[[x]][1, pos]
    }
    tmp$all.APM.mut=tmp.escape$APM.alt$all.APM.mut[match(tmp$SampleID, tmp.escape$APM.alt$SampleID)]
    return(tmp)
})      %>%     setNames(ICB.ID.xteams)

# Apply the stated eligibility and data-quality restrictions before analysis.
profile.trans.combine.escape=lapply(ICB.ID.xteams, function(x){
    tmp=profile.escape.binary[[x]]        	%>%
        dplyr::group_by(SampleID)       %>%
        dplyr::mutate(
            AntigenPresentGene.downregulation=any(HLA.A, HLA.B, HLA.C, HLA.score, CALR, na.rm=TRUE),
            Checkpoint.overexprs=any(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm=TRUE),
            ImmunoSuppressiveCell.overexprs=any(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm=TRUE),
            ImmunoSuppressiveSig.overexprs=any(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm=TRUE),
            ImmuneActivationGene.downregulation=any(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm=TRUE))         %>%
        as.data.frame()
})      %>%     setNames(ICB.ID.xteams)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
patient.center=CombineData.XTeam(ICB.ID.xteams, file.query="PatientCenter/PatientCenter.rds")
profile.ICB.clinic=lapply(ICB.ID.xteams, function(x){
    if(all(c("Purity", "Ploidy")%in%colnames(patient.center[[x]]$SampleInfo))){
        tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=c("TreatmentResponse", "SampleTime", "TreatmentHistory", "CancerType", "OS", "OS.time", "Purity", "Ploidy", "TMB"))   %>%
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

combine.features=c("Checkpoint.overexprs", "ImmunoSuppressiveCell.overexprs", "ImmunoSuppressiveSig.overexprs", "ImmuneActivationGene.downregulation", "AntigenPresentGene.downregulation")

escape.num.list=lapply(ICB.ID.xteams, function(x){
    tmp.data=Filter(function(x) !all(is.na(x)), profile.trans.combine.escape[[x]])  %>%
        dplyr::filter(SampleID%in%common.sam.list[[x]])

    cluster.data=tmp.data[, intersect(colnames(tmp.data), c(combine.features, "SampleID"))]    %>%
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

OR.RECIST=list(response=c("Complete response", "Partial response", "Very good partial response", "Response"), nonresponse=c("Stable disease", "Progressive disease", "Non-response"))

plot.data=lapply(ICB.ID.xteams, function(x){
    data=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), list(escape.num.list[[x]], profile.ICB.clinic[[x]]))  %>%
        dplyr::mutate(TreatmentResponse=case_when(
                TreatmentResponse %in% OR.RECIST$response ~ "R",
                TreatmentResponse %in% OR.RECIST$nonresponse ~ "NR",
                TRUE ~ NA),
            TreatmentResponse=factor(TreatmentResponse, levels=c("R", "NR")))  %>%
        dplyr::filter(!is.na(SampleTime) & !is.na(TreatmentResponse))
})  %>%     setNames(ICB.ID.xteams)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/BioStat/TrendTestBetweenBinaryAndOrdinalVariables/CochranArmitageTrendTestWithPermutation.R")
CA.test.result=lapply(ICB.ID.xteams, function(x){
    df=plot.data[[x]]
    CA.data=na.omit(df[, c("TreatmentResponse", "Escape.burden", "es.feature.group")])
    if(nrow(CA.data)>100){
        CA.result=CochranArmitagePermu(CA.data, binary.vars="TreatmentResponse", ordinal.var="es.feature.group", is.pCorrect=FALSE, alternative="one.sided")
        data.frame(CA.result, Data=x)
    }else{
        NULL
    }
})      %>%     do.call(what=rbind)

all.CA.result=lapply(names(plot.data), function(x){
    df=plot.data[[x]]
    CA.data=na.omit(df[, c("TreatmentResponse", "es.feature.group")])
    if(nrow(CA.data)>100){
        tmp=table(CA.data)      %>%
            as.data.frame()        %>%
            dplyr::group_by(es.feature.group)       %>%
            dplyr::mutate(total=sum(Freq), ratio=ifelse(TreatmentResponse == "R", Freq / total, NA), Data=x) %>%
            na.omit()

        re=apply(tmp, 1, function(k){
            binom.test(as.numeric(k['Freq']), as.numeric(k['total']), conf.level = 0.95)$conf.int[1:2]   %>%
                t()     %>%     as.data.frame()           %>%
                setNames(c("low.ci", "high.ci"))
        })      %>%     do.call(what=rbind)
        tmp=cbind(tmp, re)
        return(tmp)
    }
    tmp=NULL
})  %>%     do.call(what=rbind)

show.Data=unique(all.CA.result$Data)
data.shape=setNames(rep_len(c(1, 7, 8, 9, 14, 21, 22, 23, 24), length(show.Data)), show.Data)
data.color=setNames(rep_len(c("#ce3375", "#60c8b3", "#6ea1d4", "#9f86c0"), length(show.Data)), show.Data)

tmp=CA.test.result
tmp$Data=factor(tmp$Data, unique(tmp$Data))

df=all.CA.result
df$Data=factor(df$Data, levels=unique(df$Data))
p1=ggplot(df, aes(x=es.feature.group, y=ratio, colour=Data, group=Data, shape=Data)) +
        geom_errorbar(aes(ymin=low.ci, ymax=high.ci), width=.1, position=position_dodge(0.1)) +
        geom_line(linewidth=0.5)+
        geom_point()+
        scale_color_manual(values=data.color)+
        scale_shape_manual(values=data.shape)+
        labs(x="es.feature.group", y="the ratio of ICB response", title=paste0("all non-genomic features")) +
        theme_classic()+
        theme(legend.position="none")+
        geom_text(data=tmp,
            aes(x=2, y=(0.1+(0.05)*as.integer(Data)), color=Data, label=paste0("p=", signif(p_value, 2), "\n")),
            hjust=1.1, vjust=-0.5, size=3.5)
p2=ggplot(df, aes(x=es.feature.group, y=Data)) +
        geom_tile(fill=NA) +
        geom_text(aes(label=paste0(Freq, "/", total), color=Data)) +
        scale_color_manual(values=data.color)+
        labs(x="es.feature.group", y="") +
        theme_test()+
        theme(legend.position="none")
p=plot_grid(p1, p2, ncol=1, align="v", axis="lr", rel_heights=c(2, 1))
ggsave(file.path(Dir.output, paste0("2.6.[CAtest][all_transcriptomicescape_feature]cumulative_countandICBresponse.pdf")), p, width=6, height=5)

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/KM.curve.and.logrank.test-2.0.R")
pdf(file.path(Dir.output, "2.7.all_transcriptomicescape_featureburden_groupsamplessurvival_difference.pdf"), 6, 7)
surv.data=lapply(names(escape.num.list), function(x){
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
                    color=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
        )
        print(p)
    }
})
dev.off()

tmp.dir=file.path(Dir.output, "nonGenomic")
if(!dir.exists(tmp.dir)){ dir.create(tmp.dir, recursive=TRUE) }
tmp.data=c("Liu_NatureMedicine_2019", "Hugo_Cell_2016", "Freeman_CellReportsMedicine_2022", "Ravi_NatureGenetics_2023")

pdf(file.path(tmp.dir, '3.escape_groupmultivariableCoxanalysis[based_onbinarizedescape_featurecumulative_count_combinedtranscriptomic][5yr].pdf'), 6, 4)
ICB.Multivariable.Cox.result=lapply(tmp.data, function(x){
    if(!is.null(profile.ICB.clinic[[x]])){
        df=dplyr::full_join(escape.num.list[[x]], profile.ICB.clinic[[x]], by="SampleID")   %>%
            dplyr::select(SampleTime, OS.time, OS, Escape.burden, es.feature.group, PatientID, Purity, Ploidy, TMB)      %>%
            na.omit()      %>%
            dplyr::mutate(
                cluster=factor(es.feature.group),
                log10TMB=log10(TMB+1),
                OS=ifelse(OS.time > 365*5, 0, OS),
                OS.time=ifelse(OS.time > 365*5, 365*5, OS.time)
            )

        if(length(table(df$es.feature.group))>1 & min(table(df$es.feature.group))!=0){
            os_form2=coxph(as.formula("Surv(OS.time, OS) ~ Escape.burden + log10TMB + Purity + Ploidy"), data=df)
            p2=ggforest(os_form2, data=df, main=paste0(x, "_Hazard ratio"), cpositions=c(0.02,-0.15, 0.25), fontsize=0.8, noDigits=2)
            print(p2)
        }
    }
})
dev.off()
