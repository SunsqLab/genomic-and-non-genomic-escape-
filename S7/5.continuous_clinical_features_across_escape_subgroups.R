# S7: Continuous clinical features across escape subgroups
# This script defines immune escape subgroups and evaluates their molecular or
# clinical characteristics across TCGA cohorts. It performs the indicated subgroup
# comparisons or survival models and saves the resulting figures.

library(magrittr); library(ggplot2); library(openxlsx); library(writexl)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Multiple_immune_escape_subgroup_patterns/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
clinic.metric=c("TMB", "TMB.clone", "TMB.subclone", "Neoantigen.burden.clone", "Neoantigen.burden.subclone", "prop.subclonal.mut",
    "Liuwei_2025", "IE.neoantigen.refSample", "IE.neoantigen.codon", "immune.dNdS", "IE.HBMR")
profile.clinical.info=lapply(ID.xteams, function(x){
    tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=clinic.metric)         %>%
        dplyr::mutate(PatientID=substr(SampleID, 1, 12))        %>%
        dplyr::mutate(across(c(TMB, TMB.subclone, TMB.clone, Neoantigen.burden.clone, Neoantigen.burden.subclone), ~ log10(.x + 1), .names="log10({.col}+1)"))  %>%
        dplyr::select(-c(TMB, TMB.subclone, TMB.clone, Neoantigen.burden.clone, Neoantigen.burden.subclone))

    tmp$prop.subclonal.mut[tmp$Liuwei_2025=='FAIL']=NA
    return(tmp)
})      %>%     setNames(ID.xteams)
profile.clinical.info$PanCancer_TCGA=do.call(rbind, profile.clinical.info)
new.clinic.metric=setdiff(colnames(profile.clinical.info$PanCancer_TCGA), c("PatientID", "SampleID"))

escape.score=c('gene19.Score', 'Vision.Score', 'ESrna.Score', 'ECMup.Score', 'escapeSig.Score', 'TIDE')
profile.escape.score=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/2.IntegratedImmuneEscapeAnalysis/ImmuneEscapeScore.rds')
profile.escape.score.df=do.call(rbind, profile.escape.score)

immune.info=read.xlsx(file.path(Dir.output, 'Immune landscape of cancer/NIHMS958212-supplement-2.xlsx'))
immune.metric=c('Intratumor.Heterogeneity', 'Proliferation', 'Aneuploidy.Score')
tmp.clinical.info=immune.info[, c('TCGA.Participant.Barcode', immune.metric)]      %>%
    setNames(c('PatientID', immune.metric))

profile.cibersort=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/20.ImmunePhenotype/ImmuneInfiltration/ImmuneInfiltrate.rds')
cell.types=c('T cells CD8', 'NK cells activated', 'Dendritic cells activated', 'Macrophages M1')
profile.cibersort=lapply(ID.xteams, function(x){
    data=profile.cibersort[[x]]$CIBERSORT.relative
    if(x=='CRC_TCGA'){
        data$SampleID=data$SampleID
    }else{
        data$SampleID=rownames(data)
    }
    data[data$'P-value'<0.05, c('SampleID', cell.types)]
})      %>%      do.call(what=rbind)

outDir=file.path(Dir.output, 'All_escape_features_binary_distance')
profile.cluster=readRDS(file.path(outDir, "Es.features.sum.combinedTrans.20260510.cluster.rds"))

ContinMetrics.list=lapply(tail(names(profile.clinical.info), 1), function(x){

    comparisons.re=lapply(names(profile.cluster), function(plus.features){
        combine.clinical.data=dplyr::full_join(profile.clinical.info[[x]], tmp.clinical.info, by='PatientID')
        data=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(profile.cluster[[plus.features]], combine.clinical.data, profile.cibersort, profile.escape.score.df))  %>%
            dplyr::filter(!is.na(SampleID))

        clini.var=setdiff(c(new.clinic.metric, immune.metric, cell.types, escape.score), "Liuwei_2025")
        tmp=lapply(clini.var, function(tmp.clini){
            plot.data=data[, c(tmp.clini, 'cluster', 'SampleID')]       %>%
                na.omit()       %>%         unique()        %>%
                dplyr::filter(!(tmp.clini %in% c("immune.dNdS", "IE.HBMR") & !!sym(tmp.clini) > 2))

            order.cluster=sort(unique(plot.data$cluster))
            x.label=paste0(order.cluster, ' (n=', table(plot.data$cluster)[order.cluster], ')')
            plot.data$cluster=factor(plot.data$cluster, levels=order.cluster)
            levels(plot.data$cluster)=x.label

            if (length(levels(plot.data$cluster)) >= 2) {
                pairwise_comparisons=combn(levels(plot.data$cluster), 2, simplify=FALSE)
                sig.pairs=purrr::map_dfr(pairwise_comparisons, function(pair) {
                        data_sub=plot.data[plot.data$cluster %in% c(pair[1], pair[2]), ]
                        test=wilcox.test(as.formula(paste0("`", tmp.clini, "` ~ cluster")), data=data_sub, exact=F)
                        data.frame(group1=pair[1], group2=pair[2], p.value=test$p.value)
                    })      %>%
                    dplyr::filter(p.value<0.05)
                group_list=Map(c, sig.pairs$group1, sig.pairs$group2)

                cluster.color=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
                p1=ggplot(data=plot.data, aes(x=cluster, y=!!sym(tmp.clini), color=cluster)) +
                    geom_violin(scale='width', alpha= 0.8) +
                    stat_summary(fun=median, geom="point", size=0.5) +
                    stat_summary(fun.data="median_hilow", fun.args=list(conf.int=0.5), geom="pointrange", size=0.5) +
                    stat_summary(fun=median, geom="text", aes(label=round(after_stat(y), 2)), color="black", position=position_nudge(x=0.2), vjust=-0.5, size=2)+
                    scale_color_manual(values=cluster.color)+
                    labs(title=x) +
                    theme_classic() +
                    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
                p2=p1+
                    ggpubr::stat_compare_means(comparisons=group_list, label="p.signif", method='wilcox.test', tip.length=0.01)

                p=list(p1, p2)
                return(p)
            }else{
                return(NULL)
            }
        })
    })      %>%     setNames(names(profile.cluster))
})      %>%     setNames(tail(names(profile.clinical.info), 1))

tmp=lapply(names(ContinMetrics.list), function(x){
    sapply(names(ContinMetrics.list[[x]]), function(plus.features){
        tmp.dir=file.path(outDir, x, plus.features)
        if(!dir.exists(tmp.dir)){ dir.create(tmp.dir, recursive=TRUE) }

        pdf(file.path(tmp.dir, '2.continuous_clinical_variables_across_escape_subgroups.pdf'), 3, 4)
        print(ContinMetrics.list[[x]][[plus.features]])
        dev.off()
    })
})
