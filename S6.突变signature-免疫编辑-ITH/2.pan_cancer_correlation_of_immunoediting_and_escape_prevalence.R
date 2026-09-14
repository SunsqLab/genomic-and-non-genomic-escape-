# S6: Pan-cancer correlation of immunoediting and escape prevalence
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Association_between_immune_escape_and_immunoediting_scores/Question1'
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
        tmp$B2M.biallelic.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
        tmp$B2M.biallelic.inactivation=ifelse(tmp$B2M.biallelic.inactivation==1, TRUE, FALSE)
        return(tmp)
})      %>%     setNames(ID.xteams)

# Immunoediting score profiles stored as a list.
IE.methods=c("immune.dNdS", "IE.HBMR")
profile.IE=lapply(ID.xteams, function(x){
        patient.center=readRDS(paste0('/IData/DataCenter/TCGA/', x, '/PatientCenter/PatientCenter.rds'))
        GetInfor.PatientCenter(patient.center, SampleID=profile.escape.binarization[[x]]$SampleID, colNames=IE.methods)
})      %>%     setNames(ID.xteams)

# Within each cancer type, calculate the prevalence of each immune escape feature and the median immunoediting score.
all.features.list=list(
    APMalt=c("HLA.LOH", "B2M.biallelic.inactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),

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
all.features.color=setNames(c("#c82621", "#F6C141", "#fa8b69ff", "#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"), names(all.features.list))

# Calculate the proportion of immunoedited samples and immune escape feature prevalence within each cancer type.
plot.data=lapply(IE.methods, function(IE.method){
    tmp=lapply(ID.xteams, function(x){
        IE.data=profile.IE[[x]]     %>%
            dplyr::filter(!!sym(IE.method) != 'Inf' & !is.na(!!sym(IE.method)))

        if(nrow(IE.data)>29){
            IE.ratio=sum(IE.data[, IE.method] > 1, na.rm=TRUE)/nrow(IE.data)
            escape.ratio=colSums(profile.escape.binarization[[x]][, all.features], na.rm=TRUE)/colSums(!is.na(profile.escape.binarization[[x]][, all.features]))
            result=data.frame(CancerType=gsub('_TCGA', '', x), IE.method=IE.method, escape.ratio=escape.ratio, escape.feature=all.features, IE.ratio=IE.ratio)
        }else{
            result=NULL
        }
        return(result)
    })      %>%     do.call(what=rbind)
})      %>%     setNames(IE.methods)

pdf(file.path(Dir.output, "1.1.correlation_between_immunoediting_and_escape_feature_prevalence.across_cancers.pdf"), 13, 35)
result=lapply(plot.data, function(data){
    plots.p=lapply(all.features, function(m){
        data=na.omit(data[which(data$escape.feature==m), ])

        p=ggplot(data, aes(IE.ratio, escape.ratio)) +
            geom_point()+
            ggrepel::geom_text_repel(aes(label=CancerType), max.overlaps=20, segment.colour="lightgrey")+
            ggpubr::stat_cor(method="spearman", color='red', label.y=max(data$escape.ratio)+max(data$escape.ratio)/5)+
            labs(x=paste0('unedited.ratio (', data$IE.method[1], ')'))+
            theme_bw()+ labs(title=m)
    })
    print(plot_grid(plotlist=plots.p, ncol=4))
})
dev.off()

# For each feature, calculate the Spearman correlation between immunoedited sample prevalence and immune escape feature prevalence across cancers and visualize the correlations as a heatmap.
spear.data=lapply(plot.data, function(data){
    data=lapply(1:length(all.features), function(m){
        tmp=na.omit(data[which(data$escape.feature==all.features[m]), ])

        cor_test=cor.test(tmp$IE.ratio, tmp$escape.ratio, method="spearman", exact=FALSE)
        result=data.frame(feature=all.features[m], IE.method=data$IE.method[1], IE.ratio='non-edited.ratio', p.value=cor_test$p.value, cor=cor_test$estimate, group=names(all.features)[m])
    })  %>%     do.call(what=rbind)
    data$FDR=p.adjust(data$p.value, method="BH")
    return(data)
})      %>%     setNames(IE.methods)

plots=lapply(IE.methods, function(IE.method){
    data=spear.data[[IE.method]]    %>%
        dplyr::mutate(feature=ifelse(FDR<0.05, feature, ""))
    p=ggplot(data, aes(x=cor, y=log10(FDR)*(-1))) +
        geom_point(aes(color=group)) +
        geom_hline(yintercept=(-1)*log10(0.05), linetype="dashed", color="lightgrey")+
        scale_color_manual(values=all.features.color) +
        ggrepel::geom_text_repel(aes(label=feature), size=3) +
        theme_classic() +
        labs(y="-log10 (FDR)", x="Spearman's correlation", title=IE.method)
})
pdf(file.path(Dir.output, "1.2.correlation_between_escape_feature_and_immunoedited_sample_prevalence.pdf"), 5, 4)
print(plots)
dev.off()
