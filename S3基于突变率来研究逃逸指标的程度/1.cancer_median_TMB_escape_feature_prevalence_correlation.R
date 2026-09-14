# S3: Association between tumor mutational burden and immune escape feature prevalence
# For each TCGA cancer type, this script calculates the median tumor mutational
# burden (TMB) and the prevalence of individual immune escape features. Spearman
# correlation is then used to assess their associations across cancer types.
# Results are presented as feature-level scatter plots and an FDR summary plot.
library(magrittr); library(cowplot); library(ggplot2); library(patchwork)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_mutation_rate_affect_escape_feature_prevalence/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Binarized immune escape profiles.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds')
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")

# Primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
profile.escape.binarization=lapply(ID.xteams, function(x){
			tmp=lapply(profile.escape.binarization[[x]], function(y){			
						y[GetInfor.PatientCenter(patient.center[[x]], SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
					})
            tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
            tmp$biallelic.B2M.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
            tmp$biallelic.B2M.inactivation=ifelse(tmp$biallelic.B2M.inactivation==1, TRUE, FALSE)
            return(tmp)
		})      %>%     setNames(ID.xteams)

# TMB profiles.
profile.TMB=lapply(ID.xteams, function(x){
    GetInfor.PatientCenter(patient.center[[x]], colNames=c("SampleID", "TMB")) 
})      %>%     setNames(ID.xteams)

# Immune escape feature groups.
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
names(all.features)=rep(names(all.features.list), lengths(all.features.list))

# Cancer-level median TMB and immune escape feature prevalence.
plot.data=lapply(ID.xteams, function(x){
    tmp1.data=dplyr::full_join(profile.escape.binarization[[x]], profile.TMB[[x]], by='SampleID')
    tmp2.data=tmp1.data[, c(all.features, 'TMB')]
    tmp2.data=Filter(function(x) !all(is.na(x)), tmp2.data)
    tmp2.data=tmp2.data[!apply(tmp2.data, 1, function(row) all(is.na(row))), ]

    median_TMB=median(tmp2.data[, 'TMB'], na.rm=TRUE)
    tmp.features=intersect(all.features, colnames(tmp2.data))
    escape.ratio=colSums(tmp2.data[, tmp.features], na.rm=TRUE)/nrow(tmp2.data)

    result=data.frame(CancerType=gsub('_TCGA', '', x), escape.ratio, escape.feature=tmp.features, median_TMB)
})      %>%     do.call(what=rbind)

# Spearman correlations and scatter plots across cancer types.
pdf(file.path(Dir.output, "1.cancer_TMB_vs_escape_feature_prevalence.pan_cancer.pdf"), 14, 26)
plots.p=lapply(all.features, function(m){
    data=na.omit(plot.data[which(plot.data$escape.feature==m), ])

    p=ggplot(data, aes(median_TMB, escape.ratio)) + 
        geom_point(size=0.8)+
        geom_smooth(method='lm', aes(x=median_TMB, y=escape.ratio), formula=y ~ x, color="#48cae4")+
        ggrepel::geom_text_repel(aes(label=CancerType), max.overlaps=20, segment.colour="lightgrey", size=2.5)+
        ggpubr::stat_cor(method="spearman", cor.coef.name=expression(rho), size=3, color='red', label.y=max(data$escape.ratio)+max(data$escape.ratio)/5)+
        labs(x='the median TMB of each cancer', y=paste0('escape.ratio'), title=m) +
        theme_test() +
        theme(panel.border=element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"))
})
print(plot_grid(plotlist=plots.p, ncol=5))
dev.off()

# FDR-adjusted correlation summary.
data=lapply(1:length(all.features), function(m){
    tmp=na.omit(plot.data[plot.data$escape.feature==all.features[m], ])

    cor_test=cor.test(tmp$median_TMB, tmp$escape.ratio, method="spearman", exact=FALSE)
    result=data.frame(feature=all.features[m], median_TMB='median_TMB', p_value=cor_test$p.value, cor=cor_test$estimate, group=names(all.features)[m])
})  %>%     do.call(what=rbind)
data$group=factor(data$group, levels=names(all.features.list))
data$feature=factor(data$feature, levels=all.features)
data$FDR=p.adjust(data$p_value, method="BH")

p=ggplot(data, aes(x=feature, y=(-1)*log10(FDR), fill=cor))+
    geom_bar(stat="identity", colour="black", width=0.78, position=position_dodge(0.7))+
    labs(x="", y="-log10(FDR)")+
    scale_fill_gradient2(low="#053c5e", high="#a31621", mid="white", midpoint=0, limits=c(-1, 1), name="spearman cor")+
    facet_grid(~ group, space="free_x", scales="free_x")+
    geom_hline(yintercept=log10(0.05)*(-1), linetype="dashed", linewidth=0.4)+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom")
        
ggsave(file.path(Dir.output, "1.cancer_TMB_vs_escape_feature_prevalence.pan_cancer.barplot.pdf"), p, width=15, height=5)

