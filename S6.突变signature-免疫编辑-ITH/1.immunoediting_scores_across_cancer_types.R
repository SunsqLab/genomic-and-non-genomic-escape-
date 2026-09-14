# S6: Immunoediting scores across cancer types
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Association_between_immune_escape_and_immunoediting_scores/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive = TRUE) }

# Immunoediting score profiles stored as a list.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam = "PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.IE.list=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/ImmuneEditing/ObservedToExpected/IE.rds')

source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
profile.IE.list=lapply(ID.xteams, function(x){
	patient.center=readRDS(paste0('/IData/DataCenter/TCGA/', x, '/PatientCenter/PatientCenter.rds'))
# Retain primary tumor samples only.
	tmp=GetInfor.PatientCenter(patient.center, SampleID=profile.IE.list[[x]]$SampleID, colNames=c("IE.neoantigen.refSample", "IE.neoantigen.codon", "immune.dNdS", "IE.HBMR", "SampleType"))
	tmp$CancerType=gsub('_TCGA', '', x)
	tmp[which(tmp$SampleType=="Primary"), ]
})      %>%     setNames(ID.xteams)
profile.IE.list$PanCancer_TCGA=do.call(rbind, profile.IE.list)
profile.IE.list$PanCancer_TCGA$CancerType="PanCancer"
profile.IE=do.call(rbind, profile.IE.list)

tmp.data=aggregate(profile.IE$IE.neoantigen.refSample, by=list(profile.IE$CancerType), FUN = function(x) median(x, na.rm = TRUE))
profile.IE$CancerType=factor(profile.IE$CancerType, levels=c(setdiff(tmp.data$Group.1[rev(order(tmp.data$x))], 'PanCancer'), 'PanCancer'))

plots=lapply(setdiff(colnames(profile.IE), c('SampleID', 'PatientID', 'CancerType', 'SampleType')), function(IE.method){
	p1=ggplot(profile.IE, aes(x=CancerType, y=!!sym(IE.method), color=CancerType))+
		geom_boxplot()+
		geom_hline(yintercept=1, linetype="dashed", color="red")+
		theme_bw()+
		theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
})
p=plots[[1]]/plots[[2]]/plots[[3]]/plots[[4]]
ggsave(file.path(Dir.output, '0.immunoediting_scores.pdf'), p, width=15, height=10)
