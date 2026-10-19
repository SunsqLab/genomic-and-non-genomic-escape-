# S3: Association between TMB and HLA loss of heterozygosity in MSS tumors
# This script focuses on TCGA CRC, STAD, and UCEC tumors, retains MSS samples,
# and compares log10-transformed TMB between tumors with and without HLA LOH.
# The group difference is displayed as a box plot with a statistical comparison.
library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(writexl); library(RColorBrewer); require(ggpubr)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_mutation_rate_affect_escape_feature_prevalence/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Binarized immune escape profiles.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds')

# Primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData/DataCenter/TCGA/PanCancer_TCGA/PatientCenter/PatientCenter.rds')
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
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

# MSI subtype and TMB profiles for CRC, STAD, and UCEC.
ID.xteams.MSI=c('CRC_TCGA', 'STAD_TCGA', 'UCEC_TCGA')
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
profile.MSI=lapply(ID.xteams.MSI, function(x){
            GetInfor.PatientCenter(patient.center[[x]], colNames=c("MSISubtype2", "TMB"))
})      %>%     setNames(ID.xteams.MSI)

# TMB comparison by HLA LOH status in MSS tumors.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/BioStat/AssociationBetweenBinaryAndContinuousEscapeFeatures.Wilcoxon.R')
plot.data.MSI=lapply(ID.xteams.MSI, function(x){
    data=dplyr::full_join(profile.escape.binarization[[x]], profile.MSI[[x]], by='SampleID')    %>%
        dplyr::filter(MSISubtype2 %in% c('MSS'))   %>%
        dplyr::select(c(TMB, HLA.LOH))      %>%
        na.omit()   %>%
        dplyr::mutate(CancerType=x, log10.TMB=log10(TMB+1))
})          %>%     do.call(what=rbind)
p=ggboxplot(plot.data.MSI, x="CancerType", y="log10.TMB", color="HLA.LOH", palette="npg", add="jitter", add.params=list(alpha=0.5, size=0.4)) +
          stat_compare_means(aes(group=HLA.LOH), label="p.format")+
          ggtitle("MSS samples")
ggsave(file.path(Dir.output, "2.3.association_between_HLA_LOH_and_TMB_in_MSS_samples.pdf"), width=8, height=4)
