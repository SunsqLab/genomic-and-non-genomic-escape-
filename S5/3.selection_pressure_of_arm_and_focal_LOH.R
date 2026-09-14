# S5: Selection pressure of arm and focal LOH
# This script evaluates selective pressure on immune escape alterations using TCGA
# mutation or copy-number data. It performs the indicated stratified analyses and
# saves gene- or feature-level statistical summaries and figures.

library(magrittr); library(ggplot2); library(ggrepel); library(patchwork); library(dndscv); library(writexl); library(openxlsx)
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_selection_pressure_affect_escape_features/Question1"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Copy-number alteration profiles stored as a list.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R")
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c("FPPP_TCGA", "LAML_TCGA")))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.CNA=CombineData.XTeam(ID.xteams, file.query="OMICSData/CN.rds")
profile.CNA$PanCancer_TCGA=readRDS("/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/OMICSData/CN.rds")
# Retain primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.CNA=lapply(profile.CNA, function(x){
    tmp=x[GetInfor.PatientCenter(patient.center, SampleID=x$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
    tmp$group=ifelse(tmp$overlapArm<0.75, "focal", "arm")
    return(tmp)
})      %>%     setNames(names(profile.CNA))

# Test whether HLA.LOH and B2M.LOH prevalence exceeds the random background expectation.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/ImmuneEvasion/TestLOHEnrichmentInSpecificRegions/TestLOHEnrichmentInSpecificRegions.R")
# Test enrichment of focal-level LOH events.
sig.focalLOH.result=lapply(names(profile.CNA), function(x){
    print(x)
    tmp1=LOHPositiveSelection(CN.data=profile.CNA[[x]], genes=c("HLA-A", "HLA-B", "HLA-C"), LOH.range='focal')
    tmp2=LOHPositiveSelection(CN.data=profile.CNA[[x]], genes=c("B2M"), LOH.range='focal')
    tmp=rbind(tmp1, tmp2)
    data.frame(tmp, group=c("HLA.LOH", "B2M.LOH"), CancerType=gsub("_TCGA", "", x))
})      %>%     do.call(what=rbind)
sig.focalLOH.result$significance=ifelse(sig.focalLOH.result$p_value<0.05, "p<0.05", "nosig")
sig.focalLOH.result$CancerType=factor(sig.focalLOH.result$CancerType, levels=c("PanCancer", gsub("_TCGA", "", ID.xteams)))
sig.focalLOH.result$LOH.type="focal"
sig.focalLOH.result=sig.focalLOH.result[sig.focalLOH.result$group=="HLA.LOH", ]

# Test enrichment of arm-level LOH events.
sig.armLOH.result=lapply(names(profile.CNA), function(x){
    print(x)
    tmp1=LOHPositiveSelection(CN.data=profile.CNA[[x]], genes=c("HLA-A", "HLA-B", "HLA-C"), LOH.range='arm')
    tmp2=LOHPositiveSelection(CN.data=profile.CNA[[x]], genes=c("B2M"), LOH.range='arm')
    tmp=rbind(tmp1, tmp2)
    data.frame(tmp, group=c("HLA.LOH", "B2M.LOH"), CancerType=gsub("_TCGA", "", x))
})      %>%     do.call(what=rbind)
sig.armLOH.result$significance=ifelse(sig.armLOH.result$p_value<0.05, "p<0.05", "nosig")
sig.armLOH.result$CancerType=factor(sig.armLOH.result$CancerType, levels=c("PanCancer", gsub("_TCGA", "", ID.xteams)))
sig.armLOH.result$LOH.type="arm"
sig.armLOH.result=sig.armLOH.result[sig.armLOH.result$group=="HLA.LOH", ]

sig.LOH.result=rbind(sig.focalLOH.result, sig.armLOH.result)

# Plot selection pressure on LOH features across cancer types (bubble scatter plot).
p=ggplot(sig.LOH.result, aes(x=CancerType, y=LOH.type, color=significance, size=real.LOH.ratio))+
    geom_point(alpha=0.6, shape=7)+
    scale_color_manual(name="LOH significance", values=c("p<0.05"="#ff5e00", "nosig"="lightgrey"))+
    scale_size_continuous(name="the ratio of LOH samples")+
    labs(y="", title="")+
    coord_fixed(ratio=1) +
    theme_bw()+
    theme(axis.text.x=element_text(hjust=1, angle=45), axis.text.y=element_text(color="#c82621"), legend.position="bottom")
ggsave(file.path(Dir.output, paste0("1.3.enrichment_of_focal_and_arm_level_LOH_escape_features.pdf")), p, width=8, height=4)

# Scatter plot of LOH feature prevalence versus the genome-wide average LOH prevalence in each cancer type.
sig.LOH.result=list(sig.focalLOH.result, sig.armLOH.result)
names(sig.LOH.result)=c("focal", "arm")
plot.result=lapply(names(sig.LOH.result), function(y){
    sig.HLALOH.result=sig.LOH.result[[y]]    %>%
        dplyr::filter(group=="HLA.LOH" & CancerType!="PanCancer")
    p=ggplot(data=sig.HLALOH.result, aes(x=mean.LOH.genome.ratio, y=real.LOH.ratio, label=CancerType), shape=20)+
        geom_point(aes(color=significance))+
        scale_color_manual(values=c("p<0.05"="#fb8500", "no.sig"="#d6ccc2"))+
        geom_abline(intercept=0, slope=1, color="lightgrey", linetype="dashed")+
        geom_text_repel()+
        ggpubr::stat_cor(method="spearman", label.x=0) +
        labs(y="percentage of samples with HLA LOH", title=y)+
        theme_bw()

    ggsave(file.path(Dir.output, paste0("1.4.", y, "_LOH_rate_vs_observed_HLA_LOH_rate.pdf")), p, width=5, height=4)
})
