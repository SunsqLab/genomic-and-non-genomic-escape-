# S5: Positive selection of LOH by CNA burden
# This script evaluates selective pressure on immune escape alterations using TCGA
# mutation or copy-number data. It performs the indicated stratified analyses and
# saves gene- or feature-level statistical summaries and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(dndscv); library(sigminer); library(writexl); library(openxlsx)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_selection_pressure_affect_escape_features/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Clinical profiles including hypermutation status stored as a list.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.clinic=GetInfor.PatientCenter(patient.center, colNames=c("frac.cnv"))

# Copy-number alteration profiles stored as a list.
profile.CNA=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/OMICSData/CN.rds')
# Retain primary tumor samples only.
profile.CNA=profile.CNA[GetInfor.PatientCenter(patient.center, SampleID=profile.CNA$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]

# Calculate selection pressure on HLA LOH and B2M LOH features (approximately 20 minutes each).
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/ImmuneEvasion/TestLOHEnrichmentInSpecificRegions/TestLOHEnrichmentInSpecificRegions.R')
data=unique(profile.clinic[, c('SampleID', 'frac.cnv')])   %>%
    dplyr::mutate(
        frac.cnv.group=ifelse(frac.cnv<median(frac.cnv, na.rm=T), 'low.CNA', 'high.CNA'),
        median.cnv=median(frac.cnv, na.rm=T))
sam.frac.cnv=data$frac.cnv.group[match(profile.CNA$SampleID, data$SampleID)]

split.cna=split(profile.CNA, sam.frac.cnv)

sig.HLALOH.focal.result=lapply(split.cna, function(data){
    tmp=LOHPositiveSelection(CN.data=data, genes=c('HLA-A', 'HLA-B', 'HLA-C'), LOH.range='focal')
})      %>%        do.call(what=rbind)
sig.HLALOH.focal.result$cnv.group=names(split.cna)

sig.HLALOH.arm.result=lapply(split.cna, function(data){
    tmp=LOHPositiveSelection(CN.data=data, genes=c('HLA-A', 'HLA-B', 'HLA-C'), LOH.range='arm')
})      %>%        do.call(what=rbind)
sig.HLALOH.arm.result$cnv.group=names(split.cna)

sig.LOH.result=rbind(sig.HLALOH.focal.result, sig.HLALOH.arm.result)      %>%
    dplyr::mutate(
        feature="HLA.LOH",
        cnv.region=rep(c('focal', 'arm'), c(nrow(sig.HLALOH.focal.result), nrow(sig.HLALOH.arm.result))),
        group=paste(cnv.region, cnv.group, sep='_')
    )

plots=lapply(split(sig.LOH.result, sig.LOH.result$cnv.group), function(data){
    data_long=tidyr::pivot_longer(data, cols=c(real.LOH.ratio, mean.LOH.genome.ratio), names_to="ratio.type", values_to="LOH.ratio")     %>%
        dplyr::mutate(color=ifelse(p_value<0.05, cnv.region, "no.sig"))

    p=ggplot(data_long, aes(x=cnv.region, y=LOH.ratio)) +
        geom_line(aes(group=cnv.region), color="grey", linewidth=1) +
        geom_point(aes(color=color, shape=ratio.type), size=3) +
        scale_color_manual(values=c("focal"="#e4572e", "arm"="#29335c", "no.sig"="#e5e5e5"))+
        scale_shape_manual(values=c("real.LOH.ratio"=19, "mean.LOH.genome.ratio"=17))+
        labs(title=data$cnv.group[1])+
        theme_classic()+
        theme(legend.position="right")
})
p=plot_grid(plotlist=plots, ncol=1)
ggsave(file.path(Dir.output, "2.3.enrichment_of_focal_and_arm_level_LOH_by_CNV_burden.pdf"), p, width=4, height=6)
