# S5: Clonal and subclonal selection pressure by TMB group
# This script evaluates selective pressure on immune escape alterations using TCGA
# mutation or copy-number data. It performs the indicated stratified analyses and
# saves gene- or feature-level statistical summaries and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(dndscv); library(writexl); library(openxlsx)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_selection_pressure_affect_escape_features/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# TCGA pan-cancer mutation profile stored as a data frame.
profile.mut=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/OMICSData/Mutations.data.rds')
# Retain primary tumor samples only.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.mut=profile.mut[GetInfor.PatientCenter(patient.center, SampleID=profile.mut$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]

# Clinical profiles including hypermutation status stored as a list.
profile.clinic=GetInfor.PatientCenter(patient.center, colNames=c("TMB", "frac.cnv"))

# Genes corresponding to different immune escape features.
all.APM.gene=c(
        "HLA-A", "HLA-B", "HLA-C",
        "B2M", "CALR", "TAP1", "TAP2", "TAPBP", "CIITA", "RFX5", "NLRC5",
        "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
        "CANX", "CD4", "CD74", "CD8A", "CD8B", "CREB1", "CTSB", "CTSL", "CTSS", "ERAP1", "ERAP2", "FAS", "HLA-E", "HLA-F", "HLA-G", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "HSPBP1", "IFI30", "IFNG", "IRF1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "KLRD1", "LGMN", "MEX3B", "NFYA", "NFYB", "NFYC", "PDIA3", "PSMA7", "PSMB10", "PSMB11", "PSMB6", "PSMB8", "PSMB9", "PSME1", "PSME2", "PSME3", "PSMF1", "RFXANK", "RFXAP", "TNF"
    )
IFNG.pathway.gene=c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")
other.gene=c('CD58', 'IDH1')
all.features.list=list(all.APM.gene, IFNG.pathway.gene, other.gene)
all.features=unlist(all.features.list)
names(all.features)=rep(c('#c82621', '#00b4d8', '#00b4d8'), lengths(all.features.list))

# Calculate selection pressure on genomic immune escape genes using clonal and subclonal mutation data.
data("refcds_hg19", package="dndscv")
refcds.genes=unlist(lapply(RefCDS, function(ref)      ref$gene_name   ))

tmp.mut=profile.mut[, c("SampleID", "chromosome", "startPosition", "refAllele", "mutAllele", "geneSymbol", "isClonal")]   %>%
    setNames(c('sampleID', 'chr', 'pos', 'ref', 'mut', 'geneSymbol', 'isClonal'))       %>%
    dplyr::mutate(chr=gsub('chr', '', chr),
        TMB=profile.clinic$TMB[match(sampleID, profile.clinic$SampleID)])     %>%
    na.omit()       %>%
    dplyr::mutate(TMB=ifelse(TMB<10, "lowTMB", "highTMB"))
split.mut=split(tmp.mut, tmp.mut$TMB)

tmp.result=lapply(names(split.mut), function(TMB.status){
    data=split.mut[[TMB.status]]
    split.clonal=split(data, data$isClonal)
    re=lapply(names(split.clonal), function(is.clonal){
        tmp=NULL
        if(length(unique(split.clonal[[is.clonal]]$sampleID))>29){
            re=dndscv(mutations=split.clonal[[is.clonal]], refdb="hg19", gene_list=intersect(all.features, refcds.genes), outmats=T)
            tmp=geneci(dndsout=re, level=0.95)
            tmp[, c("qallsubs_cv", "qmis_cv", "qtrunc_cv")]=re$sel_cv[match(tmp$gene, re$sel_cv$gene_name), c("qallsubs_cv", "qmis_cv", "qtrunc_cv")]
            tmp$isClonal=ifelse(is.clonal=='TRUE', 'clonal', 'subclonal')
            tmp$TMB.status=TMB.status
        }
        return(tmp)
    })      %>%         do.call(what=rbind)
})      %>%     do.call(what=rbind)

sig.genes=tmp.result$gene[tmp.result$qmis_cv<0.1 | tmp.result$qtrunc_cv<0.1]     %>%     unique()
esgene.subclonal.selection=tmp.result[which(tmp.result$gene%in%sig.genes), ]
tmp=esgene.subclonal.selection[esgene.subclonal.selection$isClonal=='clonal', ]
esgene.subclonal.selection$gene=factor(esgene.subclonal.selection$gene, levels=rev(unique(tmp$gene[order(tmp$mis_mle)])))
esgene.subclonal.selection=esgene.subclonal.selection       %>%
    dplyr::mutate(
        mut.significance=isClonal,
        trunc.significance=isClonal,
        mut.significance=ifelse(!qmis_cv < 0.1, paste0("no.sig(", mut.significance, ")"), mut.significance),
        trunc.significance=ifelse(!qtrunc_cv < 0.1, paste0("no.sig(", trunc.significance, ")"), trunc.significance)
    )
mis.trunc.data=split(esgene.subclonal.selection, esgene.subclonal.selection$TMB.status)

# Examine selection pressure on each immune escape gene at clonal and subclonal levels for missense and truncating mutations.
plots.p=lapply(names(mis.trunc.data), function(TMB.status){
    p1=ggplot(data=mis.trunc.data[[TMB.status]], aes(x=gene, y=mis_mle))+
        geom_errorbar(aes(ymin=mis_low, ymax=mis_high, color=mut.significance), width=.5, position=position_dodge(0.7)) +
        geom_point(aes(color=mut.significance), position=position_dodge(0.7))+
        scale_color_manual(values=c("clonal"="#b5e2fa", "subclonal"="#0fa3b1", "no.sig(subclonal)"="lightgrey", "no.sig(clonal)"="lightgrey"))+
        labs(x='', y='selection.pressure', title=paste0(TMB.status, ' (missense mutation)'))+
        geom_hline(yintercept=1, linetype="dashed", color="red")+
        GGally::geom_stripped_cols()+
        theme_classic()+
        theme(axis.text.x=element_text(angle=45, hjust=1, color=names(all.features)[match(levels(mis.trunc.data[[TMB.status]]$gene), all.features)]),
            legend.position="right")
    p2=ggplot(data=mis.trunc.data[[TMB.status]], aes(x=gene, y=tru_mle))+
        geom_errorbar(aes(ymin=tru_low, ymax=tru_high, color=trunc.significance), width=.5, position=position_dodge(0.7)) +
        geom_point(aes(color=trunc.significance), position=position_dodge(0.7))+
        scale_color_manual(values=c("clonal"="#eddea4", "subclonal"="#f7a072", "no.sig(subclonal)"="lightgrey", "no.sig(clonal)"="lightgrey"))+
        labs(x='', y='selection.pressure', title=paste0(TMB.status, ' (truncating mutation)'))+
        geom_hline(yintercept=1, linetype="dashed", color="red")+
        GGally::geom_stripped_cols()+
        theme_classic()+
        theme(axis.text.x=element_text(angle=45, hjust=1, color=names(all.features)[match(levels(mis.trunc.data[[TMB.status]]$gene), all.features)]),
            legend.position="right")
    p=p1/p2
})

pdf(file.path(Dir.output, "3.1.clonal_and_subclonal_selection_pressure_of_escape_genes_by_TMB_group.pdf"), width=6, height=5)
print(plots.p)
dev.off()
