# S5: Positive selection of escape genes by TMB group
# This script evaluates selective pressure on immune escape alterations using TCGA
# mutation or copy-number data. It performs the indicated stratified analyses and
# saves gene- or feature-level statistical summaries and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(dndscv); library(sigminer); library(writexl); library(openxlsx)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_selection_pressure_affect_escape_features/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Clinical profiles including hypermutation status stored as a list.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R")
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.clinic=GetInfor.PatientCenter(patient.center, colNames=c("TMB", "frac.cnv"))

# Mutation profiles stored as a list.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R")
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c("FPPP_TCGA", "LAML_TCGA")))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.mut.df=CombineData.XTeam(ID.xteams, file.query="OMICSData/Mutations.data.rds")
profile.mut=lapply(ID.xteams, function(x){
    profile.mut.df[[x]][, c('SampleID', 'chromosome', 'startPosition', 'refAllele', 'mutAllele', 'geneSymbol')]
})      %>%     do.call(what=rbind)
# Retain primary tumor samples only.
profile.mut=profile.mut[GetInfor.PatientCenter(patient.center, SampleID=profile.mut$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]

# Data preparation and dNdScv calculation.
all.APM.gene=c("HLA-A", "HLA-B", "HLA-C",
        "B2M", "CALR", "TAP1", "TAP2", "TAPBP", "CIITA", "RFX5", "NLRC5",
        "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
        "CANX", "CD4", "CD74", "CD8A", "CD8B", "CREB1", "CTSB", "CTSL", "CTSS", "ERAP1", "ERAP2", "FAS", "HLA-E", "HLA-F", "HLA-G", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "HSPBP1", "IFI30", "IFNG", "IRF1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "KLRD1", "LGMN", "MEX3B", "NFYA", "NFYB", "NFYC", "PDIA3", "PSMA7", "PSMB10", "PSMB11", "PSMB6", "PSMB8", "PSMB9", "PSME1", "PSME2", "PSME3", "PSMF1", "RFXANK", "RFXAP", "TNF"
    )
IFNG.pathway.gene=c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")
other.gene=c('CD58', 'IDH1')
all.features.list=list(all.APM.gene, IFNG.pathway.gene, other.gene)
all.features=unlist(all.features.list)
names(all.features)=rep(c('#c82621', '#00b4d8', '#00b4d8'), lengths(all.features.list))

tmp.mut=profile.mut[, c('SampleID', 'chromosome', 'startPosition', 'refAllele', 'mutAllele', 'geneSymbol')]        %>%
    setNames(c('sampleID', 'chr', 'pos', 'ref', 'mut', 'geneSymbol'))       %>%
    dplyr::mutate(chr=gsub('chr', '', chr),
        TMB=profile.clinic$TMB[match(sampleID, profile.clinic$SampleID)])     %>%
    dplyr::filter(!is.na(TMB))     %>%
    dplyr::mutate(TMB=ifelse(TMB<10, "lowTMB", "highTMB"))

data("refcds_hg19", package="dndscv")
ref.genes=unlist(lapply(RefCDS, function(ref){ref$gene_name   }))
split.mut=split(tmp.mut, tmp.mut$TMB)
tmp.escape.dnds.group=lapply(split.mut, function(data){
        data=dndscv(mutations=data, refdb="hg19", outmats=T, gene_list=intersect(all.features, ref.genes))
    })      %>%     setNames(names(split.mut))

# Plot global selection pressure for all immune escape genes in hypermutated and non-hypermutated samples.
plots=lapply(names(split.mut), function(x){
        tmp=tmp.escape.dnds.group[[x]]$globaldnds   %>%
            dplyr::filter(name%in%c("wmis", "wtru", "wall"))
        p=ggplot(tmp, aes(x=name, y=mle)) +
            geom_errorbar(aes(ymin=cilow, ymax=cihigh), width=0.1, color="#99d98c") +
            geom_point(size=3, color="#34a0a4") +
            scale_x_discrete(labels=c("wmis"="Missense", "wtru"="Truncating", "wall"="All")) +
            geom_hline(yintercept=1, linetype="dashed", color="lightgrey") +
            theme_classic() +
            labs(x="", y="dN/dS ratio", title=paste0("PanCancer (", x, ")")) +
            theme(axis.text.x=element_text(angle=45, hjust=1))
})
p=plot_grid(plotlist=plots, nrow=1)
ggsave(file.path(Dir.output, "2.1.global_selection_pressure_of_escape_genes_by_TMB_group.pdf"), p, width=5, height=3)

# Plot gene-specific selection pressure for dNdScv-significant genes in hypermutated and non-hypermutated samples.
escape.dnds.group=lapply(tmp.escape.dnds.group, function(data){
        data$sel_cv
    })      %>%     setNames(names(tmp.escape.dnds.group))
sig.genes.group=unique(unlist(lapply(escape.dnds.group, function(data){
    data$gene_name[data$qallsubs_cv<0.1]
})))
max.sel=max(unlist(lapply(escape.dnds.group, function(data){
    max(data[, c('wmis_cv', 'wnon_cv')])
})))

pdf(file.path(Dir.output, "2.2.selection_pressure_by_mutation_type_and_TMB_group.pdf"), width=5, height=4)
result=lapply(names(escape.dnds.group), function(x){
    df1=escape.dnds.group[[x]][escape.dnds.group[[x]]$gene_name %in% sig.genes.group, ]     %>%
        dplyr::mutate(significance=ifelse(qallsubs_cv<0.1, "sig", "no.sig"))
    df2=reshape2::melt(df1, id.vars=c('gene_name', 'significance'), measure.vars=c('wmis_cv', 'wnon_cv'),
            variable.name='mutation.type', value.name='selection.pressure')       %>%
        dplyr::mutate(
            gene_name=factor(gene_name, levels=rev(sig.genes.group)),
            mutation.type=as.character(mutation.type),
            mutation.type=case_when(
                significance == 'no.sig' & mutation.type=='wmis_cv' ~ 'no.sig.mis',
                significance == 'no.sig' & mutation.type=='wnon_cv' ~ 'no.sig.trunc',
                mutation.type=='wmis_cv' ~ 'wmis',
                mutation.type=='wnon_cv' ~ 'wtrunc',
                TRUE ~ mutation.type
                )
        )

    p=ggplot(df2, aes(x=gene_name, y=selection.pressure, fill=mutation.type))+
        geom_bar(stat="identity", width=0.6, position=position_dodge(width=0.75))+
        ylim(c(0, ceiling(max.sel)))+
        labs(x='', title=x)+
        scale_fill_manual(values=c("wmis"='#f18f01', "wtrunc"="#048ba8", "no.sig.mis"="#BFBFC0", "no.sig.trunc"="#A3A4A5"), name="mutation.type")+
        geom_hline(yintercept=1, linetype="dashed", color="#3a405a")+
        theme_classic(base_size=12)+
        theme(axis.text.y=element_text(color=names(all.features)[match(rev(sig.genes.group), all.features)]))+
        coord_flip()
    print(p)
})
dev.off()
