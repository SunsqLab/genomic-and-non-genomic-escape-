# S5: Mutational selection pressure in escape genes
# This script evaluates selective pressure on immune escape alterations using TCGA
# mutation or copy-number data. It performs the indicated stratified analyses and
# saves gene- or feature-level statistical summaries and figures.

library(magrittr); library(ggplot2); library(ggrepel); library(patchwork); library(dndscv); library(writexl); library(openxlsx)
Dir.output="/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_selection_pressure_affect_escape_features/Question1"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Mutation profiles stored as a list.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R")
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c("FPPP_TCGA", "LAML_TCGA")))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.mut=CombineData.XTeam(ID.xteams, file.query="OMICSData/Mutations.data.rds")
profile.mut$PanCancer_TCGA=readRDS("/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/OMICSData/Mutations.data.rds")
# Retain primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.mut=lapply(profile.mut, function(x){
    x[GetInfor.PatientCenter(patient.center, SampleID=x$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
})      %>%     setNames(names(profile.mut))

all.APM.gene=c("HLA-A", "HLA-B", "HLA-C",
        "B2M", "CALR", "TAP1", "TAP2", "TAPBP", "CIITA", "RFX5", "NLRC5",
        "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
        "CANX", "CD4", "CD74", "CD8A", "CD8B", "CREB1", "CTSB", "CTSL", "CTSS", "ERAP1", "ERAP2", "FAS", "HLA-E", "HLA-F", "HLA-G", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "HSPBP1", "IFI30", "IFNG", "IRF1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "KLRD1", "LGMN", "MEX3B", "NFYA", "NFYB", "NFYC", "PDIA3", "PSMA7", "PSMB10", "PSMB11", "PSMB6", "PSMB8", "PSMB9", "PSME1", "PSME2", "PSME3", "PSMF1", "RFXANK", "RFXAP", "TNF"
    )
IFNG.pathway.gene=c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")
other.gene=c("CD58", "IDH1")
all.features.list=list(all.APM.gene, IFNG.pathway.gene, other.gene)
all.features=unlist(all.features.list)
names(all.features)=rep(c("#c82621", "#00b4d8", "#00b4d8"), lengths(all.features.list))

# Use dNdScv to calculate selection pressure on mutation-based immune escape genes.
# Calculate the sample prevalence of mutations in each immune escape gene.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/ImmuneEvasion/GeneticImmuneEscape/ConstructGeneLevelMutationProfile.R")
escape.mut.ratio=lapply(names(profile.mut), function(x){
    data=profile.mut[[x]][profile.mut[[x]]$geneSymbol %in% all.features, ]
    tmp1=ConstructMutationMatrix(mutation.data=data, mut.types="nonsilent")
    tmp2=tmp1[intersect(rownames(tmp1), all.features), ]
    rowSums(tmp2)/length(unique(profile.mut[[x]]$SampleID))
})      %>%     setNames(names(profile.mut))

# Calculate gene-specific and global selection pressure for immune escape genes within individual cancer types and pan-cancer data using dNdScv (approximately 10 minutes).
data("refcds_hg19", package="dndscv")
ref.genes=unlist(lapply(RefCDS, function(ref){ref$gene_name   }))
escape.dnds=lapply(names(profile.mut), function(x){
    tmp.mut=profile.mut[[x]][, c("SampleID", "chromosome", "startPosition", "refAllele", "mutAllele", "geneSymbol")]        %>%
        setNames(c("sampleID", "chr", "pos", "ref", "mut", "geneSymbol"))       %>%
        dplyr::mutate(chr=gsub("chr", "", chr))
    gene_list=Reduce(intersect, list(all.features, ref.genes, tmp.mut$geneSymbol))
    if(sum(tmp.mut$geneSymbol%in%gene_list)>10){
        tmp.dnds=dndscv(mutations=tmp.mut, refdb="hg19", outmats=T, gene_list=gene_list)
    }else{
        tmp.dnds=NULL
    }
    return(tmp.dnds)
})    %>%   setNames(names(profile.mut))
escape.dnds.list=escape.dnds[lengths(escape.dnds)!=0]

escape.select.df=lapply(names(escape.dnds.list), function(x){
    tmp=escape.dnds.list[[x]]$sel_cv
    tmp$CancerType=gsub("_TCGA", "", x)
    tmp$mut.ratio=escape.mut.ratio[[x]][match(tmp$gene_name, names(escape.mut.ratio[[x]]))]
    return(tmp)
})      %>%     do.call(what=rbind)

tmp.escape.select.df=escape.select.df[escape.select.df$qallsubs_cv<0.1, ]
write_xlsx(list("Table S1"=tmp.escape.select.df), path=file.path(Dir.output, "Selection_pressure_of_genomic_escape_features.xlsx"))

# Plot dNdScv results, including global and gene-specific selection pressure for immune escape genes.
# Plot global selection pressure for all immune escape genes.
data=escape.dnds.list$PanCancer_TCGA$globaldnds
p=ggplot(data, aes(x=name, y=mle)) +
    geom_errorbar(aes(ymin=cilow, ymax=cihigh), width=0.1, color="#99d98c") +
    geom_point(size=3, color="#34a0a4") +
    scale_x_discrete(labels=c("wmis"="Missense", "wnon"="Nonsense", "wspl"="Splice", "wtru"="Truncating", "wall"="All")) +
    geom_hline(yintercept=1, linetype="dashed", color="lightgrey") +
    theme_classic() +
    labs(x="", y="dN/dS ratio", title="PanCancer") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(Dir.output, "1.1.global_selection_pressure_of_escape_genes.pdf"), p, width=2, height=3)

# Generate a bubble plot of gene-specific selection pressure.
sig.genes=escape.select.df$gene_name[escape.select.df$qallsubs_cv<0.1]      %>%
    table()     %>%     sort()      %>%     names()

mut.esgene.data=escape.select.df %>%
    dplyr::filter(gene_name %in% sig.genes) %>%
    dplyr::mutate(group=ifelse(wmis_cv > 1 | wnon_cv > 1, "positive", "negative"),
        group=ifelse(qallsubs_cv > 0.1, "no.sig", group),
        gene_name=factor(gene_name, levels=sig.genes),
        CancerType=factor(CancerType, levels=c("PanCancer", gsub("_TCGA", "", ID.xteams))))

p=ggplot(mut.esgene.data, aes(x=CancerType, y=gene_name, color=group, size=mut.ratio))+
    geom_point(alpha=0.6)+
    scale_x_discrete(drop = FALSE) +
    scale_color_manual(name="mut significance", values=c("positive"="#023e8a", "negative"="#52b788", "no.sig"="lightgrey"))+
    scale_size_continuous(name="the ratio of mut samples")+
    geom_hline(yintercept=0, color="red", linetype="dotted", linewidth=1)+
    guides(color=guide_legend(nrow=2))+
    labs(x="", y="")+
    theme_bw()+
    theme(
        axis.text.x=element_text(angle=45, hjust=1),
        axis.text.y=element_text(color=names(all.features)[match(levels(mut.esgene.data$gene_name), all.features)]),
        legend.position="top"
)
ggsave(file.path(Dir.output, "1.2.gene_specific_selection_pressure_of_escape_genes.pdf"), p, width=7, height=7)
