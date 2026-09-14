# S5: Selection pressure by mutation type in escape genes
# This script evaluates selective pressure on immune escape alterations using TCGA
# mutation or copy-number data. It performs the indicated stratified analyses and
# saves gene- or feature-level statistical summaries and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(dndscv)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_selection_pressure_affect_escape_features/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Mutation profiles stored as a list.
# Retain primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.mut=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/OMICSData/Mutations.data.rds')
profile.mut=profile.mut[GetInfor.PatientCenter(patient.center, SampleID=profile.mut$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]

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

# Extract selection pressure and significance for immune escape genes.
escape.dnds=readRDS(file.path(Dir.output, 'escape.dnds.rds'))
escape.select.df=readRDS(file.path(Dir.output, 'escape.select.df.rds'))
sig.genes=escape.select.df$gene_name[escape.select.df$qallsubs_cv<0.1 & escape.select.df$CancerType=='PanCancer']      %>%
    table()     %>%     names()

# Gene-specific selection pressure and confidence intervals.
tmp=geneci(dndsout=escape.dnds$PanCancer_TCGA, gene_list=all.features, level=0.95)
sig.escape.gene.ci=tmp[tmp$gene%in%sig.genes, ]

# Calculate mutation counts for individual immune escape genes in pan-cancer data.
tmp=escape.dnds$PanCancer_TCGA$genemuts
tmp2=tmp[tmp$gene_name%in%sig.genes, ]
escape.gene.mut=reshape2::melt(tmp2, id.vars="gene_name", measure.vars=c("n_syn", "n_mis", "n_non", "n_spl"), variable.name="mutType", value.name="mutations")

# Calculate the sample prevalence of mutations in each immune escape gene for every cancer type.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/ImmuneEvasion/GeneticImmuneEscape/ConstructGeneLevelMutationProfile.R")
Pan.mut=profile.mut[profile.mut$geneSymbol %in% sig.genes, ]
tmp.mis=ConstructMutationMatrix(mutation.data=Pan.mut, mut.types="missense")
tmp.trunc=ConstructMutationMatrix(mutation.data=Pan.mut, mut.types="truncating")

tmp2.mis=tmp.mis[intersect(rownames(tmp.mis), all.features), ]
tmp2.trunc=tmp.trunc[intersect(rownames(tmp.trunc), all.features), ]
sam.ratio.mis=rowSums(tmp2.mis)/length(unique(profile.mut$SampleID))
sam.ratio.trunc=rowSums(tmp2.trunc)/length(unique(profile.mut$SampleID))

t1=data.frame(gene_name=names(sam.ratio.mis), mutType='missense', sam.ratio=sam.ratio.mis)
t2=data.frame(gene_name=names(sam.ratio.trunc), mutType='truncating', sam.ratio=sam.ratio.trunc)
escape.mut.ratio=rbind(t1, t2)

pdf(file.path(Dir.output, "1.5.selection_pressure_by_mutation_type_in_escape_genes.pdf"), width=6, height=6)
data=sig.escape.gene.ci
if(nrow(data)!=0){
    data$gene=factor(data$gene, levels=data$gene[rev(order(data$tru_mle))])
    tmp1=data.frame(data[, c('gene', 'mis_mle', 'mis_low', 'mis_high')], type="missense")        %>%     dplyr::rename(selection.pressure=mis_mle, ci_low=mis_low, ci_high=mis_high)
    tmp2=data.frame(data[, c('gene', 'tru_mle', 'tru_low', 'tru_high')], type="truncating")        %>%     dplyr::rename(selection.pressure=tru_mle, ci_low=tru_low, ci_high=tru_high)
    df1=rbind(tmp1, tmp2)

    p1=ggplot(df1, aes(x=gene, y=selection.pressure, color=type))+
        geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.2, position=position_dodge(0.4))+
        geom_point(position=position_dodge(0.4), size=2)+
        scale_color_manual(values=c("missense"="#A20056FF", "truncating"="#008280FF"))+
        geom_text(aes(label=signif(selection.pressure, 2)), angle=90, size=2, nudge_x=0.1, nudge_y=2, vjust=-1)+
        labs(title='PanCancer', y="dN/dS ratio", x="")+
        geom_hline(aes(yintercept=1), linetype="dashed", color="grey")+
        theme_bw()+
        theme(axis.text.x=element_blank())

    df2=escape.gene.mut
    df2$gene_name=factor(df2$gene_name, levels=levels(data$gene))
    p2=ggplot(df2, aes(x=gene_name, y=mutations, fill=mutType))+
        geom_bar(position="stack",stat="identity")+
        ggsci::scale_fill_npg(labels=c('n_mis'='missense', 'n_non'='nonsense', 'n_spl'="splice", "n_syn"="synonymous"))+
        labs(x="", y="The number of mutations", fill="")+
        theme_classic()+
        theme(axis.text.x=element_blank())

    df3=escape.mut.ratio
    df3$gene_name=factor(df3$gene_name, levels=levels(data$gene))
    p3=ggplot(df3, aes(x=gene_name, y=sam.ratio, fill=mutType))+
        geom_bar(position='dodge', stat="identity")+
        scale_fill_manual(values=c('missense'='#cf648c', 'truncating'='#98c9a3'))+
        labs(x="", y="The ratio of samples", fill="")+
        theme_classic()+
        theme(axis.text.x=element_text(angle=45, hjust=1, color=names(all.features)[match(levels(data$gene), all.features)]))
    print(p1/p2/p3)
}
dev.off()
