# S4: Immune infiltration escape feature associations by cancer
# This script evaluates immune infiltration, TMB, and immune escape relationships
# across TCGA cancer types. It prepares cohort-level profiles, applies the specified
# association models, and saves the resulting tables and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(writexl); library(openxlsx)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_infiltration_affect_escape_feature_prevalence/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Binarized immune escape profiles stored as a list.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R")
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds')
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
# Retain primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData3/DataCenter/IntegratedData/PanCancer_TCGA.dataset/PatientCenter/PatientCenter.rds')
profile.escape.binarization=lapply(ID.xteams, function(x){
            data=profile.escape.binarization[[x]]
			tmp=lapply(data, function(y){
						y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
					})
            tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
            tmp$biallelic.B2M.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
            tmp$biallelic.B2M.inactivation=ifelse(tmp$biallelic.B2M.inactivation==1, TRUE, FALSE)
            return(tmp)
		})      %>%     setNames(ID.xteams)

# Danaher immune cell infiltration scores stored as a list.
profile.exprs=CombineData.XTeam(ID.xteams, file.query='OMICSData/Exprs.data.rds')
profile.exprs$CRC_TCGA=do.call(cbind, CombineData.XTeam(c('COAD_TCGA', 'READ_TCGA'), file.query='OMICSData/Exprs.data.rds'))
source('/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/ImmuneEscape/ImmuneEscapeMechanisms/ImmuneCellInfiltration/ImmuneCellInfiltrationScores.R')
profile.Danaher.cell=lapply(profile.exprs, function(x){
    ImmuneInfiltraScore(exprs.data=x)
})      %>%     setNames(ID.xteams)

all.features.list1=list(
    APMalt=c("HLA.LOH", "biallelic.B2M.inactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),
    Checkalt=c("CD274.deepAmp"),
    ActMalt=c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut")
)
all.features1=unlist(all.features.list1)
names(all.features1)=rep(c("#c82621", "#F6C141", "#fa8b69ff"), lengths(all.features.list1))

all.features.list2=list(
    Checkexp=c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    Supprcell=c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    Supprsig=c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    Actdown=c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown=c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo=c("mean.neo.exprs")
)
all.features2=unlist(all.features.list2)
names(all.features2)=rep(c("#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"), lengths(all.features.list2))

cell.types=c('Cytotoxic.cell', 'CD8.Tcell', 'NK.cell', 'CD4.Tcell', 'Bcell', 'Neutrophils', 'Mast.cell')

# Within each cancer type, test whether immune infiltration differs between samples with and without each immune escape feature.
plot.data=lapply(ID.xteams, function(x){
    data=dplyr::full_join(profile.escape.binarization[[x]], profile.Danaher.cell[[x]], by='SampleID')

    tmp2=lapply(cell.types, function(cell){
        tmp1=lapply(intersect(names(data), c(all.features1, all.features2)), function(m){
            df=na.omit(data[, c(m, cell)])

            result=NULL
            if(length(unique(df[, m]))==2 && !min(table(df[, m]))<5){
                wilcox.result=wilcox.test(df[, cell] ~ df[, m], df, exact=FALSE)
                median.value=tapply(df[, cell], df[, m], median, na.rm=TRUE)

                p_value=wilcox.result$p.value
                fold_change=(median.value['TRUE']+10^(-6))/(median.value['FALSE']+10^(-6))
                log2FC=log2(fold_change)

                result=data.frame(CancerType=gsub('_TCGA', '', x), escape.feature=m, cell.type=cell, p_value=p_value, fold_change=fold_change, log2FC=log2FC)
            }
            return(result)
        })      %>%     do.call(what=rbind)
    })      %>%     do.call(what=rbind)
    tmp2$FDR=p.adjust(tmp2$p_value, method="BH")
    return(tmp2)
})          %>%     do.call(what=rbind)
plot.data$CancerType=factor(plot.data$CancerType, levels=unique(plot.data$CancerType))
plot.data$significance=ifelse(plot.data$FDR<0.01, "FDR<0.01", ifelse(plot.data$FDR<0.05, "FDR<0.05", "no.sig"))

# Plot Wilcoxon association results between selected cell infiltration levels and immune escape features within each cancer type (bubble plot).
my_colors=c("#023e8a", "#3bb273", "#fdfcdc", "#ee4266", "#540d6e")

plots=lapply(list(all.features1, all.features2), function(tmp.features){
    if("HLA.LOH"%in%tmp.features){
        pdf(file.path(Dir.output, "1.1.genomic_escape_T_cell_infiltration_associations_by_cancer.pdf"), width=8, height=3)
    }else{
        pdf(file.path(Dir.output, "1.1.non_genomic_escape_T_cell_infiltration_associations_by_cancer.pdf"), width=8, height=9)
    }

    data=plot.data[which(plot.data$escape.feature%in%tmp.features), ]
    plots.p=lapply(cell.types, function(cell){
        tmp.data=data[which(data$cell.type==cell), ]
        tmp.data$escape.feature=factor(tmp.data$escape.feature, levels=rev(intersect(tmp.features, tmp.data$escape.feature)))

        min.fd=(min(tmp.data$log2FC))
        max.fd=(max(tmp.data$log2FC))
        my_values=c(min.fd, min.fd/3, 0, (max.fd/3), max.fd)
        p=ggplot(tmp.data, aes(x=CancerType, y=escape.feature, fill=log2FC, color=log2FC, size=significance)) +
            geom_point(shape=21) +
            scale_fill_gradientn(colors=my_colors, values=scales::rescale(my_values)) +
            scale_color_gradientn(colors=my_colors, values=scales::rescale(my_values)) +
            scale_size_manual(values=c("FDR<0.01"=4, "FDR<0.05"=2.5, "no.sig"=0.1)) +
            guides(fill=guide_colorbar(direction="horizontal"))+
            labs(x='', y='', title=cell, fill="log2FC", color="log2FC") +
            theme_test() +
            theme(
                legend.position="bottom",
                axis.text.x=element_text(angle=45, hjust=1, size=10),
                axis.text.y=element_text(size=10, color=names(tmp.features)[match(levels(tmp.data$escape.feature), tmp.features)])
            )
        print(p)
    })
    dev.off()
})

tmp.data=plot.data      %>%
    dplyr::filter(escape.feature%in%c("HLA.LOH", "biallelic.B2M.inactivation", "all.APM.mut", "IFNG.pathway.mut"), cell.type=="Cytotoxic.cell", significance!="no.sig")        %>%
    dplyr::mutate(fold_change=ifelse(fold_change<1, "negative", ifelse(fold_change>1, "positive", NA)),
        fold_change=factor(fold_change, levels=c("positive", "negative")),
        label=paste0(CancerType, "\n", escape.feature))

library(ggsci)
shape.values=c(22, 23, 24, 25)         %>%         setNames(c("HLA.LOH", "biallelic.B2M.inactivation", "all.APM.mut", "IFNG.pathway.mut"))
color.values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#0073C2FF", "#EFC000FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#868686FF", "#91D1C2FF")    %>%
    setNames(c("CRC", "STAD", "UCEC", "HNSC", "LUAD", "CESC", "KIRP", "LUSC", "BRCA", "SARC"))

p=ggplot(tmp.data, aes(x=log2FC, y=(-1)*log10(FDR), shape=escape.feature, color=CancerType, fill=CancerType))+
    geom_hline(yintercept=(-1)*log10(0.05), color="lightgrey", linetype="dashed", linewidth=0.5)+
    geom_vline(xintercept=0, color="lightgrey", linetype="dashed", linewidth=0.5)+
    geom_point()+
    scale_color_manual(values=color.values)+
    scale_fill_manual(values=color.values)+
    scale_shape_manual(values=shape.values)+
    ggrepel::geom_text_repel(aes(label=label), size=1.5, color="#868686FF", max.overlaps=Inf)+
    labs(y="-log10 (FDR)", x="log2FC")+
    theme_classic(base_size=10)
ggsave(file.path(Dir.output, "1.2.T_cell_infiltration_associations_with_selected_genomic_escape_features_by_cancer.pdf"), p, width=5, height=4)

data=plot.data      %>%
    dplyr::filter(CancerType=='CRC', escape.feature%in%unlist(all.features.list1[c("APMalt", "Checkalt", "ActMalt")]))  %>%
    dplyr::group_by(cell.type)      %>%
    dplyr::mutate(
        num.sig=sum(significance!="no.sig" & fold_change>1)/length(unlist(all.features.list1[c("APMalt", "Checkalt", "ActMalt")])),
        cell.type=factor(cell.type, levels=cell.types))     %>%
    dplyr::filter(significance!='no.sig' & fold_change>1)   %>%
    dplyr::mutate(
        escape.feature=paste(escape.feature, "\n", sep='', collapse=''))      %>%
    dplyr::select(c(cell.type, num.sig, escape.feature))      %>%
    unique()    %>%
    as.data.frame()

color.values=setNames(rep(c("#ffbe0b", "#d5bdaf"), c(3, 4)), cell.types)
p=ggplot(data, aes(x=cell.type, y=num.sig, fill=cell.type))+
    geom_col()+
    scale_fill_manual(values=color.values)+
    ylim(c(0, 1))+
    labs(y="percent of significant escape features", title="CRC")+
    geom_text(aes(label=escape.feature), size=1.8, color="#6d7278", vjust=0)+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
ggsave(file.path(Dir.output, "1.3.associations_between_genomic_escape_features_and_cell_types_in_CRC.pdf"), p, width=5, height=5)
