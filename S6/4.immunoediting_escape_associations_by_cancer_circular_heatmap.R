# S6: Immunoediting escape associations by cancer circular heatmap
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

library(magrittr); library(ggplot2); library(patchwork); library(cowplot); library(circlize); library(ComplexHeatmap); library(grid)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Association_between_immune_escape_and_immunoediting_scores/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Binarized immune escape profiles stored as a list.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData/DataCenter/TCGA/PanCancer_TCGA/PatientCenter/PatientCenter.rds')
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds')
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
# Retain primary tumor samples only.
profile.escape.binarization=lapply(ID.xteams, function(x){
    tmp=lapply(profile.escape.binarization[[x]], function(y){
                y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
            })
    tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
    tmp$B2M.biallelic.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
    tmp$B2M.biallelic.inactivation=ifelse(tmp$B2M.biallelic.inactivation==1, TRUE, FALSE)
    return(tmp)
})      %>%     setNames(ID.xteams)

# Immunoediting score profiles stored as a list.
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
IE.methods=c("immune.dNdS", "IE.HBMR")
profile.IE=lapply(ID.xteams, function(x){
	tmp=GetInfor.PatientCenter(patient.center[[x]], colNames=IE.methods)
})      %>%     setNames(ID.xteams)

all.features.list=list(
    APMalt=c("HLA.LOH", "B2M.biallelic.inactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),
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
names(all.features)=rep(c("#c82621", "#F6C141", "#fa8b69ff", "#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"), lengths(all.features.list))

source("/pub5/xiaoyun/BioY/sunshangqin/Functions/BioStat/AssociationBetweenTwoCategoricalVariables.FisherExactTest.R")
fisher.result=lapply(IE.methods, function(IE.method){
    fisher.data=lapply(ID.xteams, function(x){
        data=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(profile.IE[[x]], profile.escape.binarization[[x]]))
        data=Filter(function(x)  !all(is.na(x)), data)

        tmp.result=lapply(intersect(all.features, colnames(data)), function(feature){
            tmp.data=data[, c(feature, IE.method)]      %>%     na.omit()
            tmp=tmp.data[tmp.data[, IE.method]!='Inf', ]
            tmp[, IE.method]=ifelse(tmp[, IE.method]<1, 'edited', 'unedited')

            if(length(table(tmp[, c(feature, IE.method)]))==4){
                tmp=BinaryVarFisher(tmp, feature, IE.method)
                result=data.frame(CancerType=gsub("_TCGA", "", x), escape.feature=feature, IE.method=IE.method, tmp, row.names=NULL)
            }else{
                result=data.frame(CancerType=gsub("_TCGA", "", x), escape.feature=feature, IE.method=IE.method, odds_ratio=NA, p_value=NA,
                    ci_low=NA, ci_high=NA, FALSE.edited=NA, TRUE.edited=NA, FALSE.unedited=NA, TRUE.unedited=NA, row.names=NULL)
            }
            return(result)
        })      %>%     do.call(what=rbind)
        return(tmp.result)
    })          %>%     do.call(what=rbind)
    fisher.data$escape.feature=factor(fisher.data$escape.feature, levels=all.features)
    fisher.data$OR.padjust=log10(fisher.data$odds_ratio)
    fisher.data$CancerType=factor(fisher.data$CancerType, levels=unique(fisher.data$CancerType))
    return(fisher.data)
})          %>%    setNames(IE.methods)

# Generate a circular heatmap showing Fisher-test odds ratios and p-values for every immune escape feature and IE method, marking significant results with an asterisk.
source("/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/Script/ImmuneEscape/0.8.ImmunoeditingAndImmuneEscape/PlotCircularHeatmap.R")
pdf(file.path(Dir.output, "3.1.immunoediting_escape_associations_by_IE_method.circular_heatmap.fisher.pdf"), width=4, height=4)
plots=lapply(IE.methods, function(IE.method){
    data=fisher.result[[IE.method]]      %>%
        dplyr::filter(escape.feature%in%c("HLA.mut", "B2M.mut", "all.APM.mut", "IFNG.pathway.mut"))
    data$CancerType=factor(data$CancerType, levels=unique(data$CancerType))
    data$escape.feature=factor(data$escape.feature, levels=unique(data$escape.feature))

    df_OR=data        %>%
        dplyr::select(c(CancerType, escape.feature, OR.padjust))     %>%
        tidyr::pivot_wider(names_from=escape.feature, values_from=OR.padjust)        %>%
        column_to_rownames("CancerType")    %>%
        as.matrix()
    df_OR[is.na(df_OR)]=0

    data$x.index=as.integer(data$CancerType)
    data$y.index=rev(as.integer(data$escape.feature))
    tmp.data=data[which(data$p_value<0.05), ]

    CircularHeatmap(data=df_OR, label.x.color=NULL, pvalue.data=tmp.data, track.height=0.6,
        legend.title="log10(OddsRatio)", title=IE.method)

})
dev.off()

# Generate conventional heatmaps for the remaining cancer types.
escape.group=rep(names(all.features.list), lengths(all.features.list))       %>%     setNames(unlist(all.features.list))
plots=lapply(IE.methods, function(IE.method){
    data=fisher.result[[IE.method]]      %>%
        dplyr::filter(!escape.feature%in%c("HLA.mut", "B2M.mut", "all.APM.mut", "IFNG.pathway.mut"))   %>%
        na.omit()       %>%
        dplyr::mutate(significance=case_when(p_value<0.01 ~ "**", p_value<0.05 ~ "*", TRUE ~ ""),
            escape.feature=factor(escape.feature, levels=intersect(all.features, unique(escape.feature)))
        )      %>%
        tidyr::complete(CancerType, escape.feature)

    data$CancerType=factor(data$CancerType, levels=rev(unique(data$CancerType)))
    data$group=factor(escape.group[as.character(data$escape.feature)], levels=names(all.features.list))

    p1=ggplot(data, aes(escape.feature, y=1, fill=group))+
        geom_tile() +
        theme_minimal()+
        labs(x="", y="", title=IE.method)+
        coord_fixed(ratio=1)+
        theme(axis.text.y=element_blank(), axis.text.x=element_blank(), legend.position="top")+
        scale_fill_manual(values=setNames(unique(names(all.features)), names(all.features.list)))

    p2=ggplot(data, aes(x=escape.feature, y=CancerType, fill=OR.padjust)) +
        geom_tile(color="lightgrey") +
        scale_fill_gradient2(low="#3DE4BE", mid="white", high="#FF6027", midpoint=0,
            limits=c(ceiling(min(data$OR.padjust))-1, ceiling(max(data$OR.padjust))), na.value="white") +
        geom_text(aes(label=significance), color="black", size=4, vjust=0.5, hjust=0.5) +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom") +
        coord_fixed(ratio=1)+
        labs(x="", y="", fill="log10 (OR)")
    p=p1+p2+plot_layout(ncol=1)
})
pdf(file.path(Dir.output, "3.2.immunoediting_escape_associations_by_IE_method.heatmap.fisher.pdf"), width=10, height=10)
print(plots)
dev.off()

# Use bar plots to display associations between immunoediting status and immune escape features within individual cancer types.
pdf(file.path(Dir.output, "3.3.immunoediting_escape_associations_by_cancer.bar_plot.fisher.pdf"), width=9, height=12)
plots=lapply(IE.methods, function(IE.method){
    split.data=split(fisher.result[[IE.method]], fisher.result[[IE.method]]$CancerType)
    re=lapply(c("CRC", "SKCM", "STAD", "UCEC"), function(x){
        data=split.data[[x]]         %>%
            tidyr::pivot_longer(
                cols=starts_with(c("FALSE.edited", "TRUE.edited", "FALSE.unedited", "TRUE.unedited")),
                names_to="group", values_to="sam.num")     %>%
            tidyr::separate(col=group, into=c("escape.status", "IE.status"), sep="\\.") %>%
            na.omit()
        data$escape.feature=factor(data$escape.feature, levels=unique(data$escape.feature))
        p=ggplot(data, aes(x=IE.status, y=sam.num, fill=escape.status))+
            geom_bar(position="fill", stat="identity")+
            geom_text(aes(label=sam.num), position=position_fill(vjust=0.5), size=3)+
            scale_fill_manual(values=c("TRUE"="#ef7b64", "FALSE"="#8eb1de"))+
            geom_text(aes(x=1.5, y=1.1, label=paste('p =', signif(p_value, 2))),
                    colour="black", size=3, inherit.aes=FALSE, data=data) +
            scale_y_continuous(labels=scales::percent_format())+
            facet_wrap(~escape.feature, scales="free_x")+
            labs(title=paste0(x, ' (', IE.method, ')'), y='the proportions of samples')+
            theme_classic()+
            theme(legend.position='bottom')
        print(p)
    })
})
dev.off()
