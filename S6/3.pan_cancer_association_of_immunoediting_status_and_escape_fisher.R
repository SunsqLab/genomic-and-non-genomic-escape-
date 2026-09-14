# S6: Pan-cancer association of immunoediting status and immune escape
# This script examines immunoediting, intratumor heterogeneity, or mutational
# signatures in relation to immune escape. It performs the indicated association
# analyses and saves the resulting statistical summaries and figures.

library(magrittr); library(ggplot2); library(patchwork); library(cowplot)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Association_between_immune_escape_and_immunoediting_scores/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Binarized immune escape profiles stored as a list.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds')
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
# Retain primary tumor samples only.
profile.escape.binarization=lapply(ID.xteams, function(x){
    data=profile.escape.binarization[[x]]
    tmp=lapply(data, function(y){
                y[GetInfor.PatientCenter(patient.center[[x]], SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
            })
    tmp=Reduce(function(x, y) dplyr::full_join(x, y, by="SampleID"), tmp)
    tmp$B2M.biallelic.inactivation=profile.B2M.inactive[[x]][1, match(tmp$SampleID, names(profile.B2M.inactive[[x]][1, ]))]
    tmp$B2M.biallelic.inactivation=ifelse(tmp$B2M.biallelic.inactivation==1, TRUE, FALSE)
    return(tmp)
})      %>%     setNames(ID.xteams)
profile.escape.binarization$PanCancer_TCGA=do.call(rbind, profile.escape.binarization)

# Immunoediting score profiles stored as a list.
IE.methods=c("immune.dNdS", "IE.HBMR")
profile.IE=lapply(ID.xteams, function(x){
	tmp=GetInfor.PatientCenter(patient.center[[x]], SampleID=profile.escape.binarization[[x]]$SampleID, colNames=IE.methods)
})      %>%     setNames(ID.xteams)
profile.IE$PanCancer_TCGA=do.call(rbind, profile.IE)

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

data=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(profile.IE$PanCancer_TCGA, profile.escape.binarization$PanCancer_TCGA))
plot.data.long=lapply(IE.methods, function(IE.method){
    tmp.data=Filter(function(x)  !all(is.na(x)), data)

    fisher.data=lapply(intersect(all.features, colnames(tmp.data)), function(feature){
        tmp=tmp.data[, c(feature, IE.method)]        %>%
            dplyr::filter(!!sym(IE.method) != 'Inf')        %>%
                na.omit()
        tmp[, IE.method]=ifelse(tmp[, IE.method]<1, 'edited', 'unedited')

        df=table(tmp[, c(feature, IE.method)])
        result=NULL
        if(length(df)==4){
            re=fisher.test(df)
            result=data.frame(escape.feature=feature, IE.method=IE.method, IE.status=rep(colnames(df), each=2),
                escape.status=rownames(df), Value=c(df), p.value=re$p.value, OddsRatio=re$estimate, row.names=NULL)
        }
        result$escape.feature=factor(result$escape.feature, levels=all.features)
        return(result)
    })      %>%     do.call(what=rbind)
    fisher.data$FDR=p.adjust(fisher.data$p.value, method="BH")
    return(fisher.data)
})      %>%     do.call(what=rbind)
plot.data.long$significance=ifelse(plot.data.long$FDR<0.01, "FDR<0.01", ifelse(plot.data.long$FDR<0.05, "FDR<0.05", "no.sig"))
plot.data.long$IE.method=factor(plot.data.long$IE.method, levels=rev(IE.methods))

# Use a heatmap to display associations between immune escape feature occurrence and immunoediting status.
p=ggplot(plot.data.long, aes(x=escape.feature, y=IE.method))+
    geom_point(aes(color=OddsRatio, size=significance), shape=15) +
    scale_color_gradient2(low="blue", high="red", mid="white", midpoint=1) +
    scale_size_manual(values=c("FDR<0.01"=4, "FDR<0.05"=2, "no.sig"=0.1))+
    labs(x='', y='')+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45, hjust=1, color=names(all.features)), legend.position="bottom")
ggsave(file.path(Dir.output, "2.1.association_between_immunoediting_status_and_escape_features.fisher.pdf"), p, width=10, height=3)

# Use bar plots to display associations between each immune escape feature and immunoediting status (Fisher's exact test).
pdf(file.path(Dir.output, "2.2.association_between_immunoediting_status_and_escape_features.fisher.bar_plot.pdf"), width=12, height=15)
plots=lapply(split(plot.data.long, plot.data.long$IE.method), function(data){
    p=ggplot(data, aes(x=IE.status, y=Value, fill=escape.status))+
        geom_bar(position="fill", stat="identity") +
        geom_text(aes(label=Value), position=position_fill(vjust=0.5), size=3)+
        scale_fill_manual(values=c("TRUE"="#ef7b64", "FALSE"="#8eb1de"))+
        geom_text(aes(x=1.5, y=1.1, label=paste('p =', signif(p.value, 2))),
                    colour="black", size=3, inherit.aes=FALSE, data=data) +
        scale_y_continuous(labels=scales::percent_format())+
        facet_wrap(~escape.feature, scales="free_x", ncol=7)+
        labs(title=data$IE.method[1], y='the proportions of samples')+
        theme_classic()+
        theme(legend.position='bottom')
    print(p)
})
dev.off()
