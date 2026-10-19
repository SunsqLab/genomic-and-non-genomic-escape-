# S4: Immune cell infiltration–TMB associations by cancer using Fisher's exact test
# This script evaluates immune infiltration, TMB, and immune escape relationships
# across TCGA cancer types. It prepares cohort-level profiles, applies the specified
# association models, and saves the resulting tables and figures.

library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(scatterpie); library(Cairo); library(openxlsx); library(writexl)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_immune_infiltration_affect_escape_feature_prevalence/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Danaher immune cell infiltration scores stored as a list.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
source('/pub5/xiaoyun/BioY/sunshangqin/5.Immunoediting/NewImmunoeditingMethodDesign/ImmuneEscape/ImmuneEscapeMechanisms/ImmuneCellInfiltration/ImmuneCellInfiltrationScores.R')
profile.exprs=CombineData.XTeam(ID.xteams, file.query='OMICSData/Exprs.data.rds')
profile.exprs$CRC_TCGA=do.call(cbind, CombineData.XTeam(c('COAD_TCGA', 'READ_TCGA'), file.query='OMICSData/Exprs.data.rds'))
profile.Danaher.cell=lapply(profile.exprs, function(tmp.exprs){
    ImmuneInfiltraScore(exprs.data=tmp.exprs)
})      %>%     setNames(ID.xteams)

# TMB profiles stored as data frames.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData/DataCenter/TCGA/PanCancer_TCGA/PatientCenter/PatientCenter.rds')
profile.TMB=GetInfor.PatientCenter(patient.center, colNames=c("SampleID", "TMB"))

# Compare the prevalence of hypermutated samples between high- and low-infiltration groups using Fisher's exact test.
cell.types=c('Cytotoxic.cell', 'CD8.Tcell', 'NK.cell', 'CD4.Tcell', 'Bcell', 'Neutrophils', 'Mast.cell')
fisher.result=lapply(ID.xteams, function(x){
    data=Reduce(function(x, y) dplyr::full_join(x, y, by='SampleID'), list(profile.Danaher.cell[[x]], profile.TMB))
    data=data[!is.na(data$TMB), ]

    df_grouped=data %>%
        dplyr::mutate(across(all_of(cell.types), ~ case_when(is.na(.) ~ NA_character_, . > median(., na.rm=TRUE) ~ "hi.infil", TRUE ~ "lo.infil")))     %>%
        dplyr::mutate(across(all_of(cell.types), ~ factor(., levels=c("hi.infil", "lo.infil"))),
            TMB=ifelse(TMB<10, "lowTMB", "highTMB")
        )

    tmp=lapply(cell.types, function(cell){
        df=na.omit(df_grouped[, c(cell, 'TMB')])
        df$TMB=factor(df$TMB, levels=c("highTMB", "lowTMB"))
        table_data=table(df$TMB, df[, cell])
        if(nrow(df)>0 && length(table_data)==4){
            re.fisher=fisher.test(table_data)
            OR=re.fisher$estimate
            if(min(table_data)==0){
                table_data=table_data+0.5
                OR=(table_data[1]*table_data[4])/(table_data[2]*table_data[3])
            }
            re=data.frame(CancerType=gsub("_TCGA", "", x), cell.type=cell, lowTMB.loinfil=table_data["lowTMB", "lo.infil"], lowTMB.hiinfil=table_data["lowTMB", "hi.infil"],
                highTMB.hiinfil=table_data["highTMB", "hi.infil"], highTMB.loinfil=table_data["highTMB", "lo.infil"], p_value=re.fisher$p.value, OddsRatio=OR)
        }else{
            re=NULL
        }
        return(re)
    })      %>%         do.call(what=rbind)
    tmp$FDR=p.adjust(tmp$p_value, method="BH")
    return(tmp)
})      %>%         do.call(what=rbind)
fisher.result$significance=ifelse(fisher.result$FDR<0.05 & fisher.result$OddsRatio>1, 'FDR<0.05 & OR>1', ifelse(fisher.result$FDR<0.05 & fisher.result$OddsRatio<1, 'FDR<0.05 & OR<1', 'no.sig'))
fisher.result$CancerType=factor(fisher.result$CancerType, levels=unique(fisher.result$CancerType))
fisher.result$cell.type=factor(fisher.result$cell.type, levels=unique(fisher.result$cell.type))

fisher.result[, c("lowTMB.loinfil", "lowTMB.hiinfil")]=fisher.result[, c("lowTMB.loinfil", "lowTMB.hiinfil")]/rowSums(fisher.result[, c("lowTMB.loinfil", "lowTMB.hiinfil")])
fisher.result[, c("highTMB.hiinfil", "highTMB.loinfil")]=fisher.result[, c("highTMB.hiinfil", "highTMB.loinfil")]/rowSums(fisher.result[, c("highTMB.hiinfil", "highTMB.loinfil")])
fisher.result$x=as.integer(fisher.result$CancerType)
fisher.result$y=as.integer(fisher.result$cell.type)
fisher.result$radius=ifelse(grepl('FDR<0.05', fisher.result$significance), 0.3, 0.1)
fisher.result$region=factor(1:nrow(fisher.result))

circle_border=do.call(rbind, lapply(seq_len(nrow(fisher.result)), function(i){
    theta=seq(0, 2*pi, length.out=100)
    data.frame(
        x=fisher.result$x[i] + fisher.result$radius[i] * cos(theta),
        y=fisher.result$y[i] + fisher.result$radius[i] * sin(theta),
        region=fisher.result$region[i],
        significance=fisher.result$significance[i]
    )
}))
circle_border=circle_border[circle_border$significance!="no.sig", ]

if(!exists("from_theme", mode="function")){
    from_theme=function(x){
        x_name=deparse(substitute(x))
        switch(x_name, ink="black", paper="white", accent="#3366FF", linewidth=0.5, borderwidth=0.5, NA)
    }
}

p=ggplot(fisher.result, aes(x=x, y=y, color=significance)) +
    scale_color_manual(name="significance", values=c("FDR<0.05 & OR>1"="red", "FDR<0.05 & OR<1"="#129990", "no.sig"="transparent")) +
    scatterpie::geom_scatterpie(
        data=fisher.result, aes(x=x, y=y, group=region, r=radius),
        cols=c("lowTMB.loinfil", "lowTMB.hiinfil", "highTMB.hiinfil", "highTMB.loinfil"),
        linetype="solid",
        color=NA) +
    geom_path(data=circle_border, aes(x=x, y=y, group=region, color=significance), linewidth=0.5)+
    coord_equal()+
    labs(x="", y="")+
    scale_fill_manual(name="isEscape", values=setNames(c("#FFE1E0", "#F49BAB", "#90D1CA", "#FFFBDE"), c("highTMB.loinfil", "highTMB.hiinfil", "lowTMB.loinfil", "lowTMB.hiinfil")))+
    scale_x_continuous(breaks=1:max(fisher.result$x), labels=levels(fisher.result$CancerType)) +
    scale_y_continuous(breaks=unique(fisher.result$y), labels=levels(fisher.result$cell.type)) +
    theme_classic()+
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom")+
    guides(fill=guide_legend(ncol=2), color=guide_legend(nrow=1, override.aes=list(size=1)))
ggsave(file.path(Dir.output, "2.1.associations_between_immune_cell_infiltration_and_median_dichotomized_TMB.bubble_pie_plot.pdf"), p, width=8, height=5)
