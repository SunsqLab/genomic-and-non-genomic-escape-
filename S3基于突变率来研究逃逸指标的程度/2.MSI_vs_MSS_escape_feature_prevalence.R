# S3: Associations of MSI status and TMB group with immune escape feature prevalence
# Within each TCGA cancer type, this script tests whether individual immune escape
# features are enriched in MSI versus MSS tumors or in high- versus low-TMB tumors.
# Associations are evaluated using Fisher's exact test and summarized with FDR,
# odds-ratio dot plots, and sample-count bar plots.
library(magrittr); library(cowplot); library(ggplot2); library(patchwork); library(writexl); library(RColorBrewer)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_mutation_rate_affect_escape_feature_prevalence/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Binarized immune escape profiles.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/DatasetLabels.R')
ID.xteam="PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/CombineData.XTeam.R")
profile.escape.binarization=CombineData.XTeam(ID.xteams, file.query='Results/BioImmune/50.ImmuneEscape/1.ImmuneEscapeProfiling/ProfilingImmuneEvading.Binarization.rds')
profile.B2M.inactive=CombineData.XTeam(ID.xteams, file.query="Results/BioGenomics/[Question]/[Q]B2M_biallelic_inactivation_status/B2M_biallelic_inactivation_matrix.rds")
# Primary tumor samples only.
source('/pub5/xiaoyun/BioY/sunshangqin/Functions/DataPrepare/GetInfor.PatientCenter.R')
patient.center=readRDS('/IData/DataCenter/TCGA/PanCancer_TCGA/PatientCenter/PatientCenter.rds')
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

# MSI subtype profiles for CRC, STAD, and UCEC.
ID.xteams.MSI=c('CRC_TCGA', 'STAD_TCGA', 'UCEC_TCGA')
patient.center=CombineData.XTeam(ID.xteams, file.query='PatientCenter/PatientCenter.rds')
profile.MSI=lapply(ID.xteams.MSI, function(x){
    GetInfor.PatientCenter(patient.center[[x]], SampleID=profile.escape.binarization[[x]]$SampleID, colNames="MSISubtype2")
})      %>%     setNames(ID.xteams.MSI)

# TMB profiles.
profile.TMB=lapply(ID.xteams, function(x){
    GetInfor.PatientCenter(patient.center[[x]], SampleID=profile.escape.binarization[[x]]$SampleID, colNames="TMB")
})      %>%     setNames(ID.xteams)

# Immune escape feature groups.
all.features.list=list(
    APMalt=c("HLA.LOH", "biallelic.B2M.inactivation", "HLA.mut", "B2M.mut", "all.APM.mut"),
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
names(all.features)=rep(names(all.features.list), lengths(all.features.list))

# MSI-H is classified as MSI; MSI-L is combined with MSS.
plot.data.MSI=lapply(ID.xteams.MSI, function(x){
    data=profile.escape.binarization[[x]]   %>%
        dplyr::mutate(MSI.Subtype=profile.MSI[[x]]$MSISubtype2,
            MSI.Subtype=ifelse(MSI.Subtype=='MSI-H', 'MSI', ifelse(MSI.Subtype %in% c('MSS', 'MSI-L'), 'MSS', NA)))       %>%  
        dplyr::filter(!is.na(MSI.Subtype))  %>%
        dplyr::mutate(MSI.Subtype=factor(MSI.Subtype, levels=c("MSI", "MSS")))

    tmp.data=lapply(all.features, function(m){
        df=na.omit(data[, c('MSI.Subtype', m)])
        df[, m]=factor(df[, m], levels=c("TRUE", "FALSE"))
        contingency_table=table(df$MSI.Subtype, df[, m])

        # Require more than 50 samples in each MSI group.
        if(length(contingency_table)>3 & min(table(df$MSI.Subtype))>50){       
            source("/pub5/xiaoyun/BioY/sunshangqin/Functions/BioStat/AssociationBetweenTwoCategoricalVariables.FisherExactTest.R")
            tmp=BinaryVarFisher(df, 'MSI.Subtype', m)
            result=data.frame(CancerType=gsub('_TCGA', '', x), escape.feature=m, MSI.sample=table(data$MSI.Subtype)['MSI'], MSS.sample=table(data$MSI.Subtype)['MSS'], tmp)
            return(result)
        }else{
            return(NULL)
        }
    })      %>%     do.call(what=rbind)
    tmp.data$FDR=p.adjust(tmp.data$p_value, method="BH")
    return(tmp.data)
})          %>%     do.call(what=rbind)
plot.data.MSI$escape.feature=factor(plot.data.MSI$escape.feature, levels=intersect(all.features, plot.data.MSI$escape.feature))
plot.data.MSI$feature.group=factor(names(all.features)[match(plot.data.MSI$escape.feature, all.features)], unique(names(all.features)))
plot.data.MSI$significance=factor(ifelse(plot.data.MSI$FDR<0.01, "FDR<0.01", ifelse(plot.data.MSI$FDR<0.05, "FDR<0.05", "no.sig")), levels=c("FDR<0.01", "FDR<0.05", "no.sig"))
plot.data.MSI$Odds_Ratio=factor(ifelse(plot.data.MSI$odds_ratio>1, "OR>1", "OR<1"), levels=c("OR>1", "OR<1"))

# MSI/MSS association and sample-count plots.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/DotPlot.R")
pdf(file.path(Dir.output, "2.1.association_between_MSI_subtype_and_escape_feature_prevalence.pdf"), width=15, height=4)
p1=dot_plot(t.data.frame=plot.data.MSI, x="escape.feature", y="CancerType", size="significance", 
        size.value=c("FDR<0.01"=4, "FDR<0.05"=3, "no.sig"=0.5),
        color="Odds_Ratio", facet="cols", group="feature.group")

bar.data=tidyr::pivot_longer(plot.data.MSI, cols=MSI.sample:MSS.sample,  names_to="MSI.Subtype", values_to="Value")
bar.data=unique(bar.data[, c('CancerType', 'Value', 'MSI.Subtype')])
p2=ggplot(data=bar.data, mapping=aes(x=CancerType, y=Value, fill=MSI.Subtype))+
            geom_bar(stat="identity", position=position_dodge(0.8), width=0.8)+
            scale_fill_manual(values=c("MSI.sample"="#2ec4b6", "MSS.sample"="#cbf3f0"))+
            labs(x='', y='samples number')+
            theme_classic()+
            ylim(c(0, ceiling(max(bar.data$Value))+100))+
            geom_text(aes(label=Value), size=3, position=position_dodge(width=0.9), vjust=0.25)+
            theme(legend.position='right', axis.text.y=element_blank())+
            coord_flip()
print(plot_grid(p1, p2, rel_widths=c(5, 1), align='h', axis='tb', ncol=2))
dev.off()

# High TMB is defined as TMB >= 10; low TMB is defined as TMB < 10.
plot.data.TMB=lapply(ID.xteams, function(x){
    data=profile.escape.binarization[[x]]       %>%
        dplyr::mutate(TMB=profile.TMB[[x]]$TMB, TMB=ifelse(TMB<10, 'low.TMB', 'high.TMB'))     %>%
        dplyr::filter(!is.na(TMB))      %>%
        dplyr::mutate(TMB=factor(TMB, levels=c('high.TMB', 'low.TMB')))

    tmp.data=lapply(all.features, function(m){
        df=data[, c('TMB', m)]
        df[, m]=factor(df[, m], levels=c('TRUE', 'FALSE'))
        contingency_table=table(df$TMB, df[, m])

        # Require more than 20 samples in each TMB group.
        if(length(contingency_table)>3 & min(table(df$TMB))>20){      
            source("/pub5/xiaoyun/BioY/sunshangqin/Functions/BioStat/AssociationBetweenTwoCategoricalVariables.FisherExactTest.R")    
            tmp=BinaryVarFisher(df, 'TMB', m)
            result=data.frame(CancerType=gsub('_TCGA', '', x), escape.feature=m, highTMB=table(data$TMB)['high.TMB'], lowTMB=table(data$TMB)['low.TMB'], tmp)
            return(result)
        }else{
            return(NULL)
        }
    })      %>%     do.call(what=rbind)
    tmp.data$FDR=p.adjust(tmp.data$p_value, method="BH")
    return(tmp.data)
})          %>%     do.call(what=rbind)
plot.data.TMB$escape.feature=factor(plot.data.TMB$escape.feature, levels=intersect(all.features, plot.data.TMB$escape.feature))
plot.data.TMB$feature.group=factor(names(all.features)[match(plot.data.TMB$escape.feature, all.features)], unique(names(all.features)))
plot.data.TMB$significance=factor(ifelse(plot.data.TMB$FDR<0.01, "FDR<0.01", ifelse(plot.data.TMB$FDR<0.05, "FDR<0.05", "no.sig")), levels=c("FDR<0.01", "FDR<0.05", "no.sig"))
plot.data.TMB$Odds_Ratio=factor(ifelse(plot.data.TMB$odds_ratio>1, "OR>1", "OR<1"), levels=c("OR>1", "OR<1"))

# High/low-TMB association and sample-count plots.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/PlotFunction/DotPlot.R")
pdf(file.path(Dir.output, "2.2.association_between_TMB_group_and_escape_feature_prevalence.pdf"), width=15, height=5)
p1=dot_plot(t.data.frame=plot.data.TMB, x="escape.feature", y="CancerType", size="significance", 
        size.value=c("FDR<0.01"=4, "FDR<0.05"=3, "no.sig"=0.5),
        color="Odds_Ratio", facet="cols", group="feature.group")

bar.data=tidyr::pivot_longer(plot.data.TMB, cols=highTMB:lowTMB,  names_to="TMB", values_to="Value")
bar.data=unique(bar.data[, c('CancerType', 'Value', 'TMB')])
p2=ggplot(data=bar.data, mapping=aes(x=CancerType, y=Value, fill=TMB))+
            geom_bar(stat="identity", position=position_dodge(0.8), width=0.8)+
            scale_fill_manual(values=c("highTMB"="#48cae4", "lowTMB"="#caf0f8"))+
            labs(x='', y='samples number')+
            theme_classic()+
            ylim(c(0, ceiling(max(bar.data$Value))+100))+
            geom_text(aes(label=Value), size=3, position=position_dodge(width=0.9), vjust=0.25)+
            theme(legend.position='right', axis.text.y=element_blank())+
            coord_flip()
print(plot_grid(p1, p2, rel_widths=c(5, 1), align='h', axis='tb', ncol=2))
dev.off()
