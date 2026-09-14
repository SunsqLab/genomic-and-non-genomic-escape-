# S3: Immune escape prevalence across TMB intervals within high- and low-TMB groups
# This script separates primary TCGA tumors at TMB = 10 and creates 20 log10(TMB)
# intervals within each group. Analyses are run for all cancers and after excluding
# CRC, STAD, and UCEC. Feature trends are assessed using permutation-based
# Cochran-Armitage tests and displayed as prevalence curves.
library(magrittr); library(cowplot); library(ggplot2); library(scales); library(openxlsx); library(RColorBrewer)
Dir.output='/WorkSpace/sunshangqin/Immune_Escape/[Q]Factors_affecting_escape_feature_prevalence/Does_mutation_rate_affect_escape_feature_prevalence/Question1'
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive=TRUE) }

# Cancer-specific escape profiles are combined for the pan-cancer analysis.
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
})      %>%     do.call(what=rbind)

# Composite transcriptomic escape features.
profile.trans.escape=profile.escape.binarization 	%>%		
    dplyr::group_by(SampleID) %>%
    dplyr::reframe(
        Checkpoint.overexprs=any(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm=TRUE),		
        SuppreCell.overexprs=any(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm=TRUE),
        SuppreFactor.overexprs=any(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm=TRUE),
        ActivationGene.downexprs=any(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm=TRUE),
        HLA.downregulation=any(HLA.A, HLA.B, HLA.C, HLA.score, CALR, na.rm=TRUE))     %>%
    as.data.frame()

# Pan-cancer TMB profiles.
profile.TMB=GetInfor.PatientCenter(patient.center, colNames=c("SampleID", "TMB", "CancerType")) %>%
    dplyr::mutate(
        CancerType=ifelse(CancerType %in% c("READ", "COAD"), "CRC", CancerType),
        CancerType=paste0(CancerType, "_TCGA"))   %>%
    as.data.frame()

# Immune escape feature groups; HLA.mut and B2M.mut are represented by all.APM.mut.
all.features.list=list(
    APMalt=c("HLA.LOH", "biallelic.B2M.inactivation", "all.APM.mut"),
    Checkalt=c("CD274.deepAmp"),
    ActMalt=c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
    Checkexp=c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    Supprcell=c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    Supprsig=c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    Actdown=c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown=c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo=c("mean.neo.exprs")
)
all.features=unname(unlist(all.features.list))

# TMB intervals within high- and low-TMB strata.
ID.xteams.list=list(ID.xteams, setdiff(ID.xteams, c("CRC_TCGA", "STAD_TCGA", "UCEC_TCGA")))
TMB.group.data=lapply(ID.xteams.list, function(ID.xteams){
    samples=profile.TMB$SampleID[profile.TMB$CancerType %in% ID.xteams]
    data=profile.escape.binarization    %>%
        dplyr::filter(SampleID %in% samples) %>%
        dplyr::left_join(profile.TMB %>% select(SampleID, TMB), by="SampleID")  %>%
        dplyr::mutate(TMB.split=ifelse(TMB<10, "lowTMB", "highTMB"))
    split.data=split(data, data$TMB.split)

    tmp=lapply(names(split.data), function(df){
        TMB.bucket=20
        TMB.group.data=split.data[[df]]       %>%        
            dplyr::filter(!is.na(TMB))     %>%        
            dplyr::reframe(SampleID=SampleID, TMB=TMB, log10_TMB=log10(TMB))     %>%
            dplyr::filter((log10_TMB>quantile(log10_TMB, 0.01)) & (log10_TMB<quantile(log10_TMB, 0.99)))    %>%
            dplyr::mutate(TMB.group=cut(log10_TMB, breaks=TMB.bucket, labels=F))
        index=TMB.group.data      %>%    dplyr::group_by(TMB.group)    %>%
            dplyr::summarise(max.TMB=signif(max(TMB), 2))
        index$max.TMB[1]=paste0('<=', index$max.TMB[1])
        index$max.TMB[TMB.bucket]=paste0('>', index$max.TMB[TMB.bucket-1])
        TMB.group.data$TMB.group=factor(index$max.TMB[match(TMB.group.data$TMB.group, index$TMB.group)], levels=index$max.TMB)
        return(TMB.group.data)
    })      %>%     setNames(names(split.data))
})          %>%     setNames(c('allCancer', 'nonMSICancer'))

# Cochran-Armitage trend tests within each TMB stratum.
source("/pub5/xiaoyun/BioY/sunshangqin/Functions/BioStat/TrendTestBetweenBinaryAndOrdinalVariables/CochranArmitageTrendTestWithPermutation.R")
all.CA.test.result=lapply(names(TMB.group.data), function(tmp.cancer){
    data=TMB.group.data[[tmp.cancer]]
    tmp.CA.result=lapply(names(data), function(tmp.TMB){
        df=data[[tmp.TMB]]
        CA.data=Reduce(function(x, y)      dplyr::full_join(x, y, by='SampleID'), list(profile.escape.binarization, profile.trans.escape, df))
        CA.data=CA.data %>%
            dplyr::mutate(across(all_of(all.features), ~factor(.x, levels=c(TRUE, FALSE))))            %>%
            dplyr::filter(!is.na(TMB.group))
        CA.test.result=CochranArmitagePermu(CA.data, binary.vars=all.features, ordinal.var='TMB.group', is.pCorrect=FALSE)
        CA.test.result$FDR=p.adjust(CA.test.result$p_value, method="BH")
        CA.test.result$TMB.split=tmp.TMB
        return(CA.test.result)
    })      %>%     do.call(what=rbind)
    file_path=file.path(Dir.output, "Association_between_TMB_and_immune_escape_features.xlsx")
    wb=loadWorkbook(file_path)
    sheet.name=ifelse(tmp.cancer=='allCancer', 'Table S3', 'Table S4')
    if (sheet.name %in% names(wb)) {
        removeWorksheet(wb, sheet.name)
    }
    addWorksheet(wb, sheet.name)
    writeData(wb, sheet=sheet.name, x=paste0(sheet.name, ". Association between TMB and immune escape features (CA-test) in ", tmp.cancer), startRow=1, startCol=1)
    writeData(wb, sheet=sheet.name, tmp.CA.result, startRow=2, headerStyle=createStyle(textDecoration="bold"))
    saveWorkbook(wb, file=file_path, overwrite=TRUE)
    return(tmp.CA.result)
})  %>%     setNames(names(TMB.group.data))

# Prevalence curves and trend-test summaries.
escape.feature.colors=c(
    "HLA.LOH"='#22577a', "biallelic.B2M.inactivation"='#FBB684', "all.APM.mut"="#e71d36", 
    "IFNG.pathway.mut"="#da627d", "CD58.mut"="#ffa69e", "IFNG.pathway.HD"="#7bdff2", "CD58.HD"="#489fb5", "IDH1.mut"="#E78538", "CD274.deepAmp"="#9c6644",
    
    "ENTPD1"="#e9ecef", "FGL1"="#e3d5ca", "HAVCR2"="#e3f2fd", "CD47"="#c7f9cc", "CD274"="#cbf3f0", 

    "PDCD1"="#E48E91", "PDCD1LG2"="#7E94B7", "CTLA4"="#E78538", "LAG3"="#94C665", "BTLA"="#D8B591", 
    "TIGIT"="#277da1", "NT5E"="#e0aaff", 

    'Treg'="#93e1d8", 'M2_Macrophage'="#ffa69e", 'MDSC'="#74a9cf", "Exhaust_CD8_Tcell"="#a4c3b2", "Cancer_Associated_Fibroblast"="#ffb627", 
    'TGFB1'="#a1dab4", 'Immune_supress_cytokine'="#edae49", "SERPINB9"="#5390d9", "PTGER2"="#64dfdf", "PTGER4"="#ff8fa3", "CXCL12"="#52b788", "VEGFA"="#a2d6f9", "CD36"="#d1b3c4", "SLC43A2"="#f08080",
    
    "CXCL9"="#57cc99", "CXCL10"="#00afb9", "CXCL11"="#ffc43d", "CCL4"="#f7a072", "CCL5"="#c37d92", "CGAS"="#cad5ca", "STING1"="#c4fff9",
    'HLA.A'="#006d2c", 'HLA.B'="#2ca25f", 'HLA.C'="#66c2a4", 'HLA.score'="#b2e2e2", "CALR"="#b298dc", "mean.neo.exprs"="#ecab63",
    "Checkpoint.overexprs"='#a06cd5', "Novel.Check.overexprs"="#e4d9ff", "SuppreCell.overexprs"='#0077b6', "SuppreFactor.overexprs"='#00b4d8', 
    "ActivationGene.downexprs"="#9ceaef", "HLA.downregulation"="#32bfa4")

pdf(file.path(Dir.output, "3.3.escape_feature_prevalence_across_TMB_intervals_by_TMB_group.pan_cancer.pdf"), 7, 4)
result=lapply(names(TMB.group.data), function(tmp.cancer){
    tmp.data=TMB.group.data[[tmp.cancer]]
    tmp2=all.CA.test.result[[tmp.cancer]]
    tmp2=split(tmp2, tmp2$TMB.split)
    lapply(names(tmp.data), function(tmp.TMB){
        df=tmp.data[[tmp.TMB]]
        plots=lapply(all.features.list, function(features){
            data=Reduce(function(x, y)      dplyr::full_join(x, y, by='SampleID'), list(profile.escape.binarization, profile.trans.escape, df))
            df1=data[, c(features, 'TMB.group')]     %>%     
                    dplyr::filter(!is.na(TMB.group))        %>%     
                    dplyr::group_by(TMB.group)     %>%
                    dplyr::mutate(across(everything(), ~ sum(., na.rm=TRUE)/length(TMB.group)), TMB.group.sam=length(TMB.group))      %>%
                    unique()
            df2=reshape2::melt(df1, measure.vars=setdiff(colnames(df1), c("TMB.group", "TMB.group.sam")),
                    variable.name='features', value.name='escape.ratio')

            p1=ggplot(df1, aes(x=TMB.group, y=TMB.group.sam))+
                    geom_col(alpha=0.4)+labs(x='', title=paste0('PanCancer ', tmp.cancer, '_', tmp.TMB), y='samples')+
                    theme_test() +
                    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                        panel.border=element_rect(fill=NA, color="black", linewidth=0.3, linetype="solid"))  

            p2=ggplot(df2, aes(x=TMB.group, y=escape.ratio, color=features, group=features)) +
                    geom_point(size=1)+
                    geom_smooth(method="loess", aes(color=features), linewidth=0.6, alpha=0.3, se=FALSE) +     
                    scale_y_continuous(labels=label_number(accuracy=0.01)) +
                    labs(y='Prevalence of escape')+
                    scale_color_manual(values=escape.feature.colors, guide=guide_legend(nrow=2))+   
                    theme_test() +
                    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom", axis.line=element_line(linewidth=0.3))

            tmp1=df2[grep('>', df2$TMB.group), ]
            tmp2=tmp2[[tmp.TMB]]        %>%
                dplyr::mutate(features=binary.var, z.statistic, FDR)
            df3=dplyr::full_join(tmp1, tmp2, by="features")     %>%
                na.omit()   %>%
                dplyr::mutate(label=case_when(
                    z.statistic>0 & FDR<0.05 ~ "pos.cor",
                    z.statistic<0 & FDR<0.05 ~ "neg.cor", TRUE ~ ""))  %>%
                unique()
            ymin=min(df2$escape.ratio, na.rm=TRUE)
            ymax=max(df2$escape.ratio, na.rm=TRUE)
            p3=ggplot(df3, aes(x=1, y=escape.ratio, color=features, shape=label)) +
                geom_point(size=4) +
                scale_y_continuous(limits=c(ymin, ymax)) +
                scale_shape_manual(values=c("pos.cor"="+", "neg.cor"="-")) +
                scale_color_manual(values=escape.feature.colors) +
                theme_classic() +
                theme(
                    axis.ticks=element_blank(), axis.text=element_blank(), 
                    axis.title.x=element_blank(), axis.title.y=element_blank(), 
                    axis.line=element_blank(),
                    legend.position="none")

            top=plot_grid(p1, NULL, ncol=2, rel_widths=c(1, 0.05))
            bottom=plot_grid(p2, p3, ncol=2, align="h", rel_widths=c(1, 0.05))
            p=plot_grid(top, bottom, ncol=1, rel_heights=c(1, 2.5))
            print(p)
        })
    })
})
dev.off()
