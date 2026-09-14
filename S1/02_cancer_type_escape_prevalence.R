# S1: Pan-cancer distribution of genomic and non-genomic immune escape
# In RStudio, set the working directory to this file's directory, edit the paths,
# then Source this file or run it from top to bottom. The functions are sourced below.
# From another directory: source("path/to/02_cancer_type_escape_prevalence.R", chdir = TRUE).
# Inputs keep the ORIGINAL RDS objects. No new TSV schema or preprocessing is required.
# DataCenter.rds must contain ID and Path; Path points to each original cohort folder.
# Edit the two query paths to match your existing folders; do not move the data.
DataCenter.file <- "/path/to/DataCenter.rds"
PanCancer.path <- "/path/to/PanCancer_TCGA.dataset"
escape.query <- "Results/BioImmune/ImmuneEscape"
b2m.query <- "Results/BioGenomics/B2M_biallelic_inactivation_matrix.rds"
functions.dir <- normalizePath("functions", mustWork = TRUE)

library(magrittr)
library(ggplot2)
library(cowplot)

Dir.output="results"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive = TRUE) }
Dir.output <- normalizePath(Dir.output, mustWork = TRUE)

source(file.path(functions.dir, "GetInfor.PatientCenter.R"))

patient.center = readRDS(file.path(PanCancer.path, "PatientCenter", "PatientCenter.rds"))

source(file.path(functions.dir, "DatasetLabels.R"))
ID.xteam = "PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))
source(file.path(functions.dir, "CombineData.XTeam.R"))
profile.escape.binarization = CombineData.XTeam(ID.xteams, file.query=file.path(escape.query, "ProfilingImmuneEvading.Binarization.rds"))

# Primary samples only; keep the original nested input tables.
profile.escape.binarization = lapply(profile.escape.binarization, function(x){
            tmp=lapply(x, function(y){
                        y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
                    })
            Reduce(function(x, y) dplyr::full_join(x, y, by = "SampleID"), tmp)
        })      %>%     setNames(ID.xteams)
B2M.loh <- CombineData.XTeam(ID.xteams, file.query = b2m.query)
profile.escape.binarization <- mapply(function(df, b2m.mat) {
    if(!is.null(b2m.mat) && "B2M" %in% rownames(b2m.mat)) {
        idx <- match(df$SampleID, colnames(b2m.mat))
        matched <- !is.na(idx)
        df$B2M.LOH[matched] <- b2m.mat["B2M", idx[matched]] == 1
    }
    df
}, profile.escape.binarization, B2M.loh[names(profile.escape.binarization)], SIMPLIFY = FALSE)

pan.mut <- readRDS(file.path(PanCancer.path, "OMICSData", "Mutations.data.rds"))

source(file.path(functions.dir, "ConstructMutationMatrix.R"))
IFNG.gene <- c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")
gene.mut <- ConstructMutationMatrix(mutation.data = pan.mut, mut.types = "nonsilent")
res <- sapply(IFNG.gene, function(g) if(g %in% rownames(gene.mut)) gene.mut[g,] & !is.na(gene.mut[g,]) else rep(FALSE, ncol(gene.mut)))
IFNG.pathway.mut <- data.frame(SampleID = rownames(res), setNames(as.data.frame(res), paste0(IFNG.gene, ".mut")))

profile.escape.binarization <- lapply(profile.escape.binarization, function(data){

    data[colnames(IFNG.pathway.mut)[-1]] <- IFNG.pathway.mut[match(data$SampleID, IFNG.pathway.mut$SampleID), -1]
    return(data)
})

profile.escape.Continuum <- CombineData.XTeam(ID.xteams, file.query=file.path(escape.query, "ProfilingImmuneEvading.rds"))

profile.escape.Continuum = lapply(profile.escape.Continuum, function(x) {
    tmp = lapply(x, function(y) {
        y[GetInfor.PatientCenter(patient.center, SampleID = y$SampleID, colNames = "SampleType")$SampleType %in% "Primary", ]
    })
    Reduce(function(x, y) dplyr::full_join(x, y, by = "SampleID"), tmp)
}) %>% setNames(ID.xteams)

profile.escape.Continuum <- dplyr::bind_rows(profile.escape.Continuum, .id = "CancerType")

# Compute the original sample-level sums separately within each cancer type.
profile.escape.freq <- lapply(profile.escape.binarization, function(data) {

    APM.gene <- setdiff(names(data)[grep(".mut", names(data))], c("IFNG.pathway.mut", "CD58.mut", "all.APM.mut", "IDH1.mut", "HLA.mut", paste(IFNG.gene,"mut",sep=".")))
    data$all.APM.mut <- rowSums(data[, APM.gene], na.rm = T)
    data$HLA.mut <- rowSums(data[, c("HLA-A.mut", "HLA-B.mut", "HLA-C.mut")], na.rm = TRUE)
    data$IFNG.pathway.mut <- rowSums(data[, paste0(IFNG.gene, ".mut")], na.rm = TRUE)

    tmp <- data %>%
        dplyr::group_by(SampleID) %>%
        dplyr::mutate(

        APMalt = sum(HLA.LOH, B2M.LOH, all.APM.mut, na.rm = TRUE),
        APMdown = sum(HLA.A, HLA.B, HLA.C, CALR, na.rm = TRUE),
        Checkalt = sum(CD274.deepAmp, na.rm = TRUE),
        Checkexp = sum(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm = TRUE),
        Supprcell = sum(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm = TRUE),
        ActMalt = sum(IFNG.pathway.HD, CD58.HD, IFNG.pathway.mut, IDH1.mut, CD58.mut, na.rm = TRUE),
        Supprsig = sum(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm = TRUE),
        Actdown = sum(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm = TRUE),
        Neo = sum(mean.neo.exprs, na.rm = TRUE),
        Genomic = sum(APMalt, Checkalt, ActMalt, na.rm = TRUE),
        NonGenomic = sum(APMdown, Checkexp, Supprcell, Supprsig, Actdown, Neo, na.rm = TRUE),
        All.feature = ifelse(rowSums(cbind(APMalt, APMdown, Checkalt, Checkexp, ActMalt, Supprcell, Supprsig, Actdown, Neo), na.rm = TRUE) != 0, 1, 0),

        Camouflage = sum(HLA.LOH, B2M.LOH, all.APM.mut, M2_Macrophage, Cancer_Associated_Fibroblast, TGFB1, CXCL9, CXCL10, CXCL11, CCL4, HLA.A, HLA.B, HLA.C, CALR, mean.neo.exprs, na.rm = TRUE),
        Coercion = sum(all.APM.mut, CD274.deepAmp, CD58.HD, IDH1.mut, CD58.mut, CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, TGFB1, Immune_supress_cytokine, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, CCL5, CGAS, STING1, na.rm = TRUE),
        Cytoprotection = sum(all.APM.mut, IFNG.pathway.HD, IFNG.pathway.mut, TGFB1, SERPINB9, na.rm = TRUE),

        Local = sum(HLA.LOH, B2M.LOH, all.APM.mut, CD274.deepAmp, IFNG.pathway.HD, CD58.HD, IFNG.pathway.mut, CD58.mut, CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, HLA.A, HLA.B, HLA.C, CALR, M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, TGFB1, Immune_supress_cytokine, SERPINB9, CD36, SLC43A2, na.rm = TRUE),
        Regional = sum(IDH1.mut, ENTPD1, NT5E, M2_Macrophage, MDSC, Cancer_Associated_Fibroblast, TGFB1, Immune_supress_cytokine, CXCL12, VEGFA, CXCL9, CXCL10, CXCL11, CCL5, CGAS, STING1, na.rm = TRUE),
        Global = sum(IFNG.pathway.HD, IFNG.pathway.mut, Treg, PTGER2, PTGER4, CD36, CCL4, CCL5, mean.neo.exprs, na.rm = TRUE)
    ) %>%
    as.data.frame()
}) %>% setNames(ID.xteams)

# Original feature groups and plotting order.
features = list(
    APMalt = c("HLA.LOH", "B2M.LOH", "HLA.mut", "B2M.mut", "all.APM.mut"),
    ActMalt = c("IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut"),
    Checkalt = c("CD274.deepAmp"),
    Checkexp = c("CD274", "CTLA4", "PDCD1LG2", "PDCD1", "FGL1", "LAG3", "BTLA", "TIGIT", "HAVCR2", "CD47", "ENTPD1", "NT5E"),
    Supprcell = c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast"),
    Supprsig = c("TGFB1", "Immune_supress_cytokine", "SERPINB9", "PTGER2", "PTGER4", "CXCL12", "VEGFA", "CD36", "SLC43A2"),
    Actdown = c("CXCL9", "CXCL10", "CXCL11", "CCL4", "CCL5", "CGAS", "STING1"),
    APMdown = c("HLA.A", "HLA.B", "HLA.C", "HLA.score", "CALR"),
    Neo = c("mean.neo.exprs")
)
features$Custom.Framework = names(features)
features$Custom.Framework.2 = c("Genomic", "NonGenomic")
features$Galassi.2024 = c("Camouflage", "Coercion", "Cytoprotection")
features$Roerden.2025 = c("Local", "Regional", "Global")

all.features <- unique(unlist(features))
names(all.features) <- c(rep(
    c("#c82621", "#fa8b69ff", "#F6C141", "#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"), lengths(features[1:9])),
    c("#c82621", "#fa8b69ff", "#F6C141", "#c0d666ff", "#97CC88", "#50AE94", "#8AC8E2", "#00b4d8", "#056795"),
    c("#C79DC9", "#CCD3E5"),
    rep("#DEDBEF",3),
    rep("#EDC3C8",3)
)

profile.escape.freq.group=lapply(profile.escape.freq, function(data){
    long_df=tidyr::pivot_longer(data[, all.features], cols = everything(), names_to = "escape.feature", values_to = "Group")        %>%
        dplyr::mutate(Group = ifelse(is.na(Group), 0, ifelse(Group > 2, ">2", Group)))        %>%
        dplyr::count(escape.feature, Group, name = "No.samples")      %>%
        dplyr::mutate(Ratio.samples=signif(No.samples/nrow(data), 2),
            escape.feature=factor(escape.feature, levels=all.features),
            Group=factor(Group, levels=c(">2", "2", "1", "0")),
            freq=paste0(No.samples, ' (', signif(No.samples/nrow(data), 2), ')'))
})      %>%     setNames(ID.xteams)

profile.escape.freq.group.df=lapply(ID.xteams, function(x){
            df=profile.escape.freq.group[[x]]     %>%
                    dplyr::select(-c(Ratio.samples, freq))        %>%
                    tidyr::pivot_wider(names_from = Group, values_from = No.samples, values_fill = 0)    %>%
                    dplyr::mutate(CancerType=x)     %>%
                    dplyr::filter(rowSums(select(., "1", "2", ">2"), na.rm=T)!=0)
        })        %>%        data.table::rbindlist(fill=TRUE)

profile.escape.freq.group.df$x=as.integer(factor(profile.escape.freq.group.df$CancerType, levels=names(profile.escape.freq)))-1
profile.escape.freq.group.df$y = as.integer(profile.escape.freq.group.df$escape.feature)

temp.features <- unlist(features,use.names = FALSE)
names(temp.features) <- rep(names(features), lengths(features))
profile.escape.freq.group.df$face_group <- names(temp.features)[match(profile.escape.freq.group.df$escape.feature, temp.features)]

data <- profile.escape.freq.group.df[profile.escape.freq.group.df$face_group %in% c("Custom.Framework", "Custom.Framework.2"), ]
data$face_group <- factor(data$face_group, levels = c("Custom.Framework", "Custom.Framework.2"))

Color=c("0"='#fdf9f7', "1"='#ffcad4', "2"='#ff9ebb', ">2"='#b9375e')
pie.cols <- names(Color)
# Original pie-polygon plotting helper, retained inline.
make_pie_polygons <- function(df, cols, r = 0.4, n = 80) {
    df <- as.data.frame(df)
    missing.cols <- setdiff(cols, names(df))
    if(length(missing.cols) > 0) df[missing.cols] <- 0

    pie.list <- lapply(seq_len(nrow(df)), function(i) {
        vals <- as.numeric(df[i, cols, drop = FALSE])
        vals[is.na(vals)] <- 0
        total <- sum(vals)
        if(total <= 0) return(NULL)

        props <- vals / total
        theta.end <- cumsum(props) * 2 * pi
        theta.start <- c(0, head(theta.end, -1))

        slice.list <- lapply(seq_along(cols), function(j) {
            if(vals[j] <= 0) return(NULL)
            theta <- seq(theta.start[j], theta.end[j], length.out = max(3, ceiling(n * props[j]) + 2))
            data.frame(
                pie.id = paste(i, cols[j], sep = "_"),
                Group = factor(cols[j], levels = cols),
                face_group = df$face_group[i],
                x = c(df$x[i], df$x[i] + r * cos(theta), df$x[i]),
                y = c(df$y[i], df$y[i] + r * sin(theta), df$y[i])
            )
        })
        do.call(rbind, slice.list)
    })
    do.call(rbind, pie.list)
}
pie.df <- make_pie_polygons(data, pie.cols, r = 0.4)
axis.df <- unique(as.data.frame(data)[, c("y", "escape.feature")])
p=ggplot() +
        geom_polygon(data = pie.df, aes(x = x, y = y, group = pie.id, fill = Group), color = "grey60", linewidth = 0.2)+
        theme_bw()+labs(x='', y='')+
        scale_fill_manual(name = "Group", values = Color, breaks = names(Color), limits = names(Color))+
        scale_x_continuous(breaks = 0:(length(ID.xteams)-1), labels = gsub('_TCGA', '', ID.xteams)) +
        scale_y_reverse(breaks = axis.df$y, labels = axis.df$escape.feature) +

        theme(axis.text.x = element_text(size = 12, hjust = 1, angle=45), strip.background.y = element_blank(), strip.text.y = element_blank(),
                axis.text.y = element_text(size = 12), legend.position="top")+
        facet_grid(face_group ~ ., scales = "free_y", space = "free_y")
ggsave(file.path(Dir.output, '5.escape_category_prevalence.by_cancer.pdf'), p, width=10, height=4.5)

pdf(file.path(Dir.output, '6.escape_feature_prevalence.by_cancer.pdf'), width=16, height=11)
plots.p=lapply(ID.xteams, function(x){
    p=ggplot(profile.escape.freq.group[[x]], aes(x=Group, y=Ratio.samples)) +
        geom_bar(stat="identity", fill="lightblue")+
        theme(legend.position = "none")+
        scale_fill_brewer(palette = 'Accent')+
        ggtitle(gsub('_TCGA', '', x))+
        facet_wrap(.~escape.feature, scales="free_x", nrow=4)+
        theme_bw()+
        geom_text(mapping = aes(label = freq), size=3.5)+
        theme(strip.text = element_text(size = 12))
    print(p)
})
dev.off()

temp.features <- c("APMalt", "ActMalt", "Checkalt", "Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo")
top.escape.combine <- c(
    "Checkexp, Supprcell, Supprsig", "APMalt, Checkexp, Supprcell, Supprsig", "Checkexp, Supprsig", "APMalt", "Supprsig", "Actdown, APMdown", "Supprsig, Actdown, APMdown",
    "Checkexp", "Actdown", "APMalt, Checkexp, Supprsig", "Checkexp, Supprsig, Actdown", "Supprsig, Actdown", "Checkexp, Supprcell", "APMalt, Actdown, APMdown",
    "APMalt, Supprsig, Actdown, APMdown", "Checkexp, Supprsig, Actdown, APMdown", "APMalt, Supprsig", "Checkexp, Supprcell, Supprsig, APMdown", "Supprcell, Supprsig", "Checkexp, Actdown"
)

prepare.combine.data <- function(data, features){
            data=as.data.frame(data)
            missing.features <- setdiff(features, colnames(data))
            if(length(missing.features)>0){ data[missing.features] <- 0 }
            data=data[, features, drop=FALSE]
            data=data.frame(lapply(data, as.numeric), check.names=FALSE)
            data[is.na(data)]=0
            data[data!=0]=1
            data
        }

profile.top.escape.freq.df=lapply(ID.xteams, function(x){
            data=prepare.combine.data(profile.escape.freq[[x]], temp.features)

            tmp=lapply(top.escape.combine, function(fea){
                        tmp.feature=unique(unlist(strsplit(fea, ', ')))

                        tmp1=data[which(rowSums(data, na.rm=T)==length(tmp.feature) & rowSums(data[, tmp.feature, drop=F], na.rm=T)==length(tmp.feature)),]
                        tmp2=data.frame(CancerType=gsub('_TCGA', '', x), combination=fea, count=nrow(tmp1), Ratio.escape=nrow(tmp1)/nrow(data), all.sample = nrow(data))
                    })      %>%     do.call(what=rbind)
        })      %>%     do.call(what=rbind)
profile.top.escape.freq.df$combination=factor(profile.top.escape.freq.df$combination, levels=rev(top.escape.combine))
profile.top.escape.freq.df$CancerType <- paste0(profile.top.escape.freq.df$CancerType,"(",profile.top.escape.freq.df$all.sample,")")

p=ggplot(profile.top.escape.freq.df, aes(x = CancerType, y = combination, fill = Ratio.escape)) +
        geom_tile(color = "grey90") +

        scale_fill_gradient(low = "white", high = "red") +
        scale_x_discrete() +
        theme_bw()+
        labs(fill='Proportion of samples with specific immune escape combinations', x='', y='')+
        coord_fixed()+
        theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1, size=10))

ggsave(file.path(Dir.output, "7.escape_feature_combinations.by_cancer.pdf"), p, width = 12, height = 7)

label <- c("Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo")
trans.features.list <- features[names(features) %in% label]
trans.features <- unlist(trans.features.list, use.names = FALSE)
names(trans.features) <- names(all.features)[match(trans.features, all.features)]

pdf(file.path(Dir.output, '1.escape_score_density.by_cancer.pdf'), width = 22, height = 10)

for(cancer in unique(profile.escape.Continuum$CancerType)){
    print(cancer)
    data <- profile.escape.Continuum[profile.escape.Continuum$CancerType == cancer, ]
    data <- as.data.frame(data)
    data <- data[, colSums(!is.na(data)) > 0, drop = FALSE]

    plots.p=lapply(intersect(trans.features, colnames(data)), function(y){
        plot.data <- data.frame(value=as.numeric(data[[y]]))
        plot.data <- plot.data[is.finite(plot.data$value) & !is.na(plot.data$value), , drop=FALSE]
        feature.color <- names(trans.features)[match(y, trans.features)]
        if(nrow(plot.data)>=2 && length(unique(plot.data$value))>1){
            tmp.density=max(density(plot.data$value, na.rm=T)$y)
            p=ggplot(plot.data, aes(x=value))+
                geom_boxplot(width = tmp.density*0.1, outlier.colour = "white")+
                geom_density(fill=feature.color)+
                coord_flip()+xlab('')+ggtitle(paste(cancer, y,sep='\n'))+ theme_bw()+ylab('')+
                theme(plot.title = element_text(size = 10, face="bold"), axis.text = element_text(size=10, face="bold"), legend.position='none')
        }else{
            p=ggplot(plot.data, aes(x=value))+
                geom_point(aes(y=1), color=feature.color, size=1.5, alpha=0.8)+
                xlab('')+ggtitle(paste(cancer, y,sep='\n'))+ theme_bw()+ylab('')+
                theme(plot.title = element_text(size = 10, face="bold"), axis.text = element_text(size=10, face="bold"), legend.position='none')
        }
        p
    })
    if(length(plots.p)==0){
        plots.p <- list(ggplot()+theme_void()+labs(title=cancer))
    }

    legend.labels <- unique(label)
    legend.colors <- vapply(legend.labels, function(x){
            names(all.features)[match(features[[x]][1], all.features)]
        }, character(1))
    legend.df <- data.frame(x=seq_along(legend.labels), y=rep(1, length(legend.labels)), label=legend.labels, color=legend.colors)
    plot_grid(plot_grid(plotlist=plots.p, nrow=4),
          ggplot(legend.df) +
            geom_point(aes(x=x, y=y, color=color), size=5, shape=15) +
            geom_text(aes(x=x, y=y, label=label), vjust=-1) +
            theme(plot.title = element_text(face="bold"), text=element_text(face="bold") ) +
            scale_color_identity() + theme_void() + labs(title="Legend"),
          nrow=2, rel_heights=c(4,1)) %>% print()
}
dev.off()
