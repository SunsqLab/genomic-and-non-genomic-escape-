# S1: Pan-cancer distribution of genomic and non-genomic immune escape
# In RStudio, set the working directory to this file's directory, edit the paths,
# then Source this file or run it from top to bottom. The functions are sourced below.
# From another directory: source("path/to/01_pan_cancer_escape_prevalence.R", chdir = TRUE).
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
library(UpSetR)

Dir.output="results"
if(!dir.exists(Dir.output)){ dir.create(Dir.output, recursive = TRUE) }
Dir.output <- normalizePath(Dir.output, mustWork = TRUE)

profile.escape <- readRDS(file.path(PanCancer.path, escape.query, "ProfilingImmuneEvading.rds"))

source(file.path(functions.dir, "GetInfor.PatientCenter.R"))

patient.center = readRDS(file.path(PanCancer.path, "PatientCenter", "PatientCenter.rds"))

# Primary samples only; keep the original nested PatientCenter input.
profile.escape = lapply(profile.escape, function(x){
            x[GetInfor.PatientCenter(patient.center, SampleID=x$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
        })

profile.escape.df=Reduce(function(x, y) dplyr::full_join(x, y, by = "SampleID"), profile.escape)
profile.escape.df=Filter(function(x) !all(is.na(x)), profile.escape.df)

source(file.path(functions.dir, "DatasetLabels.R"))
ID.xteam = "PanCancer_TCGA.dataset"
ID.xteams=get(ID.xteam)
ID.xteams=c(setdiff(ID.xteams, c('FPPP_TCGA', 'LAML_TCGA')))

source(file.path(functions.dir, "CombineData.XTeam.R"))
profile.escape.binarization <- CombineData.XTeam(ID.xteams, file.query = file.path(escape.query, "ProfilingImmuneEvading.Binarization.rds"))

profile.escape.binarization = lapply(profile.escape.binarization, function(x){
            tmp=lapply(x, function(y){
                        y[GetInfor.PatientCenter(patient.center, SampleID=y$SampleID, colNames="SampleType")$SampleType %in% "Primary", ]
                    })
            Reduce(function(x, y) dplyr::full_join(x, y, by = "SampleID"), tmp)
        })      %>%     setNames(ID.xteams)
B2M.loh <- CombineData.XTeam(ID.xteams, file.query = b2m.query)
# Replace matched B2M.LOH values with the biallelic-inactivation calls.
profile.escape.binarization <- mapply(function(df, b2m.mat) {
    if(!is.null(b2m.mat) && "B2M" %in% rownames(b2m.mat)) {
        idx <- match(df$SampleID, colnames(b2m.mat))
        matched <- !is.na(idx)
        df$B2M.LOH[matched] <- b2m.mat["B2M", idx[matched]] == 1
    }
    df
}, profile.escape.binarization, B2M.loh[names(profile.escape.binarization)], SIMPLIFY = FALSE)
profile.escape.binarization$PanCancer_TCGA = do.call(rbind, profile.escape.binarization)

pan.mut <- readRDS(file.path(PanCancer.path, "OMICSData", "Mutations.data.rds"))
pan.mut <- pan.mut[pan.mut$SampleID %in% profile.escape.df$SampleID,]

source(file.path(functions.dir, "ConstructMutationMatrix.R"))
IFNG.gene <- c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1")
gene.mut <- ConstructMutationMatrix(mutation.data = pan.mut, mut.types = "nonsilent")
res <- sapply(IFNG.gene, function(g) if(g %in% rownames(gene.mut)) gene.mut[g,] & !is.na(gene.mut[g,]) else rep(FALSE, ncol(gene.mut)))
IFNG.pathway.mut <- data.frame(SampleID = rownames(res), setNames(as.data.frame(res), paste0(IFNG.gene, ".mut")))

profile.escape.binarization <- lapply(profile.escape.binarization, function(data){

    data[colnames(IFNG.pathway.mut)[-1]] <- IFNG.pathway.mut[match(data$SampleID, IFNG.pathway.mut$SampleID), -1]
    return(data)
})

# Original feature groups and framework mappings.
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

data <- profile.escape.binarization$PanCancer_TCGA

APM.gene <- setdiff(names(data)[grep(".mut", names(data))], c("IFNG.pathway.mut", "CD58.mut", "all.APM.mut", "IDH1.mut", "HLA.mut", paste(IFNG.gene,"mut",sep=".")))
data$all.APM.mut <- rowSums(data[, APM.gene], na.rm = T)
data$HLA.mut <- rowSums(data[, c("HLA-A.mut", "HLA-B.mut", "HLA-C.mut")], na.rm = TRUE)
data$IFNG.pathway.mut <- rowSums(data[, paste0(IFNG.gene, ".mut")], na.rm = TRUE)

# Original within-sample sums: missing component values are ignored.
profile.escape.freq.df <- data %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(

        APMalt = sum(HLA.LOH, B2M.LOH, all.APM.mut, na.rm = TRUE),
        Checkalt = sum(CD274.deepAmp, na.rm = TRUE),
        ActMalt = sum(IFNG.pathway.HD, CD58.HD, IFNG.pathway.mut, IDH1.mut, CD58.mut, na.rm = TRUE),
        Checkexp = sum(CD274, CTLA4, PDCD1LG2, PDCD1, FGL1, LAG3, BTLA, TIGIT, HAVCR2, CD47, ENTPD1, NT5E, na.rm = TRUE),
        Supprcell = sum(M2_Macrophage, Treg, MDSC, Exhaust_CD8_Tcell, Cancer_Associated_Fibroblast, na.rm = TRUE),
        Supprsig = sum(TGFB1, Immune_supress_cytokine, SERPINB9, PTGER2, PTGER4, CXCL12, VEGFA, CD36, SLC43A2, na.rm = TRUE),
        Actdown = sum(CXCL9, CXCL10, CXCL11, CCL4, CCL5, CGAS, STING1, na.rm = TRUE),
        APMdown = sum(HLA.A, HLA.B, HLA.C, CALR, na.rm = TRUE),
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

temp.features <- features
temp.features$Custom.Framework.2 <- c(temp.features$Custom.Framework.2, "All.feature")

# NA becomes 0; the prevalence denominator is the full cohort, as in the source.
profile.escape.freq.group.df=tidyr::pivot_longer(profile.escape.freq.df[, unlist(temp.features)], cols = everything(), names_to = "escape.feature", values_to = "Group")        %>%
    dplyr::mutate(Group = ifelse(is.na(Group), 0, ifelse(Group > 2, ">2", Group))   )        %>%
    dplyr::count(escape.feature, Group, name = "No.samples")      %>%
    dplyr::filter(Group!=0)     %>%
    dplyr::mutate(Ratio.samples=signif(No.samples/nrow(profile.escape.freq.df), 2),
        escape.feature=factor(escape.feature, levels=unlist(temp.features)),
        Group=factor(Group, levels=c(">2", "2", "1")),
        freq=paste0(No.samples, ' (', signif(No.samples/nrow(profile.escape.freq.df), 2), ')'))%>%
        dplyr::arrange(escape.feature, desc(Group)) %>%
        as.data.frame()

framework_map <- stack(temp.features)
profile.escape.freq.group.df$Framework <- factor(framework_map$ind[match(profile.escape.freq.group.df$escape.feature, framework_map$values)], levels = names(temp.features))

data=profile.escape.freq.df
data=data[, grepl('.mut', names(data)) & !names(data)%in%c('HLA.mut', 'all.APM.mut', 'IFNG.pathway.mut', "Escape.mut", "CD58.mut", "IDH1.mut")]
profile.escape.gene.mut.df=data.frame(mutation=colnames(data), No.samples=colSums(data, na.rm=T), Ratio.samples=signif(colSums(data, na.rm=T)/nrow(data), 2))
profile.escape.gene.mut.df <- profile.escape.gene.mut.df %>%
    dplyr::mutate(mutation = gsub("\\.mut$","",mutation),
        class = ifelse(mutation %in% IFNG.gene, "IFNG.pathway.mut", "all.APM.mut"),
    )

trans.features.list <- features[names(features) %in% c("Checkexp", "Supprcell", "Supprsig", "Actdown", "APMdown", "Neo")]
trans.features <- unlist(trans.features.list, use.names = FALSE)
names(trans.features) <- names(all.features)[match(trans.features, all.features)]

# Score distributions and prevalence plots, in the original order.
pdf(file.path(Dir.output, '1.escape_score_density.pan_cancer.pdf'), 22, 10)
plots.p=lapply(intersect(trans.features, colnames(profile.escape.df)), function(y){
    tmp.density=max(density(profile.escape.df[,y], na.rm=T)$y)
    p=ggplot(profile.escape.df, aes(x =profile.escape.df[,y]))+
        geom_boxplot(width = tmp.density*0.1, outlier.colour = "white")+
        geom_density(fill=names(trans.features)[match(y, trans.features)])+
        coord_flip()+xlab('')+ggtitle(paste('PanCancer', y,sep='\n'))+ theme_bw()+ylab('')+
        theme(plot.title = element_text(size = 10), axis.text = element_text(size=10), legend.position='none')
})
print(plot_grid(plotlist=plots.p, nrow=4))
dev.off()

data <- profile.escape.freq.group.df[!profile.escape.freq.group.df$Framework %in% c("Custom.Framework", "Custom.Framework.2","Galassi.2024", "Roerden.2025"), ]

p = ggplot(data, aes(x = escape.feature, y = Ratio.samples)) +
        geom_bar(stat = "identity", width = 0.5, aes(fill = Group)) +
        geom_text(aes(label = paste0(Ratio.samples, "\n", "(", No.samples, ")")), size = 2, position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values = c("0" = "#fdf9f7", "1" = "#ffcad4", "2" = "#ff9ebb", ">2" = "#b9375e")) +
        ggtitle("PanCancer") +
        theme_bw() +
        facet_grid(. ~ Framework, scales = "free_x", space = "free_x") +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.text.y=element_text(size=10))

ggsave(plot = p, file.path(Dir.output, "2.escape_feature_prevalence.pan_cancer.pdf"), width = 15, height = 7)

temp.features <- features[names(features) %in% c("Custom.Framework", "Galassi.2024", "Roerden.2025")]
temp.features$other <- list(c("HLA.LOH", "B2M.LOH", "HLA.mut", "B2M.mut", "all.APM.mut", "CD274.deepAmp", "IFNG.pathway.HD", "CD58.HD", "IFNG.pathway.mut", "IDH1.mut", "CD58.mut", "APMdown", "Checkexp", "Supprcell", "Supprsig", "Actdown", "Neo"))

pdf(file.path(Dir.output, paste0("3.escape_feature_combinations.pan_cancer.pdf")), width = 10, height = 6)

# Top 20 exact intersections, not inclusive overlaps.
for (fw in 1:length(temp.features)) {
fea <- unlist(temp.features[[fw]])

data=profile.escape.freq.df[, fea]
data[which(data!=0, arr.ind=T)]=1
data[is.na(data)] = 0

p=upset(data, mb.ratio=c(0.6, 0.4), nsets=16,
        point.size = 3.2, keep.order = TRUE,  order.by = "freq",
        sets=rev(fea),
        nintersects = 20,
        sets.x.label = "Numbers of samples",
        sets.bar.color=rev(names(all.features)[match(fea, all.features)]))
print(p)
grid::grid.text(paste0('PanCancer (n=', nrow(data), ')'), x = 0.15, y = 0.95, gp = grid::gpar(fontsize = 15, fontface = "bold"))
}
dev.off()

data <- profile.escape.gene.mut.df
p=ggplot(data, aes(x =reorder(mutation, -Ratio.samples), y = Ratio.samples, fill = class))+
    geom_bar(stat ="identity",width = 0.6)+
    guides(fill = guide_legend(reverse = F))+
    labs(x='', y='Ratio of samples')+theme_classic()+
    theme(legend.position = "inside",legend.position.inside=c(0.9, 0.8), axis.text.x=element_text(angle=90, hjust=1, size=10))+
    scale_fill_manual(values = c("#C75D5D", "#4C78A8"))
ggsave(file.path(Dir.output, '4.escape_gene_mutation_frequency.pan_cancer.pdf'), p, width=18, height=4)
