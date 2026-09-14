# Calculate ssGSEA scores for M2 macrophages, Tregs, MDSCs, exhausted CD8 T cells, and CAFs.
ImmuneSuppreCell <- function(exprs.data){
	library(magrittr)
	
	M2_macrophage=c('HRH1', 'NPL', 'PDCD1LG2', 'RENBP', 'CFP', 'C1QC', 'C1QB', 'C1QA', 'CD5L', 'FCGR3A', 'ITGB5', 'MERTK', 'PILRA', 'CCL8')
	Treg=c('IL2RA', 'FOXP3', 'CTLA4', 'SLC35D1', 'GDPD3', 'CISH')
	
	MDSC_gene=c('CCR2','CD14','CD2','CD86','CXCR4','FCGR2A','FCGR2B','FCGR3A','FERMT3','GPSM3','IL18BP','IL4R','ITGAL','ITGAM','PARVG','PSAP','PTGER2','PTGES2','S100A8','S100A9')
	
	FOXP3='FOXP3'

CD8_Tcell_exhaust=c('FASLG','PBX3','CD244','CCL3','EOMES','CASP3','PLSCR1','MDFIC','CTLA4','PDCD1','IER5','RGS16','EEA1','TNFRSF9','PENK','COCH','PTPN13','NR4A2','CD160','PTGER4','CCL4','WBP5','GPR56','HSPC159','ENTPD1','SH2D2A','SEPT4','ISG20','TRIM47','CASP4','IFIH1','RBM39','LAG3','NFATC1','CA2','GAS2','MX1','ITIH5','GPD2','UBE2T','RNF11','CAPZB','TUBB2A','BUB1','JAK3','TCTA','CD9','TTC39B','RCN1','ROMO1','VAMP7','ETF1','CPA3','CD7','ART3','C14orf156','ATF1','WNK1','CIT','CCRL2','PLIN2','VPS37A','TCEA2','MYH4','TNFRSF1A','SPP1','S100A13','PON2','LCLAT1','ISG15','TANK','SHKBP1','RSAD2','TBCID22A','IRF8','GDPD5','GDF3','ITGAV','TMEM109','CPSF2','KLK6','CPT2','LMAN2','TOR3A','GPR65','MKI67','NPTXR','SNRPB2','NDFIP1','PTGER2','ZFP91','SPOCK2','C10orf58','CXCL10','SCIN','TRIM25','LAT2','CD200','DOCK7','PAWR','CHL1')

CAF_gene <- c("FAP", "COL1A1", "COL1A2", "DCN", "THY1", "COL6A1", "COL6A2", "COL6A3", "ACTA2", "COL3A1", "PDPN", "CXCL12", "PDGFRB", "POSTN", "CD34", "CFD", "FNDC1", "MMP11", "PDGFRA", "VIM", "WWTR1", "BGN", "FN1", "LUM", "MMP2", "RGS5", "SPARC")

library(GSVA)
	
    # Score each predefined gene set separately and return a sample-by-feature table.
    tmp.pheno=lapply(list(M2_macrophage, Treg, MDSC_gene, CD8_Tcell_exhaust, CAF_gene), function(pheno.gene) {
               tmp <- gsva(data.matrix(exprs.data), gset.idx.list = list(pheno.gene), method = "ssgsea")
            }) %>%      do.call(what = rbind) %>%    t()
    colnames(tmp.pheno) <- c("M2_Macrophage", "Treg", "MDSC", "Exhaust_CD8_Tcell", "Cancer_Associated_Fibroblast")

	# Add FOXP3 expression as a separate marker.
	tmp.FOXP3=as.numeric(exprs.data[FOXP3,])
	result=data.frame(SampleID=colnames(exprs.data), FOXP3=tmp.FOXP3, tmp.pheno)
	return(result)
}
