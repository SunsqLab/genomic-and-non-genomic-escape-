# Retrieve protein-coding gene coordinates from the selected Ensembl annotation.
extract.coor.gene <- function(geneSet=NULL, genome.build="hg19"){
	library(IRanges)
	library(GenomicRanges)

if(genome.build=="hg19")
	   {
	    gene.infor <- get(load(file.path("reference_data", "ensembl_gene_annotation_hg19.RData")))
	   }
	if(genome.build=="hg38")
	   {
	     gene.infor <- get(load(file.path("reference_data", "ensembl_gene_annotation_hg38.RData")))
	   }

data <- gene.infor[gene.infor$type=="gene"&gene.infor$gene_biotype=="protein_coding"]

# Prefer a unique HAVANA record for duplicated gene symbols; otherwise retain the longest annotation.
	re.gene<-names(table(data$gene_name))[as.numeric(table(data$gene_name))>1]
	
	index.all<-which(data$gene_name%in%re.gene)
	
	index.left<-c()	
	for(i in re.gene){
		re.data<-data[data$gene_name==i]
		
		if(sum(re.data$source%in%"HAVANA")>0){
			
			if(sum(re.data$source=="HAVANA")==1)
			index<-which(data$gene_name==i&data$source=="HAVANA")
		
		else{
			re.width<-as.data.frame(ranges(re.data))$width
			index<-which(data$gene_name==i&as.data.frame(ranges(data))$width==max(re.width))
		}
		}
		
		else{
			re.width<-as.data.frame(ranges(re.data))$width
			index<-which(data$gene_name==i&as.data.frame(ranges(data))$width==max(re.width))
		}
		index.left<-c(index.left,index)
	}
	
	index.remove<-setdiff(index.all,index.left)
	newdata <- data[-index.remove, ]
	all.gene.chr.infor <- as.data.frame(newdata)
	
	extract.colnames <- c("seqnames","start","end","width","strand","source","gene_id","gene_version","gene_name")
	if(is.null(geneSet))
	  {
	    result <- all.gene.chr.infor[,extract.colnames]
	  }else{
	    result <- all.gene.chr.infor[which(all.gene.chr.infor[,"gene_name"] %in% geneSet),extract.colnames]
	  }
	  
	return(result)
	
}
