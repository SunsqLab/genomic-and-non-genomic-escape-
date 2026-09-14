# Combine checkpoint-gene expression with CD274 high-level amplification status.
ImmuneCheckpoint <- function(exprs.data=NULL,
		CN.data=NULL,
		CN.Matrix.data=NULL,
		RefGenome=c('hg19','hg38')
){

source("immune_checkpoint_expression.R")
    if(!is.null(exprs.data)){

        check.exprs<- ImmuneCheckpoint.1(exprs.data)

	}else{
        check.exprs=NULL
    }

    source(file.path("functions", "single_gene_copy_number_events.R"))

# Use raw segments when available; otherwise use the discrete gene-level matrix.
	if(!is.null(CN.data)){
		CD274.amp=CNVGeneAmpDel(CN.data=CN.data, genes='CD274', RefGenome=RefGenome)
	}else if(!is.null(CN.Matrix.data)){
		tmp = CN.Matrix.data[intersect(rownames(CN.Matrix.data), 'CD274'), , drop=FALSE]
		CD274.amp=data.frame(SampleID=colnames(CN.Matrix.data), CD274.deepAmp=colSums(tmp == 2) > 0)
	}else{
		CD274.amp=NULL
	}

	tmp=list(check.exprs, CD274.amp)
	result=Reduce(dplyr::full_join, tmp[lengths(tmp)!=0])
	return(result)
}
