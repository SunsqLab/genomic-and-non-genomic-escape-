# Combine genomic and expression-derived antigen-presentation features.
AntigenPresentDefect <- function(exprs.data=NULL,
		mut.data=NULL,
		CN.data=NULL,
		CN.Matrix.data=NULL,
		RefGenome=c('hg19','hg38')
){
	if(!(is.null(mut.data) & is.null(CN.data))){
        source("antigen_presentation_genomic_alterations.R")

		HLA.APM.alt=AntigenPresentDefectAlteration(mut.data, CN.data, CN.Matrix.data=CN.Matrix.data, RefGenome=RefGenome)
	}else{
		HLA.APM.alt=NULL
	}
	
	if(!is.null(exprs.data)){
        source("antigen_presentation_expression.R")
		HLA.exprs=AntigenPresentDefectHLAexpr(exprs.data)
	}else{
		HLA.exprs=NULL
	}

	tmp=list(HLA.APM.alt, HLA.exprs)
	result=Reduce(dplyr::full_join, tmp[lengths(tmp)!=0])
	return(result)
}
