# Assess high-level SETDB1 amplification from copy-number data.
EpigeneticAlt <- function(CN.data,
		RefGenome=c('hg19','hg38')
){

source(file.path("functions", "single_gene_copy_number_events.R"))

	if(!is.null(CN.data)){
		SETDB1.amp=CNVGeneAmpDel(CN.data=CN.data, genes='SETDB1', RefGenome=RefGenome)
	}else{
		SETDB1.amp=NULL
	}

	return(SETDB1.amp)
}
