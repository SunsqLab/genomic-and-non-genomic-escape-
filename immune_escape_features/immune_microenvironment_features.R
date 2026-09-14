# Assemble expression-based cell and cytokine features with genomic immune-factor features.
source("immune_factor_expression.R")
source("immune_factor_genomic_alterations.R")
source("immune_suppressive_cell_infiltration.R")
ImmuneSuppreEnvir <- function(exprs.data=NULL,
		mut.data=NULL,
		CN.data=NULL,
		CN.Matrix.data=NULL,
		RefGenome=c('hg19', 'hg38')
){
	library(dplyr)
	
	if(!is.null(exprs.data)){
		Suppre.Cell=ImmuneSuppreCell(exprs.data)
		Immune.Factor=ImmuneSuppreCytokine(exprs.data)
	}else{
		Suppre.Cell=NULL; Immune.Factor=NULL
	}

	if(!is.null(mut.data)){
		Factor.Genomic=ImmuneCytokineGenomic(mut.data=mut.data, CN.data = CN.data, 
                                  CN.Matrix.data = CN.Matrix.data, 
                                  RefGenome=RefGenome)
	}else{
		Factor.Genomic=NULL
	}

tmp=list(Suppre.Cell=Suppre.Cell, Immune.Factor=Immune.Factor, Factor.Genomic=Factor.Genomic)

return(tmp)

}
