# Extract expression of predefined checkpoint genes; missing genes are retained as NA.
ImmuneCheckpoint.1 <- function(exprs.data
) {
	if(!is.null(exprs.data)){
        genes <- c('PDCD1','CD274','PDCD1LG2','CTLA4','FGL1','LAG3','BTLA','TIGIT','HAVCR2','CD47','ENTPD1','NT5E')

        check.exprs=exprs.data[match(genes, rownames(exprs.data)), , drop = TRUE]
		check.exprs=data.frame(SampleID=colnames(exprs.data), t(check.exprs))

colnames(check.exprs)=c('SampleID',genes)

	}else{
		check.exprs=NULL
	}
    return(check.exprs)
}
