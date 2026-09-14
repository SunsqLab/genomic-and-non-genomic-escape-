# Original complete-case Fisher workflow: select the one-sided tail by raw odds ratio, then BH.
source(file.path(functions.dir, "mutual_exclusivity.R"))

EscapePairFisher <- function(data, escape.feature){

        stat.data=lapply(names(data), function(cancer){
            df = data[[cancer]]

            mat=as.matrix(t(na.omit(df[, which(colnames(df)%in%escape.feature)])))
            fisher.mat=mat[rowSums(mat)!=0, ]

            if (!is.null(nrow(fisher.mat)) && nrow(fisher.mat)>1){
                df.fisher1=odds_ratio_matrix_fun(fisher.mat, alternative="less", cut_off = 1, p = 2, fdr =2)
                df.fisher2=odds_ratio_matrix_fun(fisher.mat, alternative="greater", cut_off = 1, p = 2, fdr = 2)
                fisher.re.df=rbind(df.fisher1, df.fisher2)

                df=fisher.re.df[, c('Neither', 'A_not_B', 'B_not_A', 'Both')]
                pos=which(df==0, arr.ind=T)[, 1]
                Hald.OR=apply(df[pos, ], 1, function(tmp.df){
                    tmp.df=tmp.df+0.5
                    result=(tmp.df[1]*tmp.df[4])/(tmp.df[2]*tmp.df[3])
                })
                fisher.re.df$odds_ratio[pos]=Hald.OR

                fisher.re.df$p_adj = p.adjust(fisher.re.df$p_value, method = "BH")
            }else{
                fisher.re.df=data.frame(geneA=rownames(fisher.mat), geneB=rownames(fisher.mat), odds_ratio=NA, p_value=NA, p_adj=NA)
            }
            rownames(fisher.re.df)=paste(fisher.re.df$geneA, fisher.re.df$geneB, sep='_')
            fisher.re.df$cancer =cancer
            return(fisher.re.df)
        })        %>%        do.call(what=rbind)

}
