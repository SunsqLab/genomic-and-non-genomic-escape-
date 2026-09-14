# Original panel-wide complete-case Spearman correlations; BH includes both matrix directions.
library(Hmisc)

EscapePairSpearman <- function(data, escape.feature){

        stat.data=lapply(names(data), function(cancer){
            df = data[[cancer]]
            df=Filter(function(x) !all(is.na(x)), df)
            share.feature=intersect(escape.feature, names(df))
            spear.mat=as.matrix(na.omit(df[, match(share.feature, colnames(df))]))

            if(nrow(spear.mat)>4){
                rcorr_result=rcorr(spear.mat, type = "spearman")
                OR=rcorr_result$r

                spearman.re.df=as.data.frame(reshape2::melt(OR, varnames = c("feaA", "feaB"), value.name = "cor"))
                spearman.re.df$p_value=c(rcorr_result$P)
                spearman.re.df$p_adj=p.adjust(spearman.re.df$p_value, method = "BH")

                spearman.re.df=na.omit(spearman.re.df)
            }else{
                spearman.re.df=data.frame(feaA=colnames(spear.mat), feaB=colnames(spear.mat), cor=NA, p_value=NA)
            }
            spearman.re.df$cancer =cancer
            return(spearman.re.df)
        })      %>%     do.call(what=rbind)
        stat.data$label=paste(stat.data$feaA, stat.data$feaB, sep=':')

        return(stat.data)
 }
