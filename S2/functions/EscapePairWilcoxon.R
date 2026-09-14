# Original pairwise Wilcoxon tests; groups >= 5 and median ratio stabilized by exp(-6).
EscapePairWilcoxon <- function(data, disc.feature, cont.feature){

        stat.data=lapply(names(data), function(cancer){
                df=data[[cancer]]
                df=Filter(function(m) !all(is.na(m)), df)
                tmp.genome=intersect(colnames(df), disc.feature)
                tmp.trans=intersect(colnames(df), cont.feature)

                wilcox.re.df=lapply(tmp.genome, function(disc_var){
                    tmp.wilcox=lapply(tmp.trans, function(cont_var){
                                wilcox.mat=na.omit(df[, c(disc_var, cont_var)])

                                tmp=tapply(wilcox.mat[,cont_var], wilcox.mat[,disc_var], median)
                                wilcox.re=data.frame(feaA=disc_var, feaB=cont_var, fold_change=NA, p_value=NA)

                                if(length(tmp)==2 & min(table(wilcox.mat[,disc_var]))>4){
                                    p_value=wilcox.test(wilcox.mat[,cont_var] ~ wilcox.mat[, disc_var], exact=F)$p.value

                                    epsilon=exp(-6)
                                    wilcox.re=data.frame(feaA=disc_var, feaB=cont_var, fold_change=(tmp['TRUE']+epsilon)/(tmp['FALSE']+epsilon), p_value=p_value)
                                }
                                return(wilcox.re)
                            })        %>%        do.call(what=rbind)
                })        %>%        do.call(what=rbind)

                wilcox.re.df$p_adj = p.adjust(wilcox.re.df$p_value, method = "BH")
                wilcox.re.df$cancer = cancer
                return(wilcox.re.df)
        })      %>%     do.call(what=rbind)
        stat.data$label=paste(stat.data$feaA, stat.data$feaB, sep=':')

        return(stat.data)
 }
