# Calculate per-feature high or low cutoffs as mean plus or minus one or two SDs.
Cutoff.ContinuousFeature.MeanSD <- function(cont.data.frame, 
        cutoff.type=c('high.group', 'low.group')[1],
		method=c("MeanSD", "Mean2SD")[1]
){
        # Estimate each feature independently and ignore missing values.
        feature.mean=colMeans(cont.data.frame, na.rm=T)
        feature.sd=apply(cont.data.frame, 2, sd, na.rm=T)
        
        if(cutoff.type=='high.group'){
            tmp.sd=ifelse(rep(method, length(feature.sd))=="MeanSD", feature.sd, 2*feature.sd)       
            group.cutoff=feature.mean + tmp.sd
        }else{
            tmp.sd=ifelse(rep(method, length(feature.sd))=="MeanSD", feature.sd, 2*feature.sd)
            group.cutoff=feature.mean - tmp.sd
        }
        
        return(group.cutoff)
}
