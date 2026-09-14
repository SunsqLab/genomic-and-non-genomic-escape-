# Binarize continuous escape features while preserving existing discrete features.
source(file.path("functions", "threshold_selection.R"))
EscapeFactor.Binarization <- function(escape.data
){

        library(magrittr)

# High suppressive expression and low immune-activating expression are treated as distinct events.
        if ("Immune.Factor" %in% names(escape.data) & !is.null(escape.data[["Immune.Factor"]])) {
            escape.data <- within(escape.data, {
                Immune.Factor.stimu <- Immune.Factor[, intersect(c(
                    "SampleID", "CXCL9", "CXCL10", "CXCL11",
                    "CCL4", "CCL5", "CGAS", "STING1"
                ), colnames(Immune.Factor))]
                Immune.Factor.suppr <- Immune.Factor[, c("SampleID", colnames(Immune.Factor)[!colnames(Immune.Factor) %in% colnames(Immune.Factor.stimu)])]
                rm(Immune.Factor)
            })
            all.features=c('Checkpoint.feature', 'Suppre.Cell', "Immune.Factor.stimu", "Immune.Factor.suppr", 'HLA.exp', 'NeoantigenExp.feature')
        }else{
            all.features=c('Checkpoint.feature', 'Suppre.Cell', "Immune.Factor", 'HLA.exp', 'NeoantigenExp.feature')
        }

tmp1.binar=lapply(all.features, function(feature){
            data=escape.data[[feature]]
			
			if(is.null(data)) return(NULL)
			pos=which(sapply(data, class)=='numeric')
			if(length(pos) == 0) return(NULL)
			cont.data.frame=data[, pos]

if("Immune.Factor.suppr" %in% names(escape.data) & !is.null(escape.data[["Immune.Factor.suppr"]])){
                cutoff.type=ifelse(feature%in%c('Checkpoint.feature', 'Suppre.Cell', 'Immune.Factor.suppr'), "high.group", "low.group")
            }else{
                cutoff.type=ifelse(feature%in%c('Checkpoint.feature', 'Suppre.Cell'), "high.group", "low.group") 
            } 

            # Apply a per-feature cohort-level mean plus or minus SD cutoff.
            cut.off=Cutoff.ContinuousFeature.MeanSD(cont.data.frame, cutoff.type=cutoff.type)

            if(cutoff.type=="high.group"){
                tmp.data=lapply(1:length(cut.off), function(y){
                    tmp=cont.data.frame[, y] >  cut.off[y]
                })      %>%     as.data.frame()
            }else{
                tmp.data=lapply(1:length(cut.off), function(y){
                    tmp=cont.data.frame[, y] <  cut.off[y]
                })      %>%     as.data.frame()
            }
            data[, pos]=tmp.data

            return(data)
        })      %>%     setNames(all.features)

# Convert mutation-count summaries to occurrence indicators.
        if("HLA.mut" %in% colnames(escape.data$APM.alt)){
            escape.data$APM.alt$HLA.mut=ifelse(escape.data$APM.alt$HLA.mut>0, TRUE, FALSE)
            escape.data$APM.alt$all.APM.mut=ifelse(escape.data$APM.alt$all.APM.mut>0, TRUE, FALSE)
        }

if ("Factor.Genomic" %in% names(escape.data)) {
            tmp2.binar = c(escape.data[c("APM.alt", "Epigenetic.feature", "Factor.Genomic")], tmp1.binar)
        }
        else{
            tmp2.binar=c(escape.data[c('APM.alt', 'Epigenetic.feature')], tmp1.binar)
        }
        return(tmp2.binar)

}
