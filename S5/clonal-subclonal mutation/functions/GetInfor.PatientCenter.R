# Original SampleID/PatientID lookup; nested PatientCenter format retained.
GetInfor.PatientCenter <- function(patient.center, PatientID=NULL, SampleID=NULL,
        colNames=NULL, ...
){

    if(is.null(PatientID) & is.null(SampleID)){
        SampleID = patient.center$SampleInfo$SampleID

    }
    if(is.null(PatientID)){
        PatientID = patient.center$SampleInfo$PatientID[match(SampleID, patient.center$SampleInfo$SampleID)]
        t.results = data.frame(PatientID=PatientID, SampleID=SampleID)
    } else if(is.null(SampleID)){
        SampleID = patient.center$SampleInfo$SampleID[match(PatientID, patient.center$SampleInfo$PatientID)]
        t.results = data.frame(PatientID=PatientID, SampleID=SampleID)
    }

    if(is.null(colNames)){ return(t.results) }

    all.col = unique(unlist(sapply(patient.center, colnames)))
    if(!all(colNames %in% all.col)){
        tmp = setdiff(colNames, all.col)
        stop("Information not found in PatientCenter: ", tmp, "...")
    }
    colNames = colNames[colNames %in% all.col]
    colNames = setdiff(colNames, c("SampleID", "PatientID"))

    if(length(colNames) == 0){

        return(t.results)
    }

    t.colNames = colNames

    temp = lapply(patient.center, function(x){
                if(any(t.colNames %in% colnames(x))){
                    temp = intersect(t.colNames, colnames(x))

                    t.colNames <<- setdiff(t.colNames, temp)

                    if("SampleID" %in% colnames(x)){
                        result = x[match(SampleID, x$SampleID), temp, drop=FALSE]

                        return(result)
                    }else{
                        result = x[match(PatientID, x$PatientID), temp, drop=FALSE]

                        return(result)
                    }
                }
            })
    names(temp) = NULL
    temp = temp[!sapply(temp, is.null)]
    t.results = cbind(t.results, do.call(cbind, temp))

    library(tidyverse)
    t.results <- t.results %>%
            dplyr::filter(...)

    return(t.results)
}
