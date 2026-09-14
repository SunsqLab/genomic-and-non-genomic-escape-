# Original cohort RDS loader. Set DataCenter.file and functions.dir in the main script.
CombineData.XTeam <- function(ID.xteams,
        file.query
){

    if(length(ID.xteams) == 1){

          source(file.path(functions.dir, "DatasetLabels.R"))
        ID.xteams=get(ID.xteams)
    }
    all.data.xteam = readRDS(DataCenter.file)
    results = lapply(ID.xteams, function(i){
                cat("Combining dataset: ", i, "\n")
                index.xteam = match(i, all.data.xteam$ID); if(is.na(i)){ stop("Dataset ID was not found in DataCenter.")}
                dataset.path = all.data.xteam$Path[index.xteam]
                t.file = file.path(dataset.path, file.query)
                if(!file.exists(t.file)){
                    return(NULL)
                }else{
                    temp = readRDS(t.file)
                }
                return(temp)
            })

    names(results) = ID.xteams

    return(results)
}
