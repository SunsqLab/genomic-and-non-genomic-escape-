# Calculate overlap coordinates and coverage percentages for GRanges or IRanges objects.
olRanges <- function(query, subject, output = "gr", ..., parallel = FALSE) {
    require(GenomicRanges)
    require(IRanges)

    if (!((class(query) == "GRanges" & class(subject) == "GRanges") | (class(query) == "IRanges" & class(subject) == "IRanges"))) {
        stop("Query and subject need to be of same class, either GRanges or IRanges!")
    }

    if (class(query) == "GRanges") {
        seqlengths(query) <- rep(NA, length(seqlengths(query)))
        seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
    }

    olindex <- as.matrix(findOverlaps(query, subject, ...))

    query <- query[olindex[, 1]]
    subject <- subject[olindex[, 2]]
    olma <- cbind(Qstart = start(query), Qend = end(query), Sstart = start(subject), Send = end(subject))

    startup <- olma[, "Sstart"] < olma[, "Qstart"]
    enddown <- olma[, "Send"] > olma[, "Qend"]
    startin <- olma[, "Sstart"] >= olma[, "Qstart"] & olma[, "Sstart"] <= olma[, "Qend"]
    endin <- olma[, "Send"] >= olma[, "Qstart"] & olma[, "Send"] <= olma[, "Qend"]

    olup <- startup & endin
    oldown <- startin & enddown
    inside <- startin & endin
    contained <- startup & enddown

    OLtype <- rep("", length(olma[, "Qstart"]))
    OLtype[olup] <- "oldown"
    OLtype[oldown] <- "olup"
    OLtype[inside] <- "contained"
    OLtype[contained] <- "inside"

    OLstart <- rep(0, length(olma[, "Qstart"]))
    OLend <- rep(0, length(olma[, "Qstart"]))
    OLstart[olup] <- olma[, "Qstart"][olup]
    OLend[olup] <- olma[, "Send"][olup]
    OLstart[oldown] <- olma[, "Sstart"][oldown]
    OLend[oldown] <- olma[, "Qend"][oldown]
    OLstart[inside] <- olma[, "Sstart"][inside]
    OLend[inside] <- olma[, "Send"][inside]
    OLstart[contained] <- olma[, "Qstart"][contained]
    OLend[contained] <- olma[, "Qend"][contained]

    OLlength <- (OLend - OLstart) + 1
    OLpercQ <- OLlength / width(query) * 100
    OLpercS <- OLlength / width(subject) * 100

    oldf <- data.frame(Qindex = olindex[, 1], Sindex = olindex[, 2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
    if (class(query) == "GRanges") {
        oldf <- cbind(space = as.character(seqnames(query)), oldf)
    }
    if (output == "df") {
        return(oldf)
    }
    if (output == "gr") {
        if (class(query) == "GRanges") {
            elementMetadata(query) <- cbind(cbind(as.data.frame(elementMetadata(query)), oldf), as.data.frame(elementMetadata(subject)))
        }
        if (class(query) == "IRanges") {
            query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start = oldf[, "Qstart"], end = oldf[, "Qend"]), strand = Rle(strand(rep("+", length(query)))), oldf)
        }
        return(query)
    }
}
