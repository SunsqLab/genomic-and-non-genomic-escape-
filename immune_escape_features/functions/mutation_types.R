# Return Level5 mutation labels for the requested mutation class.
MutationType <- function(type.mut) {
    type.all <- c(
        "In_Frame_Ins", "In_Frame_Del", "Nonstop_Mutation", "Translation_Start_Site", "Frame_Shift_Ins", "Frame_Shift_Del",
        "Splice_Site", "Nonsense_Mutation", "Missense_Mutation", "Startloss_Mutation", "Silent", "Intron", "5'UTR", "3'UTR", "IGR", "5'Flank", "3'Flank", "3'Flank", "RNA"
    )
    type.truncating <- c("Frame_Shift_Del", "Nonsense_Mutation", "Frame_Shift_Ins", "Splice_Site")

    type.coding <- c(
        "Missense_Mutation", "Nonstop_Mutation", "Nonsense_Mutation", "Splice_Site", "Startloss_Mutation", "Silent",
        "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del"
    )
    type.synonymous <- "Silent"

    type.nonsynonymous <- c("Missense_Mutation", "Nonstop_Mutation", "Nonsense_Mutation", "Startloss_Mutation", "Splice_Site")
    type.nonsilent <- c(
        "In_Frame_Ins", "In_Frame_Del", "Nonstop_Mutation", "Translation_Start_Site", "Frame_Shift_Ins", "Frame_Shift_Del",
        "Splice_Site", "Nonsense_Mutation", "Missense_Mutation", "Startloss_Mutation"
    )

    type.missense <- "Missense_Mutation"
    type.nonsense <- "Nonsense_Mutation"
    type.inframe <- c("In_Frame_Ins", "In_Frame_Del")
    type.frameshift <- c("Frame_Shift_Ins", "Frame_Shift_Del")

type.exon <- c(
        "Missense_Mutation", "Nonstop_Mutation", "Nonsense_Mutation", "Startloss_Mutation", "Silent",
        "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del"
    )
    type.nonexon <- c("Intron", "5'UTR", "3'UTR", "IGR", "5'Flank", "3'Flank", "3'Flank", "RNA")

    results <- switch(type.mut,
        truncating = type.truncating,
        coding = type.coding,
        synonymous = type.synonymous,
        nonsynonymous = type.nonsynonymous,
        nonsilent = type.nonsilent,
        missense = type.missense,
        nonsense = type.nonsense,
        inframe = type.inframe,
        frameshift = type.frameshift,
        exon = type.exon,
        nonexon = type.nonexon,
        all = type.all
    )

    return(results)
}
