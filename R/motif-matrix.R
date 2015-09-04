motifMatrix <- function(vr, group = "sampleNames", normalize = TRUE) {

    motifs = c("CA A.A", "CA A.C", "CA A.G", "CA A.T", "CA C.A", "CA C.C",
               "CA C.G", "CA C.T", "CA G.A", "CA G.C", "CA G.G", "CA G.T", "CA T.A",
               "CA T.C", "CA T.G", "CA T.T", "CG A.A", "CG A.C", "CG A.G", "CG A.T",
               "CG C.A", "CG C.C", "CG C.G", "CG C.T", "CG G.A", "CG G.C", "CG G.G",
               "CG G.T", "CG T.A", "CG T.C", "CG T.G", "CG T.T", "CT A.A", "CT A.C",
               "CT A.G", "CT A.T", "CT C.A", "CT C.C", "CT C.G", "CT C.T", "CT G.A",
               "CT G.C", "CT G.G", "CT G.T", "CT T.A", "CT T.C", "CT T.G", "CT T.T",
               "TA A.A", "TA A.C", "TA A.G", "TA A.T", "TA C.A", "TA C.C", "TA C.G",
               "TA C.T", "TA G.A", "TA G.C", "TA G.G", "TA G.T", "TA T.A", "TA T.C",
               "TA T.G", "TA T.T", "TC A.A", "TC A.C", "TC A.G", "TC A.T", "TC C.A",
               "TC C.C", "TC C.G", "TC C.T", "TC G.A", "TC G.C", "TC G.G", "TC G.T",
               "TC T.A", "TC T.C", "TC T.G", "TC T.T", "TG A.A", "TG A.C", "TG A.G",
               "TG A.T", "TG C.A", "TG C.C", "TG C.G", "TG C.T", "TG G.A", "TG G.C",
               "TG G.G", "TG G.T", "TG T.A", "TG T.C", "TG T.G", "TG T.T")

    df = as(unname(vr), "data.frame")
    if(!(group %in% colnames(df))) {
        stop(sprintf("Column '%s' not present in input object.", group))
    }

    ## form the matrix
    group_string = paste0(group, " ~ motif")
    df$motif = factor(constructMotif(df$alteration, df$context))
    y = sapply(tapply(df$motif, df[ ,group], table), c)

    if(normalize)
        y = t(t(y) / colSums(y))

    return(y)
}


constructMotif <- function(alteration, context) {

    motif = paste(alteration, context)

    return(motif)
}
