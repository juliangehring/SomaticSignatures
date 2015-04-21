motifMatrix <- function(vr, group = "sampleNames", normalize = TRUE) {

    df = as(unname(vr), "data.frame")
    if(!(group %in% colnames(df))) {
        stop(sprintf("Column '%s' not present in input object.", group))
    }

    ## form the matrix
    group_string = paste0(group, " ~ motif")
    df$motif = factor(constructMotif(df$alteration, df$context))
    y = t(acast(df, group_string, value.var = "motif", fun.aggregate = length))

    if(normalize)
        y = t(t(y) / colSums(y))
    
    return(y)
}


constructMotif <- function(alteration, context) {

    motif = paste(alteration, context)

    return(motif)
}
