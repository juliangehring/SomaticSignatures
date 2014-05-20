motifMatrix <- function(x, group = "sample", normalize = TRUE) {

    df = as(unname(x), "data.frame")
    if(!(group %in% colnames(df))) {
        stop(sprintf("Column '%s' not present in input object.", group))
    }

    group_string = paste0(group, " ~ motif")
    df$motif = factor(paste(df$alteration, df$context))
    y = t(acast(df, group_string, value.var = "motif", fun.aggregate = length))

    if(normalize)
        y = t(t(y) / colSums(y))
    
    return(y)
}
