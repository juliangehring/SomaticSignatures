motifMatrix <- function(vr, group = "sampleNames", normalize = TRUE) {

    alteration = expand.grid(ref = DNA_BASES, alt = DNA_BASES)
    alteration = subset(alteration, ref != alt & ref %in% c("C", "T"))
    alteration = sort(paste0(alteration$ref, alteration$alt))
    motifs = expand.grid(s = DNA_BASES, p = DNA_BASES, a = alteration)
    motifs = sprintf("%s %s.%s", motifs$a, motifs$p, motifs$s)

    df = as(unname(vr), "data.frame")
    if(!(group %in% colnames(df))) {
        stop(sprintf("Column '%s' not present in input object.", group))
    }

    ## form the matrix
    df$motif = factor(paste(df$alteration, df$context), levels = motifs)
    y = sapply(tapply(df$motif, df[ ,group], table), c)

    if(normalize)
        y = t(t(y) / colSums(y))

    return(y)
}
