motifMatrix <- function(vr, group = "sampleNames", normalize = TRUE) {

    voi <- if( group %in% names(mcols(vr)) ) {
        mcols(vr)[ ,group]
    } else {
        df = as(unname(vr), "data.frame")
        if( !(group %in% colnames(df)) ) {
            stop(sprintf("Column '%s' not present in input object.", group))
        }
        df[ ,group]
    }

    ## form the matrix
    motif = factor(paste(vr$alteration, vr$context),
                   levels = constructMotifs3())
    y = as(table(motif, voi), "matrix")
    dimnames(y) = unname(dimnames(y))

    if(normalize) {
        y = t(t(y) / colSums(y))
    }

    return(y)
}


constructMotifs3 <- function() {

    alteration = expand.grid(ref = DNA_BASES, alt = DNA_BASES)
    alteration = subset(alteration, ref != alt & ref %in% c("C", "T"))
    alteration = sort(paste0(alteration$ref, alteration$alt))
    motifs = expand.grid(s = DNA_BASES, p = DNA_BASES, a = alteration)
    motifs = sprintf("%s %s.%s", motifs$a, motifs$p, motifs$s)

    return(motifs)
}
