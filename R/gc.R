
gcContent <- function(regions, ref) {
    
    seq = getSeq(ref, regions)
    gc = letterFrequency(seq, "GC")
    acgt = letterFrequency(seq, "ACGT")
    r = as.vector(gc / acgt)
    r[is.nan(r)] = r

    return(r)
}
