kmerFrequency <- function(ref, n = 1e4, k = 1, ranges = as(seqinfo(ref), "GRanges")) {

    w = width(ranges) - k + 1 ## width 'k' is still good
    if(any(idx_short <- (w <= 0))) {
        ranges = ranges[!idx_short]
        w = w[!idx_short]
        message(sprintf("%d ranges dropped since shorter than %d.",
                        sum(idx_short), k))
    }   
    
    ## tabulizing is faster for fewer seqnames, but the gain is neglectible
    ## compared to the time we spend getting the sequence
    a_sample = sample.int(length(ranges), n, replace = TRUE, prob = w)
    a_start = sapply(w[a_sample], sample.int, 1) + start(ranges[a_sample])

    ## check
    gr = GRanges(seqnames(ranges[a_sample]),
        IRanges(start = a_start, width = k))

    seq = getSeq(ref, seqnames(ranges[a_sample]), start = a_start, width = k)
    freq = table(seq)
    freq = freq / sum(freq)

    return(freq)
}
