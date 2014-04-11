readMutect <- function(file, columns, strip = FALSE) {

    x = read.delim(file,
        comment.char = "#", stringsAsFactor = TRUE, quote = "")
    
    ## no 'with' statement here, as it yields a NOTE in the package building
    vr = VRanges(x$contig, IRanges(x$position, x$position),
                ref = x$ref_allele,
                alt = x$alt_allele,
                totalDepth = x$t_ref_count + x$t_alt_count,
                altDepth = x$t_alt_count,
                refDepth = x$t_ref_count,
                sampleNames = x$tumor_name,
                softFilterMatrix = matrix(x$judgement %in% "KEEP")
        )
        
    if(strip)
        return(vr)

    if(missing(columns))
        columns = setdiff(names(x), c("contig", "position", "ref_allele",
            "alt_allele", "t_ref_count", "t_alt_count", "tumor_name"))

    mcols(vr) = subset(x, select = columns)
    ## convert the context to DNAStringSet
    bs = BStringSet(x$context)
    subseq(bs, 4, 4) = "."
    vr$context = DNAStringSet(bs)

    return(vr)
}
