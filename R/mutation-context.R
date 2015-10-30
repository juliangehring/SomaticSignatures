mutationContext <- function(vr, ref, k = 3, strand = FALSE, unify = TRUE, check = FALSE) {

    ## only SNV substitutions beyond this point
    if(!all(ref(vr) %in% DNA_BASES & alt(vr) %in% DNA_BASES))
        stop("Only SNV substitutions are currently supported.")
    if(k %% 2 != 1)
        stop("'k' must be odd.")
    mid = (k + 1)/2

    gr = granges(vr) ## drop mcols
    ranges = resize(gr, k, fix = "center")
    context = getSeq(ref, ranges)

    ref_base = DNAStringSet(ref(vr))
    alt_base = DNAStringSet(alt(vr))
    
    ## check against ref
    if(check) {
        ref0 = subseq(context, mid, mid) ## reference base from 'ref', for checking
        idx_invalid = ( ref0 != ref_base )
        if(any(idx_invalid))
            warning(sprintf("References do not match in %d cases",
                            sum(idx_invalid)))
    }
    
    ## convert to the plus strand, generally not needed
    if(strand) {
        s = strand(gr)
        if(any(s == "*"))
            stop("The strand must be explicit, in order to read the correct strand.")
        idx_minus = (s == "-")
        context[idx_minus] = reverseComplement(context[idx_minus])
        s[idx_minus] = "+"
        strand(gr) = s
    }

    ## convert to alterations starting with "C/T"
    if(unify) {
        idx_complement = ref_base %in% c("A", "G")
        context[idx_complement] = reverseComplement(context[idx_complement])
        ref_base[idx_complement] = reverseComplement(ref_base[idx_complement])
        alt_base[idx_complement] = reverseComplement(alt_base[idx_complement])
    }
    
    subseq(context, mid, mid) = "."
    alteration = xscat(ref_base, alt_base)

    vr$alteration = alteration
    vr$context = context

    return(vr)
}


mutationContextMutect <- function(vr, k = 3, unify = TRUE) {

    ## context: {r|k|r}, 2*r+k = w
    w = width(vr[1]$context)
    if(w %% 2 != 1)
        stop("width of context must be odd.")
    r = (w - k)/2 ## must be odd
    context = subseq(vr$context, r + 1, w - r)

    ref_base = DNAStringSet(ref(vr))
    alt_base = DNAStringSet(alt(vr))
    
    ## convert to alterations starting with "C/T"
    if(unify) {
        idx_complement = ref_base %in% c("A", "G")
        context[idx_complement] = reverseComplement(context[idx_complement])
        ref_base[idx_complement] = reverseComplement(ref_base[idx_complement])
        alt_base[idx_complement] = reverseComplement(alt_base[idx_complement])
    }
    
    alteration = xscat(ref_base, alt_base)
    
    df = DataFrame(alteration = alteration, context = context)      
    mcols(vr) = cbind(mcols(vr), df)

    return(vr)
}
