identifySignaturesVRanges <- function(vr, group, nSigs, decomposition = nmfSignatures, ...) {

    m = motifMatrix(vr, group, normalize = TRUE)
    
    res = findSignatures(m, nSigs, decomposition, ...)
    res@decomposition = decomposition
    res@options = list(...)

    return(res)
}


identifySignatures <- function(m, nSigs, decomposition = nmfSignatures, ...) {

    res = findSignatures(m, nSigs, decomposition, ...)

    return(res)
}


findSignatures <- function(x, r, decomposition = nmfSignatures, ...) {

    if(!is.function(decomposition))
        stop("'decomposition' must be a function.")
    
    dc = decomposition(x, r, ...)

    required_names = c("m", "w", "h", "v")
    if(any(!(required_names %in% names(dc)))) {
        msg = paste0("The decomposition function must return a list with names: ",
            paste(required_names, ", "))
        stop(msg)
    }
    
    res = new("MutationalSignatures",
        signatures = dc$w,
        samples = dc$h,
        fitted = dc$v,
        observed = dc$m,
        nSignatures = r)
#    if(includeFit)
#        res@raw = dc$raw

    return(res)
}
