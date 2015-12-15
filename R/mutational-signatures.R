identifySignaturesVRanges <- function(vr, group, nSigs, decomposition = nmfDecomposition, ...) {

    m = motifMatrix(vr, group, normalize = TRUE)

    res = findSignatures(m, nSigs, decomposition, ...)
    res@decomposition = decomposition
    res@options = list(...)

    return(res)
}


identifySignatures <- function(m, nSigs, decomposition = nmfDecomposition, ...) {

    res = findSignatures(m, nSigs, decomposition, ...)

    return(res)
}


findSignatures <- function(x, r, decomposition = nmfDecomposition, ...) {

    ## check input arguments
    nm = min(dim(x))
    if(!( length(r) == 1 & all.equal(r, as.integer(r)) & r > 0 & r <= nm)) {
        msg = "'nSigs | r' must be a single, positive integer in the range [1,%d]"
        stop(sprintf(msg, nm))
    }

    if(!is.function(decomposition))
        stop("'decomposition' must be a function.")

    dc = decomposition(x, r, ...)

    ## check returned object
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
