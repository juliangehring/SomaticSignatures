identifySignaturesVRanges <- function(vr, group, nSigs, decomposition = c("nmf", "pca"), ..., includeFit = FALSE) {

    decomposition = match.arg(decomposition)

    m = motifMatrix(vr, group, normalize = TRUE)
    
    res = findSignatures(m, nSigs, decomposition, ..., includeFit = includeFit)
    res@decomposition = decomposition
    res@options = list(...)

    return(res)
}


identifySignatures <- function(m, nSigs, decomposition = c("nmf", "pca"), ..., includeFit = FALSE) {

    res = findSignatures(m, nSigs, decomposition, ..., includeFit = includeFit)

    return(res)
}


findSignatures <- function(x, r, decomposition = c("nmf", "pca"), ..., includeFit = FALSE) {

    decomposition = match.arg(decomposition)
    
    dc = switch(decomposition,
        nmf = .nmfSignatures(x, r, ..., includeFit = includeFit),
        pca = .pcaSignatures(x, r, ..., includeFit = includeFit),
        kmeans = .kmeansSignatures(x, r, ..., includeFit = includeFit)
        )

    res = new("MutationalSignatures",
        signatures = dc$w,
        samples = dc$h,
        fitted = dc$v,
        observed = dc$m,
        nSignatures = r)
    if(includeFit)
        res@raw = dc$raw

    return(res)
}


.nmfSignatures <- function(x, r, ..., includeFit = FALSE) {
    
    #args = c(list(...), defaultArgs)
    #args = args[!duplicated(names(args))]
    #y = nmf(x, r, seed = args$seed, ... = unlist(args))

    y = nmf(x, r, ...)

    w = basis(y) ## signatures x k
    h = t(coef(y)) ## samples x k

    ## order signatures
    ord = order(rowMax(t(w)), decreasing = TRUE)
    w = w[ ,ord]
    h = h[ ,ord]
    
    sig_names = paste0("S", 1:r)
    colnames(w) = colnames(h) = sig_names
    v = fitted(y)

    res = list(w = w, h = h, v = v, m = x, r = r)
    if(includeFit)
        res[["raw"]] = y

    return(res)
}


.kmeansSignatures <- function(x, r, ..., includeFit = FALSE) {

    y = kmeans(t(x), centers = r)
    
    w = t(y$centers)
    n_samples = ncol(x)
    h = matrix(0, r, n_samples)
    h[ cbind(as.vector(y$cluster), 1:n_samples) ] = 1
    h = matrix(0, n_samples, r)
    h[ cbind(1:n_samples, as.vector(y$cluster)) ] = 1
    stopifnot(all(rowSums(h) == 1))
    sig_names = paste0("S", 1:r)
    colnames(w) = colnames(h) = sig_names
    v = fitted(y)

    res = list(w = w, h = h, v = v, m = x, r = r)
    if(includeFit)
        res[["raw"]] = y

    return(res)    
}


.pcaSignatures <- function(x, r, ..., includeFit = FALSE) {
  
    y = pca(x, "svd", r, scale = "uv", ...)
    w = scores(y) ## signatures x k
    h = loadings(y) ## samples x k
    v = fitted(y) 
    
    sig_names = paste0("S", 1:r)      
    colnames(w) = colnames(h) = sig_names

    res = list(w = w, h = h, v = v, m = x, r = r)
    if(includeFit)
        res[["raw"]] = y
    
    return(res)    
}
