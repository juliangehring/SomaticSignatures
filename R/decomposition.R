nmfSignatures <- function(x, r, ..., includeFit = FALSE) {
    
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


kmeansSignatures <- function(x, r, ..., includeFit = FALSE) {

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


pcaSignatures <- function(x, r, ..., includeFit = FALSE) {
  
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
