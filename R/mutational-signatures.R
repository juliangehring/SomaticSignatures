
mutationContextMatrix <- function(x, group = "sample", normalize = TRUE) {

    d = as(mcols(x), "data.frame")
    group_string = paste0(group, " ~ motif")
    d$motif = factor(paste(d$alteration, d$context))
    y = t(acast(d, group_string, value.var = "motif", fun.aggregate = length))

    if(normalize)
        y = t(t(y) / colSums(y))
    
    return(y)
}


findSignatures <- function(x, r, method = c("nmf", "pca", "kmeans"), ...) {

    method = match.arg(method)
    
    y = switch(method,
        nmf = nmfSignatures(x, r, ...),
        pca = pcaSignatures(x, r, ...),
        kmeans = kmeansSignatures(x, r, ...)
        )
    
    return(y)
}


nmfSignatures <- function(x, r, seed = "ica", ...) {
    
    y = nmf(x, r, seed = seed, ...)

    ## extract the data
    w = basis(y) ## signatures x k
    h = t(coef(y)) ## samples x k
    sig_names = paste0("S", 1:r)
    colnames(w) = colnames(h) = sig_names
    v = fitted(y)
    res = list(w = w, h = h, v = v, raw = y)

    return(res)
}


kmeansSignatures <- function(x, r, ...) {

    y = kmeans(t(x), centers = r)
    
    ## extract the data
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
    res = list(w = w, h = h, v = v, raw = y)

    return(res)    
}


pcaSignatures <- function(x, r, ...) {
  
    #pca = prcomp(x, scale = scale) ##
    #w = pca$rotation ## signatures x k
    #h = pca$x ## samples x k
    #v = scale(h %*% t(w), pca$center, pca$scale)
    
    y = pca(x, "svd", r, scale = "uv", ...)
    w = scores(y) ## signatures x k
    h = loadings(y) ## samples x k
    v = fitted(y) 
    
    sig_names = paste0("S", 1:r)
    colnames(w) = colnames(h) = sig_names
    res = list(w = w, h = h, v = v, raw = y)
  
    return(res)    
}
