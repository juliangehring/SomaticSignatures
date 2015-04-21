## The following functions for correlationg signatures is currently not exported nor used in the package ##

retreiveSignatures <- function(vr, group, signatures, method = c("cospos", "cosine", "spearman", "pearson")) {

    method = match.arg(method)
    m = motifMatrix(vr, group, normalize = TRUE)

    cc = correlateSignatures(m, signatures, method)

    res = new("MutationalSignatures",
        signatures = signatures,
        samples = cc,
        fitted = m,
        observed = m,
        nSignatures = ncol(signatures),
        group = group)

    return(res)
}


correlateSignatures <- function(m, s, method = c("cospos", "cosine", "spearman", "pearson")) {

    method = match.arg(method)
    stopifnot(nrow(m) == nrow(s))

    cc = switch(method, ## sample x signature
        cospos = coscor(m, s, cospos),
        cosine = coscor(m, s, cossim),
        spearman = cor(m, s, method = "spearman"),
        pearson = cor(m, s, method = "pearson"))

    rownames(cc) = colnames(m)
    colnames(cc) = colnames(s)

    return(cc)
}

 
cossim <- function(x, y) {
    res = crossprod(x, y)/sqrt(crossprod(x) * crossprod(y))
    return(res)
}

cospos <- function(x, y) {
    1 - 2 * acos(cossim(x, y)) / pi
}

coscor <- function(a, b, fun = cossim) {

    na = ncol(a)
    nb = ncol(b)
    res = matrix(0, na, nb)
    for(i in 1:na) {
        for(j in 1:nb) {
            res[i,j] = fun(as.vector(a[ ,i]), as.vector(b[ ,j]))
        }
    }

    return(res)
}
