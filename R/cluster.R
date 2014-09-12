clusterSpectrum <- function(m, by = c("sample", "motif"), distance = "Cosine", ...) {

    by = match.arg(by)
    
    if(by == "motif")
        m = t(m)

    d = dist(m, method = distance)
    h = hclust(d)

    return(h)
}
