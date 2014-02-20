
mutationDistance <- function(x) {

    x = sort(x) ## must be sorted
    idx_change = start(seqnames(x)) ## where does a new chr begin
    dist = diff(c(1, start(x))) ## to keep the same length
    dist[idx_change] = start(x[idx_change]) ## for new chr: distance to start
    #stopifnot(all(dist > 0)) ## TODO
    #idx_same = (as(dist, "Rle") == 0)
    #dist[idx_same] = dist[idx_same - 1]
    x$distance = dist
    
    return(x)
}


plotRainfall <- function(x, group, size = 2, alpha = 0.5, space.skip = 0, ...) {

    y = mutationDistance(x)

    if(missing(group)) {
        p = plotGrandLinear(y, aes_string(y = "distance"),
            space.skip = space.skip, alpha = alpha, size = size, ...)
    } else {
        p = plotGrandLinear(y, aes_string(y = "distance", col = group),
            space.skip = space.skip, alpha = alpha, size = size, ...)
    }
    p = p + theme_bw() + scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              axis.text = element_text(size = 10)) +
                  xlab("Chromosome") + ylab("Distance [bp]")
    
    return(p)
}
