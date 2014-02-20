
plotVariantAbundance <- function(x, group = NULL, alpha = 0.5, size = 2) {

    df = data.frame(
        vf = altFraction(x),
        cov = totalDepth(x)
        )
    ## only keep the columns of interest, to keep size small
    ## convert to data.frame to have access to GRanges slots, e.g. 'seqnames'
    if(!is.null(group))
        df[ ,group] = as(x, "data.frame")[ ,group]

    df = df[!is.na(df$vf), ]

    p = ggplot(df) + geom_point(aes_string(x = "vf", y = "cov", col = group),
        size = size, alpha = alpha) + xlim(0, 1) +
            xlab("Variant Frequency") + ylab("Coverage") + theme_bw()

    return(p)
}
