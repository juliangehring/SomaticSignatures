assessNumberSignatures <- function(m, nSigs, decomposition = nmfDecomposition, ...,
                                     nReplicates = 1) {

    dev = lapply(nSigs, function(r, m, decomposition, ...) {
        d = lapply(1:nReplicates, function(i) {
            dev = assessOneSignature(m, r, decomposition, ...)
            dev$Replicate = i
            dev
        })
        return(do.call(rbind, d))
    }, m, decomposition, ...)

    gof = do.call(rbind, dev)

    return(gof)
}


plotNumberSignatures <- function(gof) {

    m = melt(gof, id.vars = c("NumberSignatures", "Replicate"),
        measure.vars = c("RSS", "ExplainedVariance"), variable.name = "stat")

    p = ggplot(m, aes_string(x = "NumberSignatures", y = "value", group = "NumberSignatures"))
    p = p + stat_summary(fun.y = mean, fun.ymin = min,
        fun.ymax = max, colour = "red", size = 1.2)
    p = p + geom_point(color = "black", shape = 3)
    p = p + facet_wrap(~stat, nrow = 2, scales = "free")
    p = p + theme_bw() + xlab("Number of Signatures") + ylab("Statistic")

    return(p)
}


## Goodness-of-fit functions ##

assessOneSignature <- function(m, n, decomposition, ...) {

    sigs = identifySignatures(m, n, decomposition, ...)
    gof = data.frame(
        `NumberSignatures` = n,
        `RSS` = rss(fitted(sigs), observed(sigs)),
        `ExplainedVariance` = evar(fitted(sigs), observed(sigs))
        )

    return(gof)
}
