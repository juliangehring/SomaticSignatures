
dfConvertColumns <- function(x, from = "character", to = "factor") {
    idx = sapply(x, is, from)
    x[idx] = lapply(x[idx], as, to)
    return(x)
}


setAs("character", "factor", function(from) {
    return(as.factor(from))
})
