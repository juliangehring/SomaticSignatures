dfConvertColumns <- function(x, from = "character", to = "factor") {
    idx = sapply(x, is, from)
    x[idx] = lapply(x[idx], as, to)
    return(x)
}


setAs("character", "factor", function(from) {
    return(as.factor(from))
})


showSome <- function(x, name, indent="") {
    res <- sprintf("%s%s (%d): %s\n",
                   indent,
                   name,
                   length(x),
                   paste(selectSome(x), collapse=", ")
                   )
    return(res)
}
