grangesPlain <- function(x) {
    mcols(x) = NULL
    x = as(x, "GRanges")
    return(x)
}


ucsc <- function(x) {
    suppressMessages(seqnameStyle(x) <- "UCSC") ## '<-' needed
    genome(x) = NA
    return(x)
}


ncbi <- function(x) {
    suppressMessages(seqnameStyle(x) <- "NCBI") ## '<-' needed
    genome(x) = NA
    return(x)
}


seqchar <- function(x) {
    y = as(seqnames(x), "character")
    return(y)
}
