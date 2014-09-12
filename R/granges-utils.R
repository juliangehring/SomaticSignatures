ucsc <- function(x) {
    suppressMessages(seqlevelsStyle(x) <- "UCSC") ## '<-' needed
    genome(x) = NA
    return(x)
}


ncbi <- function(x) {
    suppressMessages(seqlevelsStyle(x) <- "NCBI") ## '<-' needed
    genome(x) = NA
    return(x)
}


seqchar <- function(x) {
    y = as(seqnames(x), "character")
    return(y)
}
