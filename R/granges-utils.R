ucsc <- function(x) {
    suppressMessages(seqlevelsStyle(x) <- "UCSC") ## '<-' needed
    genome(x) = NA ## avoid mismatches in 'genome' slots for overlaps
    return(x)
}


ncbi <- function(x) {
    suppressMessages(seqlevelsStyle(x) <- "NCBI") ## '<-' needed
    genome(x) = NA ## avoid mismatches in 'genome' slots for overlaps
    return(x)
}


seqchar <- function(x) {
    y = as(seqnames(x), "character")
    return(y)
}
