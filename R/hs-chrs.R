hsToplevel <- function() {
    chrs = c(1:22, "X", "Y", "MT")
    return(chrs)
}

hsAutosomes <- function() {
    chrs = as.character(1:22)
    return(chrs)
}

hsAllosomes <- function() {
    chrs = c("X", "Y")
    return(chrs)
}

hsLinear <- function() {
    chrs = c(1:22, "X", "Y")
    return(chrs)
}
