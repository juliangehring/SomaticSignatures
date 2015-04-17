setClass("MutationalSignatures",
         slots = c(
             signatures = "matrix",
             samples = "matrix",
             observed = "matrix",
             fitted = "matrix",
             nSignatures = "numeric",
             raw = "ANY",
             decomposition = "character",
             options = "list"
             )
         )

setMethod("show",
          signature(object = "MutationalSignatures"),
          function(object) {
              cat("MutationalSignatures:\n")
              samples = rownames(samples(object))
              signatures = colnames(samples(object))
              motifs = rownames(signatures(object))
              cat(showSome(samples, "Samples", "  "))
              cat(showSome(signatures, "Signatures", "  "))
              cat(showSome(motifs, "Motifs", "  "))
              cat(sprintf("  Decomposition: %s\n", object@decomposition))
          })


setGeneric("signatures",
           function(object)
           standardGeneric("signatures")
           )

setMethod("signatures",
          signature(object = "MutationalSignatures"),
          function(object) {
              return(object@signatures)
          })


setGeneric("samples",
           function(object)
           standardGeneric("samples")
           )

setMethod("samples",
          signature(object = "MutationalSignatures"),
          function(object) {
              return(object@samples)
          })


setGeneric("observed",
           function(object)
           standardGeneric("observed")
           )

setMethod("observed",
          signature(object = "MutationalSignatures"),
          function(object) {
              return(object@observed)
          })


setGeneric("fitted",
           function(object)
           standardGeneric("fitted")
           )

setMethod("fitted",
          signature(object = "MutationalSignatures"),
          function(object) {
              return(object@fitted)
          })
