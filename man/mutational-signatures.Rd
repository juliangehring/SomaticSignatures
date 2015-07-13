\name{mutational-signatures}

\alias{mutational-signatures}
\alias{identifySignatures}
\alias{findSignatures}

\title{Estimate Somatic Signatures}

\description{

  Estimate somatic signatures from sequence motifs with a selection of
  statistical methods.

}

\usage{
identifySignatures(m, nSigs, decomposition = nmfDecomposition, ...)
}

\arguments{

  \item{m}{Motif matrix, as returned by 'motifMatrix' [required].}

  \item{nSigs}{Number of signatures [integer, required].}

  \item{decomposition}{Function to apply for the matrix
    decomposition. The methods NMF and PCA are already implemented in the
    functions 'nmfDecomposition' and 'pcaDecomposition', respectively.}

  \item{...}{Additional arguments passed to the 'decomposition'
    function.}

}

\details{
  
  'identifySignatures' estimate a set of 'r' somatic signatures, based
  on a matrix decomposition method (such as NMF, PCA).
  
}

\value{

  An object of class 'MutationalSignatures'.

}


\seealso{

  The predefined decomposition functions: \code{\link{nmfDecomposition}}
  and \code{\link{pcaDecomposition}}

  \code{\link{mutationContext}}, \code{\link{mutationContextMutect}}

  \code{\link{motifMatrix}}

  \code{\link{MutationalSignatures}} class
  
}


\examples{
data("sca_mm", package = "SomaticSignatures")

sigs = identifySignatures(sca_mm, 5)
}
