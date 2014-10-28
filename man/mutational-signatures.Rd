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
identifySignatures(m, nSigs, decomposition = nmfSignatures, ...)
}

\arguments{

  \item{m}{Motif matrix, as returned by 'motifMatrix' [required].}

  \item{nSigs}{Number of signatures [integer, required].}

  \item{decomposition}{Function to apply for the matrix
    decompositon. The methods NMF and PCA are already implemented in the
    functions 'nmfSignatures' and 'pcaSignatures', respectively.}

  \item{...}{Additional arguments passed to the 'decomposition'
    function.}

}

\details{
  
  'identifySignatures' estimate a set of 'r' somatic signatures, based
  on a matrix decompositon method (such as NMF, PCA).
  
}

\value{

  An object of class 'MutationalSignatures'.

}


\seealso{

  The predefined decomposition functions: 'nmfSignatures' and 'pcaSignatures'

  'mutationContext', 'mutationContextMutect'

  'motifMatrix'

  class 'MutationalSignatures'
  
}


\examples{
data("sca_mm", package = "SomaticSignatures")

sigs = identifySignatures(sca_mm, 5)
}
