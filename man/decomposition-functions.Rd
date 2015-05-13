\name{decomposition-signatures}

\alias{decomposition-signatures}
\alias{nmfDecomposition}
\alias{pcaDecomposition}
\alias{nmfSignatures}
\alias{pcaSignatures}

\title{Decomposition Functions for Somatic Signatures}

\description{

  Estimate somatic signatures from sequence motifs with a selection of
  statistical methods.

}

\usage{
  nmfDecomposition(x, r, ..., includeFit = FALSE)
  pcaDecomposition(x, r, ..., includeFit = FALSE)
}

\arguments{

  \item{x}{GRanges object [required]}

  \item{r}{Number of signatures [integer, required]}

  \item{...}{Additional arguments passed to 'NMF::nmf' or
    'pcaMethods::pca'.}

  \item{includeFit}{Include the fit object returned by the low-level
    decomposition function in the output.}

}

\details{
  
  The 'nmfDecomposition' and 'pcaDecomposition' functions estimate a set
  of 'r' somatic signatures using the NMF or PCA, respectively.

  In previous versions of the package, these functions were known as
  'nmfSignatures' and 'pcaSignatures', respectively.  While they are
  still available, we recommend using the new naming convention.
  
}

\value{

  The 'signature' functions return a list with the elements:

  \itemize{
    
    \item{w}{Matrix of the form 'motif x signature'}
    
    \item{h}{Matrix of the form 'sample x signature'}
    
    \item{v}{Matrix of the form 'motif x sample', containing the
      reconstruction of 'm' from 'w' and 'h'.}
    
    \item{m}{Input matrix 'm'}
    
    \item{r}{Number of signatures.}
    
    \item{fit}{Fit object returned by the low-level decomposition
      function, if 'includeFit' is true.}

  }
  
}


\seealso{
  
  \code{\link[NMF]{NMF}} package

  \code{\link[pcaMethods]{pcaMethods}} package
 
  \code{\link[stats]{prcomp}}

}
