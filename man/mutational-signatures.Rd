\name{mutational-signatures}

\alias{mutational-signatures}
\alias{identifySignatures}
\alias{extractSignatures}

\title{Estimate Somatic Signatures}

\description{

  Estimate somatic signatures from sequence motifs with a selection of
  statistical methods.

}

\usage{
identifySignatures(m, nSigs, decomposition = c("nmf", "pca"), ..., includeFit = FALSE)
}

\arguments{

  \item{m}{GRanges object [required]}

  \item{nSigs}{Number of signatures [integer, required]}

  \item{decomposition}{Method to apply (currently: 'nmf' or 'pca')}

  \item{...}{Additional arguments passed to 'NMF::nmf'}

  \item{includeFit}{Include the fit object returned by the low-level
    decomposition function in the output.}

}

\details{
  
  The 'nmfSignatures', 'pcaSignatures', and 'kmeansSignatures' functions
  estimate a set of 'r' somatic signatures using the NMF, PCA, or
  k-means clustering, respectively.
  
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

  'mutationContext', 'mutationContextMutect'

  'motifMatrix'

  'NMF' package, 'pcaMethods' package, 'prcomp', 'kmeans'

}
