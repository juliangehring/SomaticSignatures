\name{mutational-signatures}

\alias{mutational-signatures}
\alias{findSignatures}
\alias{nmfSignatures}
\alias{kmeansSignatures}
\alias{pcaSignatures}
\alias{extractSignatures}
\alias{plotSamplesObserved}
\alias{plotSignatureMap}
\alias{plotSignatures}
\alias{plotSampleMap}
\alias{plotSamples}

\title{Estimate Somatic Signatures}

\description{

  Estimate somatic signatures from sequence motifs with a selection of
  statistical methods.

}

\usage{
findSignatures(x, r, method = c("nmf", "pca", "kmeans"), includeRaw = FALSE, ...)
nmfSignatures(x, r, seed = "ica", includeRaw = FALSE, ...)
kmeansSignatures(x, r, includeRaw = FALSE, ...)
pcaSignatures(x, r, includeRaw = FALSE, ...)

plotSamplesObserved(s, group = "sample")
plotSignatureMap(s)
plotSignatures(s)
plotSampleMap(s)
plotSamples(s)
}

\arguments{

  \item{x}{GRanges object [required]}

  \item{group}{Grouping variable name [character, default: 'sample']}

  \item{r}{Number of signatures [integer, required]}

  \item{method}{Method to apply (currently: 'nmf' or 'kmeans')}

  \item{seed}{seed for NMF, default: "ica"}

  \item{includeRaw}{Include the raw object returned by the low-level
    decomposition function in the output.}

  \item{...}{Additional arguments passed to 'NMF::nmf'}

  \item{s}{results signature object}

}

\details{
  
  The 'nmfSignatures', 'pcaSignatures', and 'kmeansSignatures' functions
  estimate a set of 'r' somatic signatures using the NMF, PCA, or
  k-means clustering, respectively.

  With the plotting function, the obtained signatures and their
  occurrance in the samples can be visualized either as a heatmap
  ('plotSignatureMap', 'plotSampleMap') or a barchart ('plotSignature',
  'plotSamples').
  
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
    
    \item{raw}{Raw object returned by the low-level decomposition
      function, if 'includeRaw' is true.}

  }
  
}


\seealso{

  'mutationContext', 'mutationContextMutect'

  'NMF' package, 'pcaMethods' package, 'prcomp', 'kmeans'

}
