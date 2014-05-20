\name{mutational-signatures}

\alias{mutational-signatures}
\alias{motifMatrix}
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
motifMatrix(x, group = "sample", normalize = TRUE)

findSignatures(x, r, method = c("nmf", "pca", "kmeans"), ...)
nmfSignatures(x, r, seed = "ica", ...)
kmeansSignatures(x, r, ...)
pcaSignatures(x, r, ...)

plotSamplesObserved(s, group = "sample")
plotSignatureMap(s)
plotSignatures(s)
plotSampleMap(s)
plotSamples(s)
}

\arguments{

  \item{x}{GRanges object [required]}

  \item{group}{Grouping variable name [character, default: 'sample']}

  \item{normalize}{Normalize to frequency}

  \item{r}{Number of signatures [integer, required]}

  \item{method}{Method to apply (currently: 'nmf' or 'kmeans')}

  \item{seed}{seed for NMF, default: "ica"}

  \item{...}{Additional arguments passed to 'NMF::nmf'}

  \item{s}{results signature object}

}

\details{

  The 'motifMatrix' function transforms the metadata columns
  of a 'VRanges' object, as returned by the 'mutationContext' function,
  to a matrix of the form 'motifs x samples'.  This constitutes the
  bases for the estimation of the signatures.
  
  The 'nmfSignatures', 'pcaSignatures', and 'kmeansSignatures' functions
  estimate a set of 'r' somatic signatures using the NMF, PCA, or
  k-means clustering, respectively.

  With the plotting function, the obtained signatures and their
  occurrance in the samples can be visualized either as a heatmap
  ('plotSignatureMap', 'plotSampleMap') or a barchart ('plotSignature',
  'plotSamples').
  
}

\seealso{

  'mutationContext', 'mutationContextMutect'

  'NMF' package, 'pcaMethods' package, 'prcomp', 'kmeans'

}
