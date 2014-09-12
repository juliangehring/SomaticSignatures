\name{cluster-spectrum}

\alias{cluster-spectrum}
\alias{clusterSpectrum}

\title{Cluster Mutational Spectrum}

\description{

  Cluster the mutational spectrum by sample or motif.

}

\usage{
clusterSpectrum(m, by = c("sample", "motif"), distance = "Cosine", ...)
}

\arguments{

  \item{m}{Mutational spectrum matrix}

  \item{by}{Dimension to cluster by.}

  \item{distance}{Distance function used in the clustering.}

  \item{...}{Additional arguments passed to 'hclust'.}

}

\details{
  
  Hierachical clustering of the motif matrix aka mutational spectrum.
  
}

\value{

  An 'hclust' object.
  
}


\seealso{

  hclust

  'proxy' package

}
