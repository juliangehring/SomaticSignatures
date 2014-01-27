
\name{mutation-distribution}
  
\alias{mutationDistance}
\alias{plotRainfall}

\title{Distributions of mutational locations.}

\description{

  Summary and plotting function for characterizing the distributions of
  mutations along the genome.

}

\usage{
mutationDistance(x)

plotRainfall(x, group, size = 2, alpha = 0.5, space.skip = 0, ...)
}

\arguments{
  
  \item{x}{A 'GRanges' or 'VRanges' object [required].}

  \item{group}{The variable name for color groups [optional].}

  \item{size}{Point size [default: 2]}

  \item{alpha}{Alpha value for points [default: 0.5]}

  \item{space.skip}{Space between chromosomes, as defined by 'plotGrandLinear' [default: 0]}

  \item{...}{Additional arguments passed to 'plotGrandLinear'}
}

\value{
  
  \itemize{
    
    \item{mutationDensity}{The position-sorted GRanges 'x' with the
      additional column 'distance', specifying the distance from the
      previous mutation (or the beginning of the chromosome if it
      happens to be the first mutation on the chromosome.)  }

    \item{plotRainfall}{Object of class 'ggbio', as returned by
      'plotGrandLinear'.}

  }

}
  
\seealso{

  'ggbio::plotGrandLinear'

}

\examples{

library(GenomicRanges)
library(IRanges)

set.seed(1)
chr_len = 100
gr = GRanges(rep(1:3, each = 10),
  IRanges(start = sample.int(chr_len, 30, replace = FALSE), width = 1),
  mutation = sample(c("A", "C", "G", "T"), 30, replace = TRUE))
seqlengths(gr) = rep(chr_len, 3)

p = plotRainfall(gr)
print(p)

}
