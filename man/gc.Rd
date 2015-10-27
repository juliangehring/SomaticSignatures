\name{gcContent}
  
\alias{gcContent}

\title{GC Content}

\description{

  Compute the GC content for regions of a reference sequence.

}
  
\usage{
  gcContent(regions, ref)
}

\arguments{

  \item{regions}{GRanges object with the regions for which the GC
    content should be computed.}

  \item{ref}{Reference sequence object, as a 'BSgenome' or 'FaFile'
    object.}

}

\value{

  A numeric vector with the GC content [0,1] for each region.

}
  
\examples{
library(BSgenome.Hsapiens.1000genomes.hs37d5)

regs = GRanges(c("1", "2"), IRanges(1e7, width = 100))

gc = gcContent(regs, BSgenome.Hsapiens.1000genomes.hs37d5)
}
