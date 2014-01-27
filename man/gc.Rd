
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

\seealso{

  Inspired by the 'getGCcontent' function of the 'exomeCopy' package.
  
}

\examples{
library(BSgenome.Hsapiens.UCSC.hg19)

regs = GRanges(c("chr1", "chr2"), IRanges(1e7, width = 100))

gc = gcContent(regs, BSgenome.Hsapiens.UCSC.hg19)
}
