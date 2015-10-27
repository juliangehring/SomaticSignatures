\name{kmerFrequency}

\alias{kmerFrequency}

\title{Kmer Frequency}

\description{

  Estimate the occurance frequency of k-mers in a reference sequence.

}

\usage{
kmerFrequency(ref, n = 1e4, k = 1, ranges = as(seqinfo(ref), "GRanges"))
}


\arguments{

  \item{ref}{A 'BSgenome' or 'FaFile' object matching the respective
    reference sequence [required].}

  \item{n}{The number of samples to draw [integer, default: 1e4].}

  \item{k}{The 'k'-mer size of the context, including the variant
    position [integer, default: 3].}

  \item{ranges}{Ranges in respect to the reference sequence to sample
    from [GRanges, default: take from the 'seqinfo' slot].}  }

\value{

  A named vector, with names corresponding the the k-mer and value to
  the frequency.

}

\details{

  The k-mer frequency is estimated by random sampling of 'n' locations
  across the specified 'ranges' of the reference sequence.
  
}

\examples{
  library(BSgenome.Hsapiens.1000genomes.hs37d5)
  
  kmer_freq = kmerFrequency(BSgenome.Hsapiens.1000genomes.hs37d5, 1e2, 3)
}
