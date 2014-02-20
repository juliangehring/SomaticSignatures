
\name{GRanges-converters}

\alias{grangesPlain}
\alias{ncbi}
\alias{ucsc}
\alias{seqchar}

\title{GRanges converter functions}

\description{

  A set of utilities functions to convert and extract data in 'GRanges'
  objects.

}

\usage{
grangesPlain(x)
ncbi(x)
ucsc(x)
seqchar(x)
}

\arguments{

  \item{x}{A 'GRanges' object or one inheriting from the 'GRanges' class
  [required].}

}

\details{

  \itemize{

    \item{grangesPlain}{Extracts only the 'GRanges' information by
      dropping the metadata columns of the object.  The 'seqinfo' slot
      is kept.}

    \item{ncbi, ucsc}{Shorthand for converting the seqnames notation to
      'UCSC' (e.g. 'chr1', 'chrM') or 'NCBI' (e.g. '1', 'MT'') notation,
      respectively.  This also sets the 'genome' slot in the 'seqinfo'
      field to 'NA'.}

    \item{seqchar}{Extracts the 'seqnames' as a character vector.}

  }
}

\value{

  For 'grangesPlain': A GRanges object without metadata.

  For 'ncbi', 'ucsc': An object of the same class as the input.

  For 'seqchar': A character vector with 'seqnames'.

}

\seealso{
  
  'GenomicRanges' package: 'seqnameStyle', 'seqnames', 'mcols'

}

\examples{

mutect_path = system.file("examples", "mutect.tsv", package = "SomaticSignatures")
vr1 = readMutect(mutect_path, strip = TRUE)

## extract the GRanges
gr = grangesPlain(vr1)

## convert back and forth
gr_ncbi = ncbi(gr)
gr_ucsc = ucsc(gr_ncbi)

identical(gr, gr_ucsc)

## extract the seqnames as a character vector
seq_chars = seqchar(gr)

}

\keyword{manip}
\keyword{utilities}
