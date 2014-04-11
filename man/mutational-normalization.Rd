\name{mutational-normalization}

\alias{mutational-normalization}
\alias{normalizeMotifs}

\title{Normalize Somatic Motifs}

\description{

  Normalize somatic motifs, to correct for biases between samples.
  
}

\usage{
normalizeMotifs(x, norms)
}

\arguments{

  \item{x}{Matrix, as returned by 'mutationContextMatrix' [required]}

  \item{norms}{Vector with normalization factors [required].  The names
    must match the base sequence names in 'x'.}

}

\value{

  A matrix as 'x' with normalized counts.
  
}

\seealso{

  mutationContextMatrix
  
}
