\name{motif-functions}

\alias{motifMatrix}

\title{Group somatic moticfs}

\description{

  Tabulate somatic motifs by a grouping variable.

}

\usage{
motifMatrix(x, group = "sample", normalize = TRUE)
}

\arguments{

  \item{x}{GRanges object [required]}

  \item{group}{Grouping variable name [character, default: 'sample']}

  \item{normalize}{Normalize to frequency}

}

\details{

  The 'motifMatrix' function transforms the metadata columns of a
  'VRanges' object, as returned by the 'mutationContext' function, to a
  matrix of the form 'motifs x groups'.  This constitutes the bases for
  the estimation of the signatures.  By default (with 'normalize' set to
  TRUE), the counts are transformed to frequencies, such that the sum of
  frequencies of each group equal 1.  Otherwise (with 'normalize' set to
  FALSE), the counts for each mofis in a group is returned.
  
}

\seealso{

  'mutationContext', 'mutationContextMutect'

}
