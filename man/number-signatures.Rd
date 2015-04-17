\name{numberSignatures}

\alias{assessNumberSignatures}
\alias{plotNumberSignatures}

\title{Number of Signatures}

\description{

  Assessment of the number of signatures in the data.

}

\usage{
assessNumberSignatures(m, nSigs, decomposition = nmfDecomposition, ..., nReplicates = 1)

plotNumberSignatures(gof)
}

\arguments{

  \item{m}{Mutational spectrum matrix, same as used for
    'identifySignatures'.}

  \item{nSigs}{Vector of integers with the numbers of signatures that
    should be tested.  See the 'nSigs' arugment for
    'identifySignatures'.}

  \item{decomposition}{Function to apply for the matrix
    decompositon. See the 'decomposition' argument for
    'identifySignatures'.}

  \item{...}{Additional arguments passed to the 'decomposition'
    function.  See the '...' argument for 'identifySignatures'.}

  \item{nReplicates}{How many replicates should be used for assing a
    single value of 'nSigs'?  For decomposition methods with random
    seeding, values greater than 1 are reasonable.}

  \item{gof}{Data frame, as returned of 'assessNumberSignatures'}.

}

\details{
  
  Compute the decompositon for a given number of signatures, and assess
  the goodness of the reconstruction between the observed and fitted
  mutational spectra M and V, respectively.  The residual sum of squares
  (RSS)

  \deqn{RSS = \sum_{i,j} (M_{ij} - V_{ij})^2}

  and the explained variance

  \deqn{evar = 1 - \frac{RSS}{\sum_{i,j} V_{ij}^2} }{evar = 1 - RSS/sum
    V_{ij}^2}
  
  are used as summary statistics which can generally applied to all
  decomposition approaches.

  The 'plotNumberSignature' function visualizes the results of the
  'assessNumberSignatures' step.  Statistics of the indivdual replicates
  are shown as black crosses, whereas the mean and the range across the
  replicates are depicted in red.

  In practice, these summary statisics should not be trusted blindly,
  but rather interpreted together with biological knowledge and
  scientifc reasoning.  For a discussion of the interpretation of these
  statistics with special focus on the NMF decomposition, please refer
  to the references listed below.
  
}

\value{

  - assessNumberSignatures: A data frame with the RSS and KL statistic
    for each run

  - plotNumberSignatures: A ggplot object

}


\seealso{

  identifySignatures

  The 'rss' and 'evar' functions of the 'NMF' package.
  
}


\references{

  Hutchins LN, Murphy SM, Singh P and Graber JH (2008):
  'Position-dependent motif characterization using non-negative matrix
  factorization.'  Bioinformatics,
  \url{http://dx.doi.org/10.1093/bioinformatics/btn526}
  
}


\examples{

  data("sca_mm", package = "SomaticSignatures")
  
  nSigs = 2:8
  stat = assessNumberSignatures(sca_mm, nSigs, nReplicates = 3)

  plotNumberSignatures(stat)
}
