\name{scaSNVRanges}

\alias{scaSNVRanges}

\title{SNV VRanges from SCA dataset}

\description{

  Create VRanges for somatic SNV calls in the SomaticCancerAlterations
  dataset.

}

\usage{
  scaSNVRanges(chrs = hsAutosomes())
}

\arguments{

  \item{chrs}{Chromosomes to include in the results.  Defaults to human
    autosomes.}

}

\value{

  A 'VRanges' object with somatic SNV calls.
  
}


\seealso{
  
  'SomaticCancerAlterations' package
  
}
