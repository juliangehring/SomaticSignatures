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

\note{

  While the 'scaSNVRanges' is provided for a convinient access to the
  data of the 'SomaticCancerAlterations' package, we encourage you to
  develop an understanding about the underlying data and its conversion
  to a 'VRanges' object.  

}
  

\value{

  A 'VRanges' object with somatic SNV calls.
  
}


\seealso{
  
  SomaticCancerAlterations package
  
}
