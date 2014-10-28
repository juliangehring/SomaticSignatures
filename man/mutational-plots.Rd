\name{mutational-plots}

\alias{mutational-plots}
\alias{plotVariantAbundance}

\title{Mutational Plots}

\description{

  Plots for variant analysis

}

\usage{
plotVariantAbundance(x, group = NULL, alpha = 0.5, size = 2)
}

\arguments{

  \item{x}{A VRanges object [required].}

  \item{group}{Grouping variable, refers to a column name in 'x'.  By
    default, no grouping is performed.}

  \item{alpha}{Alpha value for data points.}

  \item{size}{Size value for data points.}

}

\details{

  The 'plotVariantAbundance' shows the variant frequency in relation to
  the total coverage at each variant position.  This can be useful for
  examining the support of variant calls.
  
}

\value{

  A 'ggplot' object.

}
