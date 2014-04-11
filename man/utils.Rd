\name{variants-utils}

\alias{variants-utils}
\alias{dfConvertColumns}

\title{Utility functions}

\description{
  Utility functions
}

\usage{
dfConvertColumns(x, from = "character", to = "factor")
}

\arguments{

  \item{x}{A 'data.frame' to convert [required].}

  \item{from}{The class of the columns to be converted [default:
    'character'].}

  \item{to}{The class of the columns to be converted to [default: 'factor'].}

}

\details{

  The 'dfConvertColumns' converts all columns of a data frame with class
  'from' to the class 'to'.
  
}

\value{

  A 'data.frame' object.

}
