\name{readMutect}

\alias{readMutect}

\title{readMutect}

\description{

  Import 'mutect' calls.

}

\usage{
readMutect(file, columns, strip = FALSE)
}

\arguments{

  \item{file}{Location of the mutect tsv files [character, required]}

  \item{columns}{Names of columns to import from the file [character
    vector, optional, default: missing].  If missing, all columns will
    be imported.}

  \item{strip}{Should additional columns be imported?  [logical,
    default: FALSE]. If TRUE, return only the bare 'VRanges' object.}

}

\value{

  A 'VRanges' object, with each row corresponding to one variant in the
  original file.
  
}

\details{

  The 'readMutect' functions imports the mutational calls of a '*.tsv'
  file returned by the 'mutect' caller to a 'VRanges' object.  For a
  description of the information of the columns, please refer to the
  mutect documentation.

}

\references{

  Cibulskis, Kristian, Michael S. Lawrence, Scott L. Carter, Andrey
  Sivachenko, David Jaffe, Carrie Sougnez, Stacey Gabriel, Matthew
  Meyerson, Eric S. Lander, and Gad
  Getz. "Sensitive Detection of Somatic Point Mutations in Impure and Heterogeneous Cancer Samples."
  Nature Biotechnology advance online publication (February 10,
  2013). doi:10.1038/nbt.2514.

  \url{http://www.broadinstitute.org/cancer/cga/mutect_run}

}

\examples{

mutect_path = system.file("examples", "mutect.tsv", package = "SomaticSignatures")
vr1 = readMutect(mutect_path)
vr2 = readMutect(mutect_path, strip = TRUE)

}

\keyword{IO}
