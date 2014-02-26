\name{SomaticSignatures}

\alias{SomaticSignatures}
\alias{SomaticSignatures-package}

\title{SomaticSignatures package}

\description{

  Identifying somatic signatures of single nucleotide variants.

}

\details{

  The 'SomaticSignatures' package offers the framework for identifying
  mutational signatures of single nucleotide variants (SNVs) from
  high-throughput experiments.  In the concept of mutational signatures,
  a base change resulting from an SNV is regarded in term of motifs
  which embeds the variant in the context of the surrounding genomic
  sequence.  Based on the frequency of such motifs across samples,
  mutational signatures and their occurrance in the samples can be
  estimated.  An introduction into the methodology and a use case are
  illustrated in the vignette of this package.

}

\author{

  Julian Gehring, with contributions from Bernd Fischer and Wolfgang
  Huber.

  Maintainer: Julian Gehring, EMBL Heidelberg <julian.gehring@embl.de>

}

\references{

Nik-Zainal, Serena, Ludmil B. Alexandrov, David C. Wedge, Peter Van Loo,
Christopher D. Greenman, Keiran Raine, David Jones, et al. "Mutational
Processes Molding the Genomes of 21 Breast Cancers." Cell 149, no. 5
(May 25, 2012): 979-993. doi:10.1016/j.cell.2012.04.024.

Alexandrov, Ludmil B., Serena Nik-Zainal, David C. Wedge, Samuel
A. J. R. Aparicio, Sam Behjati, Andrew V. Biankin, Graham R. Bignell, et
al. "Signatures of Mutational Processes in Human Cancer." Nature 500,
no. 7463 (August 22, 2013): 415-421. doi:10.1038/nature12477.

Gaujoux, Renaud, and Cathal Seoighe. "A Flexible R Package for
Nonnegative Matrix Factorization." BMC Bioinformatics 11, no. 1 (July 2,
2010): 367. doi:10.1186/1471-2105-11-367.

Stacklies, Wolfram, Henning Redestig, Matthias Scholz, Dirk Walther, and
Joachim Selbig. "pcaMethods - A Bioconductor Package Providing PCA Methods
for Incomplete Data." Bioinformatics 23, no. 9 (May 1, 2007):
1164-1167. doi:10.1093/bioinformatics/btm069.

}

\examples{

library(SomaticSignatures)

vignette(package = "SomaticSignatures")

}

\keyword{package}
