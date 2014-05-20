\name{mutationContext}

\alias{mutationContext}
\alias{mutationContextMutect}
\alias{mutationContextH5vc}

\title{mutationContext functions}

\description{
  
  Extract the sequence context surrounding a SNV from a genomic
  reference.
  
}

\usage{
mutationContext(vr, ref, k = 3, strand = FALSE, unify = TRUE, check = TRUE)
mutationContextMutect(vr, k = 3, unify = TRUE)
mutationContextH5vc(vc, ms, unify = TRUE)
}
  
\arguments{
  
  \item{vr}{'VRanges' object, with 'ref' and 'alt' columns filled
    [required].  For 'mutationContextMutect', an object as returned by
    the 'readMutect' function.}
  
  \item{ref}{A 'BSgenome' or 'FaFile' object representing the reference
    sequence [required].  More generally, any object with a defined
    'getSeq' method can be used.}
  
  \item{k}{The 'k'-mer size of the context, including the variant
    position [integer, default: 3].  The variant will be located at the
    middle of the k-mer which requires 'k' to be odd.}
  
  \item{strand}{Should all variants be converted to the 'plus'
    strand? [logical, default: FALSE].}
  
  \item{unify}{Should the alterations be converted to have a C/T base
    pair as a reference alleles? [logical, default: TRUE]}
  
  \item{check}{Should the reference base of 'vr' be checked against
    'ref' [logical, default: TRUE]?  In case the two references do not
    match, a warning will be printed.}
  
  \item{vc}{A 'DataFrame' object as returned from a variant calling
    analysis by 'h5vc::h5dapply'.  See the 'details' section for more
    information.}
  
  \item{ms}{A 'DataFrame' object as returned by
    'h5vc::mutationSpectrum'.  See the 'details' section for more
    information.}

}

\details{
  
  The somatic motifs of a SNV, composed out of (a) the base change and
  (b) the sequence context surrounding the variant, is extracted from a
  reference sequence with the 'mutationContext' function.
  
  For mutect variant calls, all relevant information is already
  contained in the results and somatic motifs can constructed by using
  the 'mutationContextMutect' function, without the need for the
  reference sequence.

  For h5vc variant calls, the information is merged from the outputs of
  the 'h5dapply' and 'mutationSpectrum' functions of the 'h5vc' package.
  A detailed example is shown in the vignette of the package.
  
}

\value{
  
  The original 'VRanges' object 'vr', with the additional columns
  
  \item{alteration}{DNAStringSet with 'ref|alt'.}
  
  \item{context}{DNAStringSet with '..N..' of length 'k', where N
    denotes the variant position.}
  
}

\seealso{

  'readMutect' for 'mutationContextMutect'

  'mutationSpectrum' from the 'h5vc' package for 'mutationContextH5vc'
  
}

\examples{

  mutect_path = system.file("examples", "mutect.tsv", package = "SomaticSignatures")
  vr1 = readMutect(mutect_path)
  ct1 = mutationContextMutect(vr1)


  library(h5vc)
  library(h5vcData)

  tally_file = system.file("extdata", "example.tally.hfs5", package = "h5vcData")
  data("example.variants", package = "h5vcData")
  ms = mutationSpectrum(variantCalls, tally_file, "/ExampleStudy", context = 1)
  vc = variantCalls

  vr2 = mutationContextH5vc(vc, ms)
  ct2 = motifMatrix(vr2, "sampleNames")

}
 
\keyword{manip}
