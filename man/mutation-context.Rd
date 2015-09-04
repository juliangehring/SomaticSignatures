\name{mutationContext}

\alias{mutationContext}
\alias{mutationContextMutect}
\alias{mutationContextH5vc}

\title{mutationContext functions}

\description{
  
  Extract the sequence context surrounding SNVs from a genomic
  reference.
  
}

\usage{
mutationContext(vr, ref, k = 3, strand = FALSE, unify = TRUE, check = FALSE)
mutationContextMutect(vr, k = 3, unify = TRUE)
mutationContextH5vc(vc, ms, unify = TRUE)
}

\arguments{

  \item{vr}{'VRanges' with SNV substitutions, with 'ref' and 'alt'
    columns filled [required].  Each element of 'ref' and 'alt' have be
    a single base from the DNA bases (A,C,G,T).  For
    'mutationContextMutect', an object as returned by the 'readMutect'
    function.}

  \item{ref}{A 'BSgenome', 'FaFile' or 'TwoBitfile' object representing
    the reference sequence [required].  More generally, any object with
    a defined 'getSeq' method can be used.}

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
  genomic sequence with the 'mutationContext' function.

  Different types of classes that represent the genomic sequence can
  used togther with the 'mutationContext' function: 'BSgenome',
  'FastaFile' and 'TwoBitFile' objects are supported through
  Bioconductor by default. See the vignette for examples discussing an
  analysis with non-referene genomes.
  
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

  \code{\link{readMutect}} for mutationContextMutect
  
  'showMethods("getSeq")' for genomic references that can be used
}

\examples{

  mutect_path = system.file("examples", "mutect.tsv", package = "SomaticSignatures")
  vr1 = readMutect(mutect_path)
  ct1 = mutationContextMutect(vr1)

}
 
\keyword{manip}
