\name{signature-plots}

\alias{plotObservedSpectrum}
\alias{plotFittedSpectrum}
\alias{plotMutationSpectrum}

\alias{plotSignatureMap}
\alias{plotSignatures}
\alias{plotSampleMap}
\alias{plotSamples}

\title{Plot Mutational Signatures}

\description{

  Visualize estimated signatures, sample contribution, and mutational
  spectra.

}

\usage{
plotObservedSpectrum(s, colorby = c("sample", "alteration"))
plotFittedSpectrum(s, colorby = c("sample", "alteration"))

plotMutationSpectrum(vr, group, colorby = c("sample", "alteration"), normalize = TRUE)

plotSignatureMap(s)
plotSignatures(s, normalize = FALSE, percent = FALSE)

plotSampleMap(s)
plotSamples(s, normalize = FALSE, percent = FALSE)
}

\arguments{

  \item{s}{MutationalSignatures object [required]}

  \item{vr}{VRanges object}
  
  \item{colorby}{Which variable to use for the coloring in the spectra
    representation.}

  \item{normalize}{Plot relative constributions (TRUE) instead of
    absolute (FALSE) ones.}

  \item{percent}{Display the results as fraction (FALSE) or percent
    (TRUE)}.

  \item{group}{Charactering string that represents the variable name
  used for grouping.}

}

\details{
  
  With the plotting function, the obtained signatures and their
  occurrance in the samples can be visualized either as a heatmap
  ('plotSignatureMap', 'plotSampleMap') or a barchart ('plotSignature',
  'plotSamples').

  Since the plotting is based on the 'ggplot2' framework, all properties
  of the plots can be fully controlled by the user after generating the
  plots.  Please see the examples for some customizations and the
  'ggplot2' documentation for the entire set of options.

}

\value{

  A 'ggplot' object, whose properties can further be changed
  
}


\seealso{

  See the \code{\link{ggplot2}} package for customizing the plots.

}


\examples{
data("sca_sigs", package = "SomaticSignatures")

plotSamples(sigs_nmf)

plotSignatures(sigs_nmf, normalize = TRUE)

## customize the plots ##
p = plotSamples(sigs_nmf)

library(ggplot2)
## (re)move the legend
p = p + theme(legend.position = "none")
## change the axis labels
p = p + xlab("Studies")
## add a title
p = p + ggtitle("Somatic Signatures in TGCA WES Data")
## change the color scale
p = p + scale_fill_brewer(palette = "Blues")
## decrease the size of x-axis labels
p = p + theme(axis.text.x = element_text(size = 9))

p
}
