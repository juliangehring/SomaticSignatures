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

plotMutationSpectrum(vr, group, colorby = c("sample", "alteration"))

plotSignatureMap(s)
plotSignatures(s, normalize = TRUE, percent = FALSE)

plotSampleMap(s)
plotSamples(s, normalize = TRUE, percent = FALSE)
}

\arguments{

  \item{s}{MutationalSignatures object [required]}

  \item{vr}{VRanges object}
  
  \item{colorby}{Which variable to use for the coloring in the spectra
    representation.}

  \item{normalize}{Relative contribution}

  \item{percent}{Display the results as fraction (FALSE) or percent
    (TRUE)}.

  \item{group}{Grouping variable}

}

\details{
  
  With the plotting function, the obtained signatures and their
  occurrance in the samples can be visualized either as a heatmap
  ('plotSignatureMap', 'plotSampleMap') or a barchart ('plotSignature',
  'plotSamples').
  
}

\value{

  A ggplot object
  
}


\seealso{

  'ggplot2' package

}
