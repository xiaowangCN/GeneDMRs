% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Correlation_plot.R
\name{Correlation_plot}
\alias{Correlation_plot}
\title{Plot the methylation correlation.}
\usage{
Correlation_plot(inputmethfile_QC, unmeth_exclude = TRUE)
}
\arguments{
\item{inputmethfile_QC}{refers to the input file with methylation levels, with default inputmethfile after quality control.}

\item{unmeth_exclude}{refers to whether to exclude the unmethylated sites or regions, with default TRUE.}
}
\value{
Outputs a correlation figure.
}
\description{
This function outputs the correlation plot for the methylation level of different samples or groups based on R package corrplot.
}
\examples{
Correlation_plot(inputmethfile_QC)
Correlation_plot(siteall)
Correlation_plot(regiongenealls)
Correlation_plot(genefeatureall_cpgfeature)
Correlation_plot(genefeatureall_cpgfeature, unmeth_exclude = FALSE)

}
