% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TF_Imp.R
\docType{data}
\name{TF_Imp}
\alias{TF_Imp}
\title{Contains the candidate co-regulators and the number of rules
associated with them.}
\format{An object of class \code{"data.frame"}}
\usage{
data("data_man")
}
\description{
Within the \code{data_man} data collection, the dataset TF_Imp has 3 columns
and 12 rows:
the fist column contains the transcription factors (\code{IMP.TF}),
the Importance Indexes associated with each trascription factor
are listed in the second column (\code{IMP.imp}) and the third column
contains the number of rules found for each transcription factor
(\code{IMP.nrules}).
}
\examples{
# TF_Imp is found in the data_man collection of datasets:
data("data_man")
head(TF_Imp$IMP.TF)
}
\keyword{datasets}
