% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DELTA.R
\docType{data}
\name{DELTA}
\alias{DELTA}
\title{Contains the delta variations of support,
confidence and lift.}
\format{An object of class \code{"list"}}
\usage{
data("data_man")
}
\description{
DELTA is a list of 12 elements and each element has two columns
representing support (\code{diff_supp_Z}) and confidence (\code{diff_conf_Z})
respectively. It is included in the 
\code{data_man} collection.
}
\examples{
# DELTA is found in the data_man collection of datasets:
data("data_man")
head(DELTA[[1]])
}
\keyword{datasets}
