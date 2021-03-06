% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IComp.R
\name{IComp}
\alias{IComp}
\title{Computes the Importance Index of a transcription factor
in a set of association rules.}
\usage{
IComp(TFi, rules_TF, rules_noTF, figures)
}
\arguments{
\item{TFi}{string or string vector: the transcription factor (or combination
of transcription factors) whose importance distribution is evaluated.}

\item{rules_TF}{a set of rules in which \code{TFi} is present.}

\item{rules_noTF}{a set of rules obtained from rules_TF removing
the transcription factor (or combination of transcription factors)
in TFi (obtained with the function \code{\link{rulesTF0}}).}

\item{figures}{logical; if \code{figures = TRUE}, graphics with support
and confidence distributions of the rules in \code{rulesTF} and
\code{rulesTF0} are returned.}
}
\value{
A list of four elements: the \code{imp} element is a vector of
doubles with the importances of TFi in each rule in \code{rulesTF};
the \code{delta} element of the list is a list with variations of
distributions of the two measures of support and 
confidence. This output is used in the function \code{\link{IPCA}}
for the Principal Component Analysis of such distributions. The \code{rwi}
element is a data.frame with the rules in \code{rulesTF}
in which the transcription factor \code{TFi} is present 
and the \code{rwo} element is a data.frame with rules in \code{rwi}
obtained removing the transcription factor \code{TFi}. Furthermore, if the
input argument \code{figures} is set to TRUE, also the plots of the distributions
of support and confidence of the rules before and after removing the
transcription factor TFi are provided.
}
\description{
Given an association rule and a transcription factor \code{TFi},
it is evaluated the contribution of \code{TFi} in the rule for the prediction
of the presence of the item in the right-hand-side of the rule.
Since this contribution is evaluated based on the variations of
support and confidence of the rule, the user can visualize
such variations by setting the parameter \code{figures = TRUE}.
}
\examples{
# Load r_FOSL2 and r_noFOSL2 from the data_man collection of datasets:
data('data_man')

# The Importance Indexes of FOSL2=1 in the set of rules r_FOSL2 are given by:
IComp('FOSL2=1', r_FOSL2, r_noFOSL2, figures=TRUE)
}
