#' Computes the Importance Index of a transcription factor
#' in a set of association rules.
#'
#' Given an association rule and a transcription factor \code{TFi},
#' it is evaluated the contribution of \code{TFi} in the rule for the prediction
#' of the presence of the item in the right-hand-side of the rule.
#' Since this contribution is evaluated based on the variations of
#' support and confidence of the rule, the user can visualize
#' such variations by setting the parameter \code{figures = TRUE}.
#'
#' @param TFi string or string vector: the transcription factor (or combination
#' of transcription factors) whose importance distribution is evaluated.
#' @param rules_TF a set of rules in which \code{TFi} is present.
#' @param rules_noTF a set of rules obtained from rules_TF removing
#' the transcription factor (or combination of transcription factors)
#' in TFi (obtained with the function \code{\link{rulesTF0}}).
#' @param figures logical; if \code{figures = TRUE}, graphics with support
#' and confidence distributions of the rules in \code{rulesTF} and
#' \code{rulesTF0} are returned.
#'
#' @return A list of four elements: the \code{imp} element is a vector of
#' doubles with the importances of TFi in each rule in \code{rulesTF};
#' the \code{delta} element of the list is a list with variations of
#' distributions of the two measures of support and 
#' confidence. This output is used in the function \code{\link{IPCA}}
#' for the Principal Component Analysis of such distributions. The \code{rwi}
#' element is a data.frame with the rules in \code{rulesTF}
#' in which the transcription factor \code{TFi} is present 
#' and the \code{rwo} element is a data.frame with rules in \code{rwi}
#' obtained removing the transcription factor \code{TFi}. Furthermore, if the
#' input argument \code{figures} is set to TRUE, also the plots of the distributions
#' of support and confidence of the rules before and after removing the
#' transcription factor TFi are provided.
#' @export
#' @importFrom stats var
#' @importFrom graphics par matplot legend title
#'
#' @examples
#' # Load r_FOSL2 and r_noFOSL2 from the data_man collection of datasets:
#' data('data_man')
#'
#' # The Importance Indexes of FOSL2=1 in the set of rules r_FOSL2 are given by:
#' IComp('FOSL2=1', r_FOSL2, r_noFOSL2, figures=TRUE)


IComp <- function (TFi, rules_TF, rules_noTF, figures)
{
  if (length(TFi) > 1) {
    TFi <- TFi[[1]]
    for (i in 2:length(TFi)) TFi <- paste(TFi, TFi[[i]],
                                          sep = ",")
  }
  both <- cbind(rules_TF, rules_noTF)
  if (all(is.na(both)) == TRUE) {
    return(list(imp = NA, delta = NA))
  }
  else {
    diff_supp <- both[, 3] - both[, 8]
    diff_conf <- both[, 4] - both[, 9]
    diffs <- data.frame(both[, 1], both[, 6], diff_supp,
                        diff_conf)
    colnames(diffs)[1:2] <- colnames(both)[c(1, 6)]
  }
  imp_Z_rule_0 <- apply(diffs[, 3:4], 1, sum)
  m <- max(rules_TF$support)
  n <- dim(diffs)[1]
  if (figures == TRUE) {
    TF_i <- paste0(strsplit(TFi, "=")[[1]][1], "=0")
    par(mfrow = c(1, 2))
    matplot(both[, c(3, 8)], type = "l", lwd = 2, ylab = "support", xlab = "rules ID", xaxt="n", ylim = c(0, 1), cex.lab=1.3)
    axis(1, at=1:dim(rules_TF)[1], labels=row.names(rules_TF), las=1)
    legend("topright", c(paste("with", TFi), paste("with", TF_i)), col = c("black","red"), lty = c(1, 2), lwd = c(2, 2), cex=0.8)
    matplot(both[, c(4, 9)], type = "l", lwd = 2, ylab = "confidence",xlab = "rules ID", xaxt="n", ylim = c(0, 1), cex.lab=1.3)
    axis(1, at=1:dim(rules_TF)[1], labels=row.names(rules_TF), las=1)
    legend("topright", c(paste("with", TFi),paste("with", TF_i)), col = c("black","red"), lty = c(1, 2), lwd = c(2, 2), cex=0.8)
  }
  d_Z <- diffs[, 3:4]
  rwi <- rules_TF
  rwo <- rules_noTF
  return(list(imp = imp_Z_rule_0, delta = d_Z, rwi = rwi, rwo = rwo))
}
