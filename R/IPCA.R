#' Principal Components Analysis for distributions of variations of support,
#' confidence and lift obtained removing a transcription factor
#' from a set of rules.
#'
#' The function is used for validation of the defined Importance Index.
#' It is defined as the linear combination of variations of
#' standardized support, confidence and lift measures.
#' The Principal Component Analysis lets the user evaluate
#' if in the 1D reference system defined by such linear combination,
#' it is possible to describe the variability of the data.
#'
#' @param delta_list list of variations of standardized distributions
#' of support, confidence and lift measures, obtained using the \code{IComp}.
#' @param IMP the importance matrix with the mean Importance Index
#' of every candidate co-regulator transcription factor and the number of rules
#' in which each of them appears.
#' @return Variance explained by every principal component (\code{summary}),
#' scores (i.e., the coordinates) of data in delta_list in the reference system
#' defined by the principal components (\code{scores}) and loadings
#' (i.e., the coefficinets) of the linear combination that defines
#' each principal component (\code{loadings}).
#' The plots of the variability, the cumulate percentage of variance explained
#' by each principal component and loadings of every principal component are
#' also returned.
#'
#' @export
#' @importFrom stats princomp
#' @importFrom graphics layout barplot plot title abline box
#'
#' @examples
#' # Load IMP, DELTA and I from the data:man collection of datasets:
#' data("data_man")
#'
#' colnames(IMP)
#' I <- data.frame(IMP$TF, IMP$imp, IMP$nrules)
#' i.pc <- IPCA(DELTA, I)
#' names(i.pc)


IPCA <- function(delta_list, IMP){
  DZ_2 <- unlist(delta_list)[which(!is.na(unlist(delta_list)))]
  Delta_Z_0 <- data.frame(TF=rep(IMP[,1], IMP[,3]),matrix(DZ_2,ncol=3,
                                                          byrow=TRUE))
  colnames(Delta_Z_0)[2:4] <- c('delta z_s','delta z_c','delta z_l')
  Delta_Z <- Delta_Z_0[,2:4]
  pc.Z <- princomp(Delta_Z, cor = FALSE, scores = TRUE)
  layout(matrix(c(1,2,3,3,4,4,5,5), 4, 2, byrow = TRUE),
         heights=c(1,0.6,0.6,0.6))
  barplot(pc.Z$sdev^2, las=2, main='Variances of the principal components',
          ylab='Variances', ylim=c(0,90), cex.main=1)
  plot(cumsum(pc.Z$sdev^2)/sum(pc.Z$sde^2), type='b', axes=TRUE,
       xlab='number of components', ylab='contribute to total variance')
  box()
  title(main='Cumulate variances of the principal components', cex.main=1)
  scores.Z <- pc.Z$scores
  load.Z <- pc.Z$loadings
  for(i in 1:3)
  {
    barplot(load.Z[,i], ylim = c(-1, 1), ylab=paste('Comp.', i, sep=''))
    abline(h=0)
    if (i == 1) title(main='Loadings of the principal components')
  }
  return(list("summary"=summary(pc.Z),"scores"=scores.Z,"loadings"=load.Z))
}