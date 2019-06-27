#' Generates an heatmap visualization of the mean Importance Index
#' of pairs of transcription factors.
#'
#' A square matrix (TFs x TFs) is built. The element (i,j) of such matrix
#' contains the mean Importance Index (II) of the couple of transcription
#' factors (TFi, TFj). This is represented in a heatmap visualization,
#' where the scale color from blue to white indicates low mean
#' importance of the couple and from white to red indicates 
#' high mean importance of the couple.
#'
#' @param TFs a string vector with names of transcription factors.
#' @param I a vector of mean II of pairs of transcription factors in TFs.
#'
#' @return The Importance index of pairs of Transcription Factors as a heatmap
#' @export
#' @import fields
#' @import gplots 
#' @importFrom graphics par axis title
#' @examples
#' # Load p_TFs and I_c_2 from the data_man collection of datasets:
#' data('data_man')
#'
#' # Heatmap visualization of the mean importances of transcription factors in p
#' # and their combinations in two elements:
#' heatI(p_TFs, I_c_2)

heatI <- function (TFs, I) 
{
  estract <- function(x) {
    vec <- vapply(TFs, function(TFs) gregexpr(TFs, x)[[1]][1] != 
                    -1, numeric(1))
    return(vec)
  }
  l <- length(TFs)
  aa <- t(sapply(I[, 1], estract))
  matrix_imp <- matrix(0, l, l)
  pos <- apply(aa == 1, 1, function(x) which(x)[1:2])
  matrix_imp[cbind(pos[1, which(!is.na(pos[2, ]))], pos[2, 
                                                        which(!is.na(pos[2, ]))])] <- I[which(!is.na(pos[2, ])), 
                                                                                        "imp"]
  matrix_imp[cbind(pos[2, which(!is.na(pos[2, ]))], pos[1, 
                                                        which(!is.na(pos[2, ]))])] <- I[which(!is.na(pos[2, ])), 
                                                                                        "imp"]
  rownames(matrix_imp) = TFs
  colnames(matrix_imp) = TFs
  im <-heatmap.2(matrix_imp, col=bluered(20),dendrogram = 'row',cexRow = 1.6, 
                 cexCol = 1.6,key.xlab = "", key.title = "", margins =c(8, 8),
                 key.par=list(mar=c(3,0.2,5,1),cex.axis=0.9),
                 trace='none',density.info='none', symm=F,symkey=F,
                 symbreaks = FALSE, scale="none")
  return(im)
}
