#' Boxplots visualization of the Importance Index distribution
#' of a set of transcription factors.
#'
#'  For a set of candidate co-regulator transcription factors,
#'  the Importance Index distribution of each transcription factor
#'  is plotted in boxplots. The shape of every boxplot depends on:
#'  the dimension of the distribution, which is equal to the number of rules
#'  in which each transcription factor appears (the higher is such number,
#'  the larger is a boxplot), the variability of the distribution
#'  (the higher is the variability of the Importance Index distribution,
#'  the longer is a boxplot). Moreover, the higher is the median of
#'  the Importance Index distribution for a candidate co-regulator
#'  transcription factor, the higher the boxplot is aligned with respect to
#'  the y-axis.
#'
#' @param I a list of Importance Index distributions (IIDs)
#' @param TFs string vector with the names of the transcription factors
#'
#' @return A boxplot representing the transcription factor IIDs
#' @export
#' @importFrom graphics boxplot par title
#'
#' @examples
#' # Load IMP_Z and p from the data_man collection of datasets:
#' data('data_man')
#'
#' # Plot the Importance Index distributions of the transcription factors in p:
#' distribViz(IMP_Z,p_TFs)

distribViz <- function(I, TFs) {
    lungh <- sapply(I, function(x) {
        length(x)
    })
    ll = unlist(I)
    df <- data.frame(importance = ll, TF = rep(TFs, lungh))
    bp <- boxplot(df$importance ~ df$TF, xlab="",varwidth = TRUE,par(mar = c(8,
                                                                      3, 3, 1)),cex.lab = 0.8, las = 2)
    title(main = "Transcription Factor Importance index distributions")
    return(bp)
}
