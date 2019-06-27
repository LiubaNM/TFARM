#' Splits a set of transcription factors in 'present' or 'absent'
#' transcription factors in a set of rules.
#'
#' The function is used to find candidate transcription factor co-regulator
#' with a given transcription factor, by finding the transcription factors that
#' contribute to the prediction of the given trascription factor presence in a
#' set of association rules.
#'
#'
#' @param TFs a string vector: set of transcription factors to check their
#' presence (TF=1) or absence (TF=0) in the set of rules.
#' @param rules data.frame: rules and their quality measures.
#' @param type logical parameter: if \code{type = TRUE}, only rules with
#' all present transcription factors in the left-hand-side are considered
#' (i.e., the left-hand-side of the extracted rules is for example
#' TF1=1, TF2=1, TF3=1).
#'
#' @return A list of two string vectors: the list \code{pres} contains
#' all the transcription factors in \code{TFs} that are present in \code{rules},
#' and the list \code{abs} contains all the transcription factors in \code{TFs}
#' that are absent in \code{rules}.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' # Load r_TEAD4 from the data_man collection of datasets:
#' data('data_man')
#' # Load MCF7_chr1:
#' data('MCF7_chr1')
#'
#' # Transcription factors present in at least one of the regions
#' # in the considered dataset:
#' c <- names(elementMetadata(MCF7_chr1))
#'
#' names(presAbs(c, r_TEAD4, TRUE))
#'
#' # Transcription factors present in at least one of the association rules:
#' p_TFs <- presAbs(c, r_TEAD4, TRUE)$pres
#' p_TFs


presAbs = function(TFs, rules, type) {
    TFs_1 <- sapply(TFs, paste0, "=1")
    TFs_0 <- sapply(TFs, paste0, "=0")
    TFs_v <- c(TFs_1, TFs_0)
    # To find TFs present in at least one rule
    all_io <- paste(as.vector(rules[, 1]), collapse = " ")
    ALL_pres <- sapply(TFs_v, function(x) {
        length(grep(x, all_io))
    })
    if (type == "TRUE") {
        pp <- as.vector(TFs_1[ALL_pres > 0])
        aa <- as.vector(TFs_1[!TFs_1 %in% pp])
    } else {
        # Present TFs:
        pp <- TFs_v[ALL_pres > 0]
        # Absent TFs:
        aa <- TFs_v[!TFs_v %in% pp]
    }
    return(list(pres = pp, abs = aa))
}



