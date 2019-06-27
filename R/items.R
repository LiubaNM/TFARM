#' Extracts items in the left-hand-side of an association rule.
#'
#' The function is used in the elaboration of the left-hand-side of
#' association rules to search for the items.
#'
#' @param itemset object of class rules.
#' @return A string vector with the items in the left-hand-side
#' of the rule in itemset.
#' @export
#' @import stringr
#'
#' @examples
#' items('{TAF1=1,EP300=1,MAX=1}')
#' # the output is: 'TAF1=1','EP300=1','MAX=1'

items <- function(itemset) {
    items.0 <- gsub("[:{:]", "", itemset)
    items.1 <- gsub("[:}:]", "", items.0)
    items <- strsplit(items.1, ",")
    return(items[[1]])
}
