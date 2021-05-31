#' Extract candidate trait instruments
#'
#' This functions takes the output from the outlierAssociations function and
#' produces an instrument for each of the significant candidate traits previously
#' identified. The resulting SNP list can be exported using the exportSNPs function
#' to calculate associations with exposures and/or outcomes if these have been
#' calculated from individual-level data.
#'
#'
#' @param dat Output from outlierAssociations() function
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
#'
#' @return List
#' @export
#'
#' @examples
candidateInstruments = function(dat, include_outliers = FALSE) {

  output = dat
  sign_assoc<- subset(output$search, sig)

  message("Finding instruments for candidate traits")

  #Extracting instruments for all candidate traits identified (significantly associated with outliers)
  output$candidate_instruments <- extract_instruments(unique(sign_assoc$id.outcome))

  if(nrow(output$candidate_instruments) == 0)
  {
    message("No instruments available for the candidate traits")
    return(output)
  }

  #Removing outliars from candidate trait SNP list
  if(!include_outliers)
  {
    message("Removing outlier SNPs from candidate trait instrument lists")
    output$candidate_instruments <- dplyr::group_by(output$candidate_instruments, id.exposure) %>%
      dplyr::do({
        x <- .
        y <- subset(sign_assoc, id.outcome == x$id.exposure[1])
        x <- subset(x, !SNP %in% y$SNP)
        x
      })
  }

  if(nrow(output$candidate_instruments) == 0)
  {
    message("No instruments available for the candidate traits")
    return(output)
  }

  message(length(unique(output$candidate_instruments$id.exposure)), " traits with at least one instrument")
  return(output)
}

#' Export candidate SNP list
#'
#' Exports list of SNPs which form part of the instrument for the candidate traits
#' identified. Each line represents one SNP. This file can be read by programs such
#' as bgenix (in BASH) which aid in SNP extraction from larger genomic datasets.
#'
#' @param dat output from candidateInstruments function
#' @param filename filename for SNP list
#'
#' @return
#' @export
#'
#' @examples
exportSNPs = function(dat, filename) {
  stopifnot(is.character(filename))

  SNP_list = unique(dat$candidate_instruments$SNP)
  readr::write_tsv(data.frame(SNP_list),
            file = filename,
            col_names = FALSE)
}
