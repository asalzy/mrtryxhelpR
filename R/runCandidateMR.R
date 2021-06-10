#' Run MRTryx Candidate MR analyses
#'
#' This function takes the output from the candidateInstruments() function
#' and formated data for the candidate-outcome associations and candidate-exposure
#' associations. These formated datasets are produced from the TwoSampleMR::format_data()
#' or the importPlink2() function.
#'
#' @param dat output from the CandidateInstruments function
#' @param candidate_outcome_associations formated data (formated using TwoSampleMR::format_data() or importPlink2()) for candidate SNP - outcome associations
#' @param candidate_exposure_associations formated data (formated using TwoSampleMR::format_data() or importPlink2()) for candidate SNP - exposure associations
#' @param mr_method default is IVW
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
#'
#' @return List
#' @export
#'
#' @examples
runCandidateMR = function(dat,candidate_outcome_associations,
                          candidate_exposure_associations,
                          mr_method = "mr_ivw",
                          include_outliers = FALSE,
                          use_proxies = FALSE) {

  ######
  #This section originally looks up for the instrumenting SNPs identified for the
  #Candidate trait in the outcome GWAS. However, this function allows users to import
  #their own associations for the Candidate-SNP outcome associations

  output = dat

  #Making sure that the phenotype ID is the same for the exposure and outcomes when the new
  #associations are imported
  if(unique(candidate_outcome_associations$id.outcome) != unique(output$dat$id.outcome)) {
    candidate_outcome_associations$id.outcome = unique(output$dat$id.outcome)
  }

  if(unique(candidate_exposure_associations$id.outcome) != unique(output$dat$id.exposure)) {
    candidate_exposure_associations$id.outcome = unique(output$dat$id.exposure)
  }



  message("Looking up candidate trait instruments for ", output$dat$outcome[1])
  output$candidate_outcome <- candidate_outcome_associations
  if(is.null(output$candidate_outcome))
  {
    message("None of the candidate trait instruments available for ", output$dat$outcome[1])
    return(output)
  }
  message(nrow(output$candidate_outcome), " instruments extracted for ", output$dat$outcome[1])

  #Harmonizing candidate trait - outcome dataset
  output$candidate_outcome_dat <- suppressMessages(TwoSampleMR::harmonise_data(output$candidate_instruments,
                                                                               output$candidate_outcome))
  output$candidate_outcome_dat <- subset(output$candidate_outcome_dat, mr_keep)
  if(nrow(output$candidate_outcome_dat) == 0)
  {
    message("None of the candidate trait instruments available for ", output$dat$outcome[1], " after harmonising")
    return(output)
  }

  message("Performing MR of ", length(unique(output$candidate_outcome_dat$id.exposure)), " candidate traits against ", output$dat$outcome[1])

  output$candidate_outcome_mr <- suppressMessages(mr(output$candidate_outcome_dat, method_list=c("mr_wald_ratio", mr_method)))


  ######
  #This section originally looks up for the instrumenting SNPs identified for the
  #Candidate trait in the exposure GWAS. However, this function allows users to import
  #their own associations for the Candidate-SNP expososure associations

  message("Looking up candidate trait instruments for ", output$dat$exposure[1])
  output$candidate_exposure <- candidate_exposure_associations
  if(is.null(output$candidate_exposure))
  {
    message("None of the candidate trait instruments available for ", output$dat$exposure[1])
    return(output)
  }
  message(nrow(output$candidate_exposure), " instruments extracted for ", output$dat$exposure[1])

  output$candidate_exposure_dat <- suppressMessages(harmonise_data(output$candidate_instruments, output$candidate_exposure))
  output$candidate_exposure_dat <- subset(output$candidate_exposure_dat, mr_keep)
  if(nrow(output$candidate_exposure_dat) == 0)
  {
    message("None of the candidate trait instruments available for ", output$dat$exposure[1], " after harmonising")
    return(output)
  }

  message("Performing MR of ", length(unique(output$candidate_exposure_dat$id.exposure)), " candidate traits against ", output$dat$exposure[1])

  output$candidate_exposure_mr <- suppressMessages(mr(output$candidate_exposure_dat, method_list=c("mr_wald_ratio", mr_method)))

  #####
  #This section looks at exposure-candidate trait relationships. Since all SNPs
  #for candidate traits of interest can be found on MR-Base we can still query this and
  #don't have to import our own associations

  message("Looking up exposure instruments for ", length(unique(out2$id.outcome)), " candidate traits")
  output$exposure_candidate <- extract_outcome_data(unique(output$dat$SNP), unique(out2$id.outcome), proxies=use_proxies)
  if(is.null(output$exposure_candidate))
  {
    message("None of the candidate trait instruments available for ", output$dat$exposure[1])
    return(output)
  }
  message(nrow(output$candidate_exposure), " instruments extracted")

  temp <- subset(output$dat, select=c(SNP, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, id.exposure, exposure))

  if(!include_outliers)
  {
    message("Removing outlier SNPs from candidate trait outcome lists")
    n1 <- nrow(output$exposure_candidate)
    output$exposure_candidate <- dplyr::group_by(output$exposure_candidate, id.outcome) %>%
      dplyr::do({
        x <- .
        y <- subset(out2, id.outcome == x$id.outcome[1])
        x <- subset(x, !SNP %in% y$SNP)
        x
      })
    message("Removed ", n1 - nrow(output$exposure_candidate), " outlier SNPs")
  }


  output$exposure_candidate_dat <- suppressMessages(harmonise_data(temp, output$exposure_candidate))
  output$exposure_candidate_dat <- subset(output$exposure_candidate_dat, mr_keep)
  if(nrow(output$exposure_candidate_dat) == 0)
  {
    message("None of the candidate traits have the ", output$dat$exposure[1], " instruments after harmonising")
    return(output)
  }

  message("Performing MR of ", output$dat$exposure[1], " against ", length(unique(output$exposure_candidate_dat$id.outcome)), " candidate traits")

  output$exposure_candidate_mr <- suppressMessages(mr(output$exposure_candidate_dat, method_list=c("mr_wald_ratio", mr_method)))

  return(output)

}
