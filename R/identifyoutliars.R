#' Identify outliars
#'
#' This function takes the output from TwoSampleMR::harmonise_data() function and a list object including
#' identified outliars
#'
#' @param dat Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.
#' @param outliers Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers.
#' @param outlier_correction Defualt = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#' @param outlier_threshold If outlier_correction = "none" then the p-value threshold for detecting outliers is by default 0.05.
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
#' @export
#' @return List
#' dat  Cleaned dat input
#' radialmr  Results from RadialMR analysis
#' outliers  List of outliers used
#' id_list  List of GWAS IDs used
#' search  Result from search of outliers against GWAS IDs
#' candidate_instruments  Instruments for candidate traits
#' candidate_outcome  Extracted instrument SNPs from outcome
#' candidate_outcome_dat  Harmonised candidate - outcome dataset
#' candidate_outcome_mr  MR analysis of candidates against outcome
#' candidate_exposure   Extracted instrument SNPs from exposure
#' candidate_exposure_dat  Harmonised candidate - exposure dataset
#' candidate_exposure_mr  MR analysis of candidates against exposure


identifyOutliers = function(dat, outliers="RadialMR", outlier_correction="none",
                             outlier_threshold=ifelse(outlier_correction=="none", 0.05/nrow(dat), 0.05)) {


  stopifnot(outlier_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  stopifnot(outlier_threshold > 0 & outlier_threshold < 1)

  #Initialise output object
  #This is the tryxscan file which is taken by functions in the mrtryx package
  output <- list()

  #Ensure harmonised dataset only includes one exposure and outcome, otherwise use first one
  if(length(unique(dat$id.exposure)) > 1 | length(unique(dat$id.outcome)) > 1)
  {
    message("Warning! Multiple exposure/outcome combinations found")
    message("Only using first exposure / outcome combination")
  }
  dat <- subset(dat, id.exposure == id.exposure[1] & id.outcome == id.outcome[1] & mr_keep)
  output$dat <- dat

  #If the user wants to use the RadialMR package to identify outliars, make
  #sure package is available
  if(outliers[1] == "RadialMR")
  {

    radialor <- RadialMR::ivw_radial(RadialMR::format_radial(dat$beta.exposure,
                                                             dat$beta.outcome, dat$se.exposure,
                                                             dat$se.outcome, dat$SNP), alpha=1, weights=3)

    # apply outlier_correction method with outlier_threshold to radial SNP-Q statistics
    radialor$data$rdpadj <- p.adjust(radialor$data$Qj_Chi, outlier_correction)
    radial <- radialor
    radial$outliers <- subset(radialor$data, radialor$data$rdpadj < outlier_threshold, select=c(SNP, Qj, Qj_Chi, rdpadj))
    colnames(radial$outliers) = c("SNP", "Q_statistic", "p.value", "adj.p.value")
    rownames(radial$outliers) <- 1:nrow(radial$outliers)

    # apply outlier_correction method with outlier_threshold to radial SNP-Q statistics


    if(radial$outliers[1] == "No significant outliers")
    {
      message("No outliers found")
      message("Try changing the outlier_threshold parameter")
      return(NULL)
    }
    outliers <- as.character(radial$outliers$SNP)
    message("Identified ", length(outliers), " outliers")
    output$radialmr <- radial
    if(length(outliers) == 0)
    {
      message("No outliers found. Exiting")
      return(output)
    }
  } else {
    nout <- length(outliers)
    outliers <- subset(dat, SNP %in% outliers)$SNP
    message(length(outliers), " of ", nout, " of specified outliers present in dat")
    if(length(outliers) == 0)
    {
      message("No outliers found. Exiting")
      return(output)
    }
  }
  output$outliers <- outliers

  return(output)
}


