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
#' @param search_correction Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#' @param search_threshold If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.
#' @param id_list The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes().
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
#' @param mr_method Method to use for candidate trait - exposure/outcome analysis. Default is mr_ivw. Can also provide basic MR methods e.g. mr_weighted_mode, mr_weighted_median etc. Also possible to use "strategy1" which performs IVW in the first instance, but then weighted mode for associations with high heterogeneity.



identify_outliars = function(dat, outliers="RadialMR", outlier_correction="none",
                             outlier_threshold=ifelse(outlier_correction=="none", 0.05/nrow(dat), 0.05),
                             use_proxies=FALSE, search_correction="none",
                             search_threshold=ifelse(search_correction=="none", 5e-8, 0.05),
                             id_list="default", include_outliers=FALSE, mr_method="mr_ivw") {

  #Check function inputs
  stopifnot(search_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  stopifnot(search_threshold > 0 & search_threshold < 1)

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

  stopifnot(length(mr_method) == 1)
  stopifnot(mr_method %in% mr_method_list()$obj | mr_method == "strategy1")


  #If the user wants to use the RadialMR package to identify outliars, make
  #sure package is available
  if(outliers[1] == "RadialMR")
  {
    message("Using RadialMR package to detect outliers")
    cpg <- require(RadialMR)
    if(!cpg)
    {
      stop("Please install the RadialMR package\ndevtools::install_github('WSpiller/RadialMR')")
    }
    cpg <- require(dplyr)
    if(!cpg)
    {
      stop("Please install the RadialMR package\ndevtools::install_github('WSpiller/RadialMR')")
    }



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

