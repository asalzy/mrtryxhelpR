#' Outlier Associations
#'
#' Function finds associations of outliers with candidate traits from MR-Base (https://www.mrbase.org/) and
#' adds these to an output list readable by the mrtryx package
#'
#' @param dat Output from IdentifyOutliers function
#' @param id_list The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes().
#' @param search_correction Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#' @param search_threshold If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.
#' @param min_sample_size When selecting GWAS studies from MR-Base, state the minimum sample size acceptable. Default is 5000.
#' @param min_SNP_numberWhen When selecting GWAS studies from MR-Base, state the minimum SNP number acceptable. Default is 5000.
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
#'
#' @return List
#' @export
#'
#' @examples
outlierAssociations = function(dat, id_list = "default", search_correction="none",
                               search_threshold=ifelse(search_correction=="none", 5e-8, 0.05),
                               min_sample_size = 5000,
                               min_SNP_number = 100000,
                               use_proxies = FALSE){

  output = dat
  #Check function inputs
  stopifnot(search_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  stopifnot(search_threshold > 0 & search_threshold < 1)

  if(id_list[1] == "default")
  {
    ao <- suppressMessages(available_outcomes())

    #Filter searched GWAS based on sample size, SNP number and excluding exposure
    #and outcome variables. Removing duplicated traits.
    ids <- subset(ao) %>%
      dplyr::arrange(desc(sample_size)) %>%
      dplyr::filter(!duplicated(trait),
             mr == 1,
             !grepl("ukb-a", id),
             !grepl("ebi-a", id),
             !id %in% c(dat$id.exposure[1], dat$d.outcome[1]),
             nsnp > min_SNP_number,
             sample_size >min_sample_size)

    id_list <- ids$id
    message("Using default list of ", nrow(ids), " traits")
  }
  output$id_list <- id_list

  output$search <- extract_outcome_data(output$outliers, output$id_list, proxies=use_proxies)
  padj <- p.adjust(output$search$pval.outcome, search_correction)
  output$search$sig <- padj < search_threshold
  out2 <- subset(output$search, sig)
  if(nrow(out2) == 0)
  {
    message("Outliers do not associate with any other traits. Try relaxing the search_threshold")
    return(output)
  }
  message("Found ", length(unique(out2$id.outcome)), " candidate traits associated with outliers at p < ", search_threshold)
  return(output)
}
