% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OutliarAssociations.R
\name{outlierAssociations}
\alias{outlierAssociations}
\title{Outlier Associations}
\usage{
outlierAssociations(
  dat,
  id_list = "default",
  search_correction = "none",
  search_threshold = ifelse(search_correction == "none", 5e-08, 0.05),
  min_sample_size = 5000,
  min_SNP_number = 1e+05,
  use_proxies = FALSE
)
}
\arguments{
\item{dat}{Output from IdentifyOutliers function}

\item{id_list}{The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes().}

\item{search_correction}{Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").}

\item{search_threshold}{If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.}

\item{min_sample_size}{When selecting GWAS studies from MR-Base, state the minimum sample size acceptable. Default is 5000.}

\item{use_proxies}{Whether to use proxies when looking up associations. FALSE by default for speed.}

\item{min_SNP_numberWhen}{When selecting GWAS studies from MR-Base, state the minimum SNP number acceptable. Default is 5000.}
}
\value{
List
}
\description{
Function finds associations of outliers with candidate traits from MR-Base (https://www.mrbase.org/) and
adds these to an output list readable by the mrtryx package
}
