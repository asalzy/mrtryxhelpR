#' Adjust and analyse the tryx results
#'
#' Similar to tryx.analyse, but when there are multiple traits associated with a single variant then we use a LASSO-based multivariable approach
#'
#' @param tryxscan Output from \code{tryx.scan}
#' @param lasso Whether to shrink the estimates of each trait within SNP. Default=TRUE.
#' @param plot Whether to plot or not. Default is TRUE
#' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
#' @param proxies Look for proxies in the MVMR methods. Default = FALSE.
#'
#' @export
#' @return List
#' - adj_full: data frame of SNP adjustments for all candidate traits
#' - adj: The results from adj_full selected to adjust the exposure-outcome model
#' - Q: Heterogeneity stats
#' - estimates: Adjusted and unadjested exposure-outcome effects
#' - plot: Radial plot showing the comparison of different methods and the changes in SNP effects ater adjustment
tryxAnalyseMV <- function(tryxscan, lasso=TRUE, plot=TRUE, id_remove=NULL, proxies=FALSE)
{
  adj <- tryxAdjustmentMV(tryxscan, lasso=lasso, id_remove=id_remove, proxies=proxies)
  dat <- subset(adj$dat, mr_keep)
  dat$orig.ratio <- dat$orig.beta.outcome / dat$beta.exposure
  dat$orig.weights <- sqrt(dat$beta.exposure^2 / dat$orig.se.outcome^2)
  dat$orig.ratiow <- dat$orig.ratio * dat$orig.weights
  dat$ratio <- dat$beta.outcome / dat$beta.exposure
  dat$weights <- sqrt(dat$beta.exposure^2 / dat$se.outcome^2)
  dat$ratiow <- dat$ratio * dat$weights

  est_raw <- summary(lm(orig.ratiow ~ -1 + orig.weights, data=dat))
  est_adj <- summary(lm(ratiow ~ -1 + weights, data=dat))
  est_out1 <- summary(lm(orig.ratiow ~ -1 + orig.weights, data=subset(dat, !SNP %in% tryxscan$outliers)))
  est_out2 <- summary(lm(orig.ratiow ~ -1 + orig.weights, data=subset(dat, !SNP %in% adj$mvo$SNP)))

  estimates <- data.frame(
    est=c("Raw", "Outliers removed (all)", "Outliers removed (candidates)", "Outliers adjusted"),
    b=c(coefficients(est_raw)[1,1], coefficients(est_out2)[1,1], coefficients(est_out1)[1,1], coefficients(est_adj)[1,1]),
    se=c(coefficients(est_raw)[1,2], coefficients(est_out2)[1,2], coefficients(est_out1)[1,2], coefficients(est_adj)[1,2]),
    pval=c(coefficients(est_raw)[1,4], coefficients(est_out2)[1,4], coefficients(est_out1)[1,4], coefficients(est_adj)[1,4]),
    nsnp=c(nrow(dat), nrow(subset(dat, !SNP %in% tryxscan$outliers)), nrow(subset(dat, !SNP %in% adj$mvo$SNP)), nrow(dat)),
    Q = c(sum(dat$orig.qi), sum(subset(dat, !SNP %in% tryxscan$outliers)$qi), sum(subset(dat, !SNP %in% adj$mvo$SNP)$qi), sum(dat$qi)),
    int=0
  )
  estimates$Isq <- pmax(0, (estimates$Q - estimates$nsnp - 1) / estimates$Q)

  analysis <- list(
    estimates=estimates,
    mvo=adj$mvo,
    dat=adj$dat
  )

  if(plot)
  {
    dato <- subset(dat, SNP %in% tryxscan$outliers)
    datadj <- subset(dat, SNP %in% adj$mvo$SNP)
    datadj$x <- datadj$orig.weights
    datadj$xend <- datadj$weights
    datadj$y <- datadj$orig.ratiow
    datadj$yend <- datadj$ratiow

    mvog <- tidyr::separate(adj$mvo, exposure, sep="\\|", c("exposure", "temp", "id"))
    mvog$exposure <- gsub(" $", "", mvog$exposure)
    mvog <- group_by(mvog, SNP) %>% summarise(label = paste(exposure, collapse="\n"))
    datadj <- merge(datadj, mvog, by="SNP")

    labs <- rbind(
      tibble(label=dato$SNP, x=dato$orig.weights, y=dato$orig.ratiow, col="grey50"),
      tibble(label=datadj$label, x=datadj$weights, y=datadj$ratiow, col="grey100")
    )


    p <- ggplot(dat, aes(y=orig.ratiow, x=orig.weights)) +
      geom_abline(data=estimates, aes(slope=b, intercept=0, linetype=est)) +
      geom_label_repel(data=labs, aes(label=label, x=x, y=y, colour=col), size=2, segment.color = "grey50", show.legend = FALSE) +
      geom_point(data=dato, size=4) +
      # geom_label_repel(data=labs, aes(label=label, x=weights, y=ratiow), size=2, segment.color = "grey50") +
      geom_point(data=datadj, aes(x=weights, y=ratiow)) +
      geom_point(aes(colour=SNP %in% dato$SNP), show.legend=FALSE) +
      geom_segment(data=datadj, colour="grey50", aes(x=x, xend=xend, y=y, yend=yend), arrow = arrow(length = unit(0.01, "npc"))) +
      labs(colour=NULL, linetype="Estimate", x="w", y="beta * w")
    # xlim(c(0, max(dat$weights))) +
    # ylim(c(min(0, dat$ratiow, temp$ratiow), max(dat$ratiow, temp$ratiow)))
    analysis$plot <- p
  }
  return(analysis)

}


#' MR-TRYX Multivariate Adjustment
#'
#'  This function takes the output from the tryx::sig() function and runs multivariateMR to
#'  estimate adjusted candidate outlier - outcome associations.
#'
#' @param tryxscan Result from trys::sig() function
#' @param lasso Whether to shrink the estimates of each trait within SNP. Default=TRUE.
#' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
#' @param proxies Look for proxies in the MVMR methods. Default = FALSE.
#'
#' @return List
#' @export
#'
#' @examples
tryxAdjustmentMV <- function(tryxscan, lasso=TRUE, id_remove=NULL, proxies=FALSE)
{
  #Similar to adjustment function from tryx package.

  if(!any(tryxscan$search$sig))
  {
    return(tibble())
  }

  #Only keeping significant
  l <- list()
  sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)
  sige <- subset(tryxscan$candidate_exposure_mr, sig & !id.exposure %in% id_remove)
  sigo <- subset(tryxscan$candidate_outcome_mr, sig & !id.exposure %in% id_remove)

  #id of exposure and outcome of interest
  id.exposure <- tryxscan$dat$id.exposure[1]
  id.outcome <- tryxscan$dat$id.outcome[1]


  #outliers and candidate trait significantly associated with outcome
  sigo1 <- subset(sig, id.outcome %in% sigo$id.exposure) %>% dplyr::group_by(SNP) %>% dplyr::mutate(snpcount=n()) %>% dplyr::arrange(dplyr::desc(snpcount), SNP)
  snplist <- unique(sigo1$SNP)

  mvo <- list()

  #dat is the original exposure - outcome harmonised dataset
  #This sections adds heterogeneity statistics and whether the SNPs are outliers
  #betas, se and Q stats are separated into original (before adjustment - orig.beta) and new
  #following adjustment
  dat <- tryxscan$dat
  dat$qi <- cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
  dat$Q <- sum(dat$qi)
  dat$orig.beta.outcome <- dat$beta.outcome
  dat$orig.se.outcome <- dat$se.outcome
  dat$orig.qi <- dat$qi
  dat$outlier <- FALSE
  dat$outlier[dat$SNP %in% tryxscan$outliers] <- TRUE

  #Loop over all outlier SNPs where a candidate trait is significantly associated with the
  #outcome
  for(i in 1:length(snplist))
  {
    message("Estimating joint effects of the following trait(s) associated with ", snplist[i])
    temp <- subset(sigo1, SNP %in% snplist[i])
    candidates <- unique(temp$outcome)
    message(paste(candidates, collapse="\n"))
    if(lasso & length(candidates) == 1)
    {
      message("Only one candidate trait for SNP ", snplist[i], " so performing standard MVMR instead of LASSO")
    }

    #Here we are extracting the instrument for the main exposure and the candidate traits
    #(excluding those in id-remove) in preparation for the multivariate MR
    #The id-remove outcomes are removed as part of the temp variable.
    main_exposure_snps = dat %>% dplyr::pull(SNP)

    candidates_snps = tryxscan$candidate_instruments %>%
      dplyr::filter(id.exposure %in% unique(temp$id.outcome)) %>%
      dplyr::pull(SNP)

    #This is all the instruments for all the candidate traits and exposures. These
    #have not been clumped as this will be done by the mv_extract_local function
    unclumped_snp_list = c(unique(candidates_snps, main_exposure_snps))

    #This is the SNP - candidate trait associations
    mvexp_candidate = suppressMessages(TwoSampleMR::extract_outcome_data(snps = unclumped_snp_list,
                                           outcomes = unique(temp$id.outcome))) %>%
      TwoSampleMR::convert_outcome_to_exposure()

    #This is the SNP - main exposure associations
    #We are combinding the candidate SNP - Exposure associations and the Exposure
    #SNP - Exposure associations

    dat_exposure = dat %>% dplyr::select(SNP, dplyr::contains("exposure"))

    mvexp_exposure = tryxscan$candidate_exposure_dat %>%
      dplyr::filter(id.exposure%in% unique(temp$id.outcome)) %>%
      dplyr::select(SNP, contains("outcome")) %>%
      TwoSampleMR::convert_outcome_to_exposure()%>%
      dplyr::select(names(dat_exposure)) %>%
      rbind(dat_exposure)

    #Create a separate directory for all the multivariate files
    #This is removed at the end of the loop
    dir.create(file.path(getwd(), "MultivariateMR"))

    #Since the mv_extract_exposures_local function from TwoSampleMR can only read
    #in from a file and not an R data_frame, we have to create a new text file with
    #all the SNP instrument - candidate trait/exposure associations for each candidate
    #trait associated with that outlier

    export_MRTRYX = function (df) {
      readr::write_delim(df, paste("MultivariateMR//",unique(df$id.exposure), sep = ""))
      message("Exported associations for following phenotype: ", unique(df$Phenotype), " ID: ", unique(df$id.exposure))
    }

    mvexp_candidate %>%
      plyr::rename(c("effect_allele.exposure" = "effect_allele",
                     "other_allele.exposure" = "other_allele",
                     "beta.exposure" = "beta", "eaf.exposure" = "eaf",
                     "se.exposure" = "se", "pval.exposure" = "pval",
                     "exposure" = "Phenotype")) %>%
      split(.$Phenotype) %>%
      purrr::map(export_MRTRYX)

    mvexp_exposure %>%
      plyr::rename(c("effect_allele.exposure" = "effect_allele",
                     "other_allele.exposure" = "other_allele",
                     "beta.exposure" = "beta", "eaf.exposure" = "eaf",
                     "se.exposure" = "se", "pval.exposure" = "pval",
                     "exposure" = "Phenotype"))%>%
      split(.$Phenotype) %>%
      purrr::map(export_MRTRYX)


    #This is a character vector indicating the location of all the SNP - Exposure/Candidate Trait
    #Associations.
    path_id = list.files("./MultivariateMR", full.names = TRUE)
    mvexp = suppressMessages(TwoSampleMR::mv_extract_exposures_local(path_id,
                                                                     id_col = "id.exposure"))

    mvout = tryxscan$candidate_outcome_dat %>%
      dplyr::filter(SNP %in% mvexp$SNP) %>%
      dplyr::select(SNP, contains("outcome")) %>%
      dplyr::distinct()


    mvdat <- suppressMessages(TwoSampleMR::mv_harmonise_data(mvexp, mvout))

    #Deleting the multivariateMR folder
    unlink("./MultivariateMR", recursive = TRUE)
    if(lasso & length(candidates) > 1)
    {
      message("Performing shrinkage")
      b <- glmnet::cv.glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
      c <- coef(b, s = "lambda.min")
      keeplist <- unique(c(rownames(c)[!c[,1] == 0], tryxscan$dat$id.exposure[1]))
      message("After shrinkage keeping:")
      message(paste(sort(keeplist),mvdat$expname %>% dplyr::filter(id.exposure %in% keeplist) %>% dplyr::arrange(id.exposure) %>% dplyr::pull(exposure), collapse="\n"))

      mvexp2 <- subset(mvexp, id.exposure %in% keeplist)
      remsnp <- dplyr::group_by(mvexp2, SNP) %>% dplyr::summarise(mp = min(pval.exposure)) %>% dplyr::filter(mp > 5e-8) %$% as.character(SNP)
      mvexp2 <- subset(mvexp2, !SNP %in% remsnp)
      mvout2 <- subset(mvout, SNP %in% mvexp2$SNP)
      mvdat2 <- suppressMessages(TwoSampleMR::mv_harmonise_data(mvexp2, mvout2))
      mvo[[i]] <- TwoSampleMR::mv_multiple(mvdat2)$result
    } else {
      mvo[[i]] <- TwoSampleMR::mv_multiple(mvdat)$result
    }
    mvo[[i]] <- subset(mvo[[i]], !exposure %in% tryxscan$dat$exposure[1])
    temp2 <- with(temp, tibble(SNP=SNP, exposure=outcome, snpeff=beta.outcome, snpeff.se=se.outcome, snpeff.pval=pval.outcome))
    mvo[[i]] <- merge(mvo[[i]], temp2, by="exposure")
    boo <- with(subset(dat, SNP == snplist[i]),
                bootstrap_path(
                  beta.outcome,
                  se.outcome,
                  mvo[[i]]$snpeff,
                  mvo[[i]]$snpeff.se,
                  mvo[[i]]$b,
                  mvo[[i]]$se
                ))
    dat$beta.outcome[dat$SNP == snplist[i]] <- boo[1]
    dat$se.outcome[dat$SNP == snplist[i]] <- boo[2]
  }
  mvo <- bind_rows(mvo)

  dat$qi[dat$mr_keep] <- cochrans_q(dat$beta.outcome[dat$mr_keep] / dat$beta.exposure[dat$mr_keep], dat$se.outcome[dat$mr_keep] / abs(dat$beta.exposure[dat$mr_keep]))
  return(list(mvo=mvo, dat=dat))
}

cochrans_q <- function(b, se)
{
  xw <- sum(b / se^2) / sum(1/se^2)
  qi <- (1/se^2) * (b - xw)^2
  return(qi)
}

bootstrap_path <- function(gx, gx.se, gp, gp.se, px, px.se, nboot=1000)
{
  nalt <- length(gp)
  altpath <- tibble(
    p = rnorm(nboot * nalt, gp, gp.se) * rnorm(nboot * nalt, px, px.se),
    b = rep(1:nboot, each=nalt)
  )
  altpath <- group_by(altpath, b) %>%
    summarise(p = sum(p))
  res <- rnorm(nboot, gx, gx.se) - altpath$p
  pe <- gx - sum(gp * px)
  return(c(pe, sd(res)))
}

