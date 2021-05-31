#' Import Plink 2 Associations
#'
#' This function takes the output from the plink 2 -glm output and returns a
#' data frame readable by the TwoSampleMR package. Users can choose to include
#' .freq files with details concerning allele frequencies and these will be
#' included in the output data frame which aids in the process of harmonisation
#' when using the harmonise_data() function from the TwoSampleMR package.
#'
#' @param association_file output from plink2 -glm function
#' @param variable_type choose from c("exposure", "outcome")
#' @param eaf_file .afreq output from plink2 detailing allele frequencies (optional)
#' @param phenotype string detailing variable of interest (optional)
#'
#' @return data frame
#' @export
#'
#' @examples
importPlink2 = function(association_file, variable_type, eaf_file, phenotype) {

  stopifnot(is.character(association_file),
            is.character(variable_type))

  #Reading plink assoc. file and selecting columns relevant for TwoSampleMR
  #format_data() import. If a A1 frequency column is present, this is imported.
  #
  #Double check to make sure the effect allele and other allele are alligned properly
  assoc_file = read.delim(association_file) %>%
    dplyr::select(chr = contains("CHROM"), pos = POS, SNP = ID, effect_allele = A1, other_allele = REF, ALT,
                  samplesize = OBS_CT, beta = BETA, se = SE, pval = P, eaf = contains("A1_FREQ")) %>%
    dplyr::mutate(other_allele = ifelse(other_allele == effect_allele, ALT, other_allele))

  #If a separate eaf file is called, this is imported (only if a frequency column was not
  #found in the assoc file)
  if (!missing(eaf_file)){

    #Case when eaf has already been imported from .assoc file, ignore eaf file
    if ("eaf" %in% colnames(assoc_file)) {
      message("EAF found in association file, ignorning EAF file")
      break
    }

    eafrequency_file = read.delim(eaf_file)
    assoc_file = dplyr::left_join(assoc_file, eafrequency_file, by = c("SNP" = "ID")) %>%
      dplyr::rename(eaf = ALT_FREQS)
  }

  #If user doesn't provide any method to include EAF, message is shown
  if (!"eaf" %in% colnames(assoc_file)) {

  message("No Effect Allele Frequency collumn in association file or no separate
          EAF file included. Including EAF aids during harmonisation.")
  }

  #If phenotype value is provided, this is included in the formatted data

  if(!missing(phenotype)) {

    assoc_file$Phenotype = phenotype

  }

  formatted_data = TwoSampleMR::format_data(dat = assoc_file,
                                            type = variable_type)
  return(formatted_data)
}
