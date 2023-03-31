# DESCRIPTION: Check heterozygosity in the pruned and combined genotypes
# AUTHOR: Julia Romanowska
# DATE CREATED: 2021-11-09
# DATE MODIFIED: 2021-11-29

#' @param het data calculated from PLINK
#' @return ggplot: histogram of values with mean +/- 3*SD marked
check_het <- function(het){
  het_mean <- mean(het$F)
  het_sd <- sd(het$F)
  het_3sd <- c(het_mean - 3*het_sd, het_mean + 3*het_sd)
  het_median <- median(het$F)
    
  het %>%
    ggplot(aes(F)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = het_3sd, color = "red", lty = 2) +
    geom_vline(xintercept = het_mean, color = "red", lty = 1, size = 1.3) +
    geom_vline(xintercept = het_median, color = "black") +
    labs(title = "Heterozygosity histogram",
         subtitle = "with mean +/- 3*SD shown in red and median in black")
}

filter_het <- function(het){
  het_mean <- mean(het$F)
  het_sd <- sd(het$F)
  het_3sd <- c(het_mean - 3*het_sd, het_mean + 3*het_sd)
  
  # Filter ---
  # exclude persons with heterozygosity farther than 3SD from the sample mean
  het_ok <- het %>%
    filter(
      F >= het_3sd[1] &
      F <= het_3sd[2]
    )
  het_ok
  
  # Save ---
  out_file <- here("DATA", "cntrl_all_chr_QC_sample_OK.txt")
  write_delim(
    het_ok %>%
      select(FID, IID),
    file = out_file,
    delim = "\t"
  )
  return(out_file)
}
