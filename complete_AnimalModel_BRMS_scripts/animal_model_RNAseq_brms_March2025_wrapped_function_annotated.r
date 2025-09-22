## ----------------------------------------------------------------------
## Nhu L.T. Tran
## March 2025
## Purpose: Estimate genetic variance components (VA, VD, VM, VR)
##          from gene expression data using Bayesian animal models.
## 
## Model is fitted using `brms` with a skew-normal likelihood.
## 
## This script is designed to run on a cluster, where 200 genes are
## processed per job script using a script generator.
##
## Features:
## - Automatically skips genes already processed (via individual file check)
## - Saves per-gene output as well as batch-level summaries
## ----------------------------------------------------------------------


# Set temporary directory for R compilation
Sys.setenv(TMPDIR = "/scratch/ntran5/tmp/")

# Main function to run animal model across a batch of genes
animal_model_job_function = function() {

# Store results from each gene
out1 = list() 

for (iter in st:end) {
  gene = colnames(counts)[iter]
  output_gene_file = paste0(outdir, "/output_", gene, ".csv")
  
  # Skip if result already exists (to allow parallel job execution)
  if (individual_files && file.exists(output_gene_file)) {
    next  # Skip if already done
  }

  # Prepare input data for brms
  temp = data.frame(counts = as.numeric(counts[, iter]))
  temp$animal = as.character(rownames(counts))
  temp$dom = temp$animal
  temp$dam = factor(dam_info$dam[match(temp$animal, dam_info$animal)])
  
  # Define the Bayesian animal model formula
  myformula = brmsformula(counts ~ 1 +
                          (1 | gr(animal, cov = Amat)) +
                          (1 | gr(dom, cov = Dmat)) +
                          (1 | dam))
  # Fit the model using brms
  fit = brm(myformula,
            data = temp,
            family = skew_normal(),
            data2 = list(Amat = Amat, Dmat = Dmat),
            chains = 4, cores = 4, iter = 20000, warmup = 5000, thin = 10,
            threads = threading(4),control = list(adapt_delta = 0.99),
            refresh = 0)
  # Extract posterior draws of variance components and square to get variances
  Va_posterior = unlist(VarCorr(fit, summary = FALSE)$animal$sd)^2
  Vd_posterior = unlist(VarCorr(fit, summary = FALSE)$dom$sd)^2
  Vm_posterior = unlist(VarCorr(fit, summary = FALSE)$dam$sd)^2
  Vr_posterior = unlist(VarCorr(fit, summary = FALSE)$residual$sd)^2
  
  # Total phenotypic variance (per draw)
  Vp_posterior = Va_posterior + Vd_posterior + Vm_posterior + Vr_posterior

  # Calculate ratios (genetic variance components)
  ha = Va_posterior / Vp_posterior
  hd = Vd_posterior / Vp_posterior
  hm = Vm_posterior / Vp_posterior
  hr = Vr_posterior / Vp_posterior
  
  # Already normalized because Vp = Va + Vd + Vm + Vr
  #draws = data.frame(ha, hd, hm, hr)

  # (Re)normalize per draw just to be sure
  #draws = draws / rowSums(draws)

  # Posterior summaries
  median_ha = median(ha); median_hd = median(hd)
  median_hm = median(hm); median_hr = median(hr)
  
  mean_ha = mean(ha); mean_hd = mean(hd)
  mean_hm = mean(hm); mean_hr = mean(hr)
  # 95% Credible Interval
  ha_CI = bayestestR::hdi(ha, ci = 0.95); hd_CI = bayestestR::hdi(hd, ci = 0.95)
  hm_CI = bayestestR::hdi(hm, ci = 0.95); hr_CI = bayestestR::hdi(hr, ci = 0.95)
  
  # Remove to save memory
  rm(draws) 

  # Rhat diagnostics
  RA = summary(fit)$random$animal$Rhat
  RD = summary(fit)$random$dom$Rhat
  RM = summary(fit)$random$dam$Rhat
  RR1 = summary(fit)$spec_pars$Rhat[1]
  RR2 = summary(fit)$spec_pars$Rhat[2]

  # Save result row as a named data frame
  outp = data.frame(
    gene = gene,
    h2_add_med = median_ha, h2_dom_med = median_hd,
    h2_mat_med = median_hm,h2_resid_med = median_hr,
    h2_add_mean = mean_ha, h2_dom_mean = mean_hd,
    h2_mat_mean = mean_hm, h2_resid_mean = mean_hr,
    h2_add_CI_low = ha_CI[2], h2_add_CI_high = ha_CI[3],
    h2_dom_CI_low = hd_CI[2], h2_dom_CI_high = hd_CI[3],
    h2_mat_CI_low = hm_CI[2], h2_mat_CI_high = hm_CI[3],
    h2_resid_CI_low = hr_CI[2], h2_resid_CI_high = hr_CI[3],
    Rhat_animal = RA, Rhat_dom = RD, Rhat_dam = RM,
    Rhat_resid_sigma = RR1, Rhat_resid_alpha = RR2
  )
  # Store result in list
  out1[[iter]] = outp
  # Save per-gene output if specified
  if (individual_files) {
    write.table(outp, file = output_gene_file,
                sep = "\t", row.names = FALSE, col.names = TRUE)
  }
 
  # Memory cleanup
  rm(fit)
  gc()
}

# Combine and save batch output
if (length(out1) > 0) {
  out1_df = do.call(rbind, out1)
  output_file = paste0(outdir, "/output_batch_", st, "_", end, ".tsv") 
  write.table(out1_df, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  message("All genes in this batch were already processed. No output written.")
}


}

# Wrap and run the function
animal_model_job_function()
