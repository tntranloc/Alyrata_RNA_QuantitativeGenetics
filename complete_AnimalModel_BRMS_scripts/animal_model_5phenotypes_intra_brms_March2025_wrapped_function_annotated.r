## ----------------------------------------------------------------------
## Simplified genetic variance‐component script for 5 phenotypes
## Author: Nhu L.T. Tran
## Date: June 2025
## Purpose: Run Bayesian animal model on 5 phenotypes locally,
##          then plot broad‐sense heritability (H2) for each.
## ----------------------------------------------------------------------

# 1) Setup ----------------------------------------------------------------
library(brms)
library(bayestestR)
library(ggplot2)


# Where your data live:
setwd("/Users/nhutran/Documents/PhD/epidom/runs/BRMS_AnimalModel/brms_phenotypes/brms_f1_all/intra")
#   - 'expr' is a data.frame (rows = individuals, cols = phenotypes)
#   - rownames(expr) are animal IDs
#   - 'dam_info' has columns: animal, dam
#   - Amat, Dmat are your additive & dominance relationship matrices
# Adjust paths as needed:
expr = read.csv("input_phenotypes.csv", row.names = 1, check.names = F)
dam_info = read.csv("dam_info.csv", header = T)
load("Amat.RData")  # loads Amat
load("Dmat.RData")  # loads Dmat

# For illustration, assume you already have:
# expr, dam_info, Amat, Dmat in your workspace.

# List your five phenotypes (must match colnames(expr))
pheno_list = c("leaf_thickness","rosette_area","flowering_time","number_of_flowers","leaf_serration")

# 2) Model settings -------------------------------------------------------
n_chains = 4
n_iter   = 20000   
n_warmup = 5000
thin     = 10

# 3) Loop + fit -----------------------------------------------------------
rownames(expr) = as.character(rownames(expr))
colnames(dam_info) = c("animal","dam")
dam_info$animal = as.character(dam_info$animal)
dam_info$dam  = as.character(dam_info$dam)


results = lapply(pheno_list, function(ph) {
  
  message("Fitting phenotype: ", ph)
  
  # Prepare data
  d = data.frame(
    y      = expr[[ph]],
    animal = rownames(expr)
  )
  d$dom = d$animal
  d$dam = factor(dam_info$dam[match(d$animal, dam_info$animal)])
  
  # Formula
  f = brmsformula(
    y ~ 1 +
      (1 | gr(animal, cov = Amat)) +
      (1 | gr(dom,    cov = Dmat)) +
      (1 | dam),
    family = skew_normal()
  )
  
  # Fit
  fit = brm(
    formula    = f,
    data       = d,
    data2      = list(Amat = Amat, Dmat = Dmat),
    chains     = n_chains, cores = n_chains,
    iter       = n_iter, warmup = n_warmup, thin = thin,
    control    = list(adapt_delta = 0.99),
    refresh    = 0
  )
  
  # Extract posterior SDs → variances
  vc = VarCorr(fit, summary = FALSE)
  Va = as.numeric(vc$animal$sd)^2
  Vd = as.numeric(vc$dom$sd)^2
  Vm = as.numeric(vc$dam$sd)^2
  Vr = as.numeric(vc$residual$sd)^2
  
  # Rhat
  RA = summary(fit)$random$animal$Rhat
  RD = summary(fit)$random$dom$Rhat
  RM = summary(fit)$random$dam$Rhat
  RR1 = summary(fit)$spec_pars$Rhat[1]
  RR2 = summary(fit)$spec_pars$Rhat[2]
  
  # Broad‐sense H2 per draw
  Vg = Va + Vd
  H2 = Vg / (Vg + Vm + Vr)
  Va_scaled = Va / Vg
  Vd_scaled = Vd / Vg
  
  # Summarize
  med_Va = median(Va)
  med_Vd = median(Vd)
  med_Vm = median(Vm)
  med_Vr = median(Vr)

  med_H2 = median(H2)
  med_Va_scaled = median(Va_scaled)
  med_Vd_scaled = median(Vd_scaled)
  ci  = hdi(H2, ci = 0.95)
  
  
  data.frame(
    phenotype      = ph,
    Va = med_Va,
    Vd = med_Vd,
    Vm = med_Vm,
    Vr = med_Vr,
    Va_scaled = med_Va_scaled,
    Vd_scaled = med_Vd_scaled,
    h2_med         = med_H2,
    h2_CI_lower    = ci[2],
    h2_CI_upper    = ci[3],
    Rhat_animal = RA, Rhat_dom = RD, Rhat_dam = RM,
    Rhat_resid_sigma = RR1, Rhat_resid_alpha = RR2
  )
})

res_df = do.call(rbind, results)
write.csv(res_df, "animal_model_results_5phenotypes_intra_Jun2025.csv", row.names = F)

# load the data back
res_df = read.csv("/Users/nhutran/Documents/PhD/epidom/runs/BRMS_AnimalModel/brms_phenotypes/brms_f1_all/intra/animal_model_results_5phenotypes_intra_Jun2025.csv")
head(res_df)

# 4) Plot -----------------------------------------------------------------
ggplot(res_df, aes(x = phenotype, y = h2_med)) +
  geom_col() +
  geom_errorbar(aes(ymin = h2_CI_lower, ymax = h2_CI_upper), width = 0.2) +
  labs(
    x    = "Phenotype",
    y    = "Broad‐sense H² (median ± 95% CI)",
    title= "Genetic Variance Proportions Across 5 Phenotypes"
  ) +
  theme_minimal(base_size = 14)



library(tidyr)
library(ggplot2)

# 1) Pivot to long form
df_long = res_df %>%
  pivot_longer(
    cols      = c("Va_scaled", "Vd_scaled"),
    names_to  = "component",
    values_to = "scaled_value"
  )

# 2) Plot grouped bars
df_long$component = factor(df_long$component, levels = c("Va_scaled", "Vd_scaled"))
df_long$phenotype = factor(df_long$phenotype, levels = c("leaf_serration","flowering_time","leaf_thickness","rosette_area","number_of_flowers"))

#or
df_long$component = factor(df_long$component, levels = c("Vd_scaled", "Va_scaled"))
df_long$phenotype = factor(df_long$phenotype, levels = c("number_of_flowers","rosette_area","leaf_thickness","flowering_time","leaf_serration"))


theme_set(
  theme_minimal(base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#dddddd", linetype = "dashed"),
      axis.line = element_line(color = "#444444", size = 0.8)
    )
)
ggplot(df_long, aes(x = phenotype, y = scaled_value, fill = component)) +
  #geom_col(position = position_dodge(width = 0.75), width = 0.6) +
  geom_col(width = 0.6, alpha= 0.8) +
  scale_fill_manual(
    values = c("Va_scaled" = "#82D08C", "Vd_scaled" = "#FFAE55"),
    labels = c("Additive", "Dominance"),
    name   = "Genetic Variance Component"
  ) +
  labs(
    x     = "Phenotype",
    y     = "Variance over Vg",
    title = "Additive vs. Dominance Proportions Across Phenotypes"
  ) +
  #theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )


#######################################

genetic_variance = read.table("/Users/nhutran/Documents/PhD/epidom/runs/BRMS_AnimalModel/combined_F1_intra_mean_median_animalmod_March2025.csv")  # Genetic variance data (genes x components)
colnames(genetic_variance)=c("Gene","Va_median", "Vd_median", "Vm_median", "Vr_median","Va_mean", "Vd_mean", "Vm_mean", "Vr_mean","CI_Va_low", "CI_Va_high", "CI_Vd_low", "CI_Vd_high", "CI_Vm_low", "CI_Vm_high", "CI_Vr_low", "CI_Vr_high", "RA", "RD", "RM", "RR1", "RR2")

genetic_variance$Va_scaled = genetic_variance$Va_median / (genetic_variance$Va_median + genetic_variance$Vd_median)
genetic_variance$Vd_scaled = genetic_variance$Vd_median / (genetic_variance$Va_median + genetic_variance$Vd_median)

library(ggplot2)

special_x = c(0.699,0.330,0.844,0.238,0.872)
special_x_vd = 1 - special_x

ggplot(genetic_variance, aes(x = Vd_scaled)) +
  geom_histogram(bins = 30, fill = "#FFAE55", color = "chocolate", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(Vd_scaled, na.rm = TRUE)),
             color = "black", linetype = "dashed", size = 0.8) +
  # add small red ticks at the bottom
  geom_rug(data = data.frame(x = special_x_vd), 
           aes(x = x), linewidth = 1.5,
           sides = "b",                   # bottom only
           color = "brown",
           length = unit(0.04, "npc")) +  # tick length ~2% of plot height
  #theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Scaled Genetic Variance (Vd_scaled)",
    x     = "Scaled Genetic Variance (Vd_scaled)",
    y     = "Frequency"
  ) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    axis.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


ggplot(genetic_variance, aes(x = Va_scaled)) +
  geom_histogram(bins = 30, fill = "#82D08C", color = "springgreen4", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(Vd_scaled, na.rm = TRUE)),
             color = "black", linetype = "dashed", size = 0.8) +
  # add small red ticks at the bottom
  geom_rug(data = data.frame(x = special_x), 
           aes(x = x), linewidth = 1.5,
           sides = "b",                   # bottom only
           color = "brown",
           length = unit(0.04, "npc")) +  # tick length ~2% of plot height
  #theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Scaled Genetic Variance (Va_scaled)",
    x     = "Scaled Genetic Variance (Va_scaled)",
    y     = "Frequency"
  ) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    axis.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

### Plot all components
# 0) get the ratio
res_df$Va_ratio = res_df$Va / (res_df$Va + res_df$Vd + res_df$Vm + res_df$Vr)
res_df$Vd_ratio = res_df$Vd / (res_df$Va + res_df$Vd + res_df$Vm + res_df$Vr)
res_df$Vm_ratio = res_df$Vm / (res_df$Va + res_df$Vd + res_df$Vm + res_df$Vr)
res_df$Vr_ratio = res_df$Vr / (res_df$Va + res_df$Vd + res_df$Vm + res_df$Vr)

# 1) Pivot to long form
df_long = res_df %>%
  pivot_longer(
    cols      = c("Va_ratio", "Vd_ratio", "Vm_ratio", "Vr_ratio"),
    names_to  = "component",
    values_to = "scaled_value"
  )
# 2) Plot grouped bars
df_long$component = factor(df_long$component, levels = c("Vr_ratio","Vm_ratio","Vd_ratio", "Va_ratio" ))
df_long$phenotype = factor(df_long$phenotype, levels = c("leaf_serration","flowering_time","leaf_thickness","rosette_area","number_of_flowers"))


theme_set(
  theme_minimal(base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#dddddd", linetype = "dashed"),
      axis.line = element_line(color = "#444444", size = 0.8)
    )
)


ggplot(df_long, aes(x = phenotype, y = scaled_value, fill = component)) +
  #geom_col(position = position_dodge(width = 0.75), width = 0.6) +
  geom_col(width = 0.6, alpha= 0.8) +
  scale_fill_manual(
    values = c("Va_ratio" = "#82D08C", "Vd_ratio" = "#FFAE55", "Vr_ratio" = "#C5C6C7", "Vm_ratio" = "#74D8E1"),
    name   = "Genetic Variance Component"
  ) +
  labs(
    x     = "Phenotype",
    y     = "Variance",
    title = "Additive vs. Dominance Proportions Across Phenotypes"
  ) +
  #theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(angle = 0, hjust = 0.5)
  )

