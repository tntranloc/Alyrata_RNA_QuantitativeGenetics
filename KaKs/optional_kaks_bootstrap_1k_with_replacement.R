# Load the data
kaks = read.table("kaks_results_aggregated.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filtering
kaks = subset(kaks, kaks$P_value <= 0.05)
kaks = subset(kaks, kaks$Ka.Ks <= 1)
kaks = kaks[!is.na(kaks$`Ka.Ks`),]

kaks$Lyrata_ID = gsub("\\.t1$", "", kaks$Lyrata_ID)

# Load the gene groups of interets
newhighD = read.csv("high_hd_genes_intra.csv", header = T)
newhighA = read.csv("high_ha_genes_intra.csv", header = T)
newhighR = read.csv("high_hr_genes_intra.csv", header = T)
intermediate = read.csv("rest_genes_intra.csv", header = T)

# subset by gene groups
kaks_A = subset(kaks, kaks$Lyrata_ID %in% newhighA$Gene)
kaks_D = subset(kaks, kaks$Lyrata_ID %in% newhighD$Gene)
kaks_R = subset(kaks, kaks$Lyrata_ID %in% newhighR$Gene)
kaks_i = subset(kaks, kaks$Lyrata_ID %in% intermediate$Gene)

# randomly sample 4000 genes
allgenes = rbind(newhighA, newhighD, newhighR, intermediate)
set.seed(40)
samplegenes = allgenes[sample(nrow(allgenes), 4000),]
kaks_rand = subset(kaks, kaks$Lyrata_ID %in% samplegenes)

# rearranging df
kaks_A$Group = "highA"
kaks_D$Group = "highD"
kaks_R$Group = "highR"
kaks_i$Group = "iVg"
kaks_rand$Group = "random4k"

kaks_all = rbind(kaks_A, kaks_D, kaks_R, kaks_i, kaks_rand)
kaks_all$Group = as.factor(kaks_all$Group)
levels(kaks_all$Group)

# run bootstrap 
setwd("/.../kaks")
write.table(kaks_A, file = "kaks_table_highA_genes.txt", sep = "\t", row.names = FALSE, quote = F, col.names = T)
write.table(kaks_D, file = "kaks_table_highD_genes.txt", sep = "\t", row.names = FALSE, quote = F, col.names = T)
write.table(kaks_i, file = "kaks_table_iVg_genes.txt", sep = "\t", row.names = FALSE, quote = F, col.names = T)
write.table(kaks_R, file = "kaks_table_highR_genes.txt", sep = "\t", row.names = FALSE, quote = F, col.names = T)
write.table(kaks_rand, file = "kaks_table_random4k_genes.txt", sep = "\t", row.names = FALSE, quote = F, col.names = T)

# Set the directory containing your files
file_directory = "/.../kaks" 

# List all files in the directory
group_files = list.files(path = file_directory, pattern = "*_genes.txt", full.names = TRUE)

# Print the files to verify
print(group_files)

# Set bootstrap iterations
bootstrap_iterations = 1000

# Function to calculate mean pN/pS and bootstrap
process_group = function(file) {
  # Read the group file
  group_data = read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Ensure columns pN and pS are numeric
  group_data$Ka = as.numeric(group_data$Ka)
  group_data$Ks = as.numeric(group_data$Ks)
  group_data = group_data[!is.na(as.numeric(group_data$Ka)) & !is.na(as.numeric(group_data$Ks)),]
  # Replace zeros in pS with a small constant (e.g., 0.001) to avoid division by zero
  group_data$Ks[group_data$Ks == 0] <- 0.001
  print(head(group_data))  # Show the first few rows of the data
  print(str(group_data))   # Inspect the structure and data types of the columns
  # Calculate original mean pN, pS, and ratio
  mean_Ka = mean(group_data$Ka, na.rm = TRUE)
  mean_Ks = mean(group_data$Ks, na.rm = TRUE)
  original_ratio = mean_Ka / mean_Ks
  
  # Bootstrap to calculate mean pN, pS, and pN/pS ratio
  bootstrap_ratios = replicate(bootstrap_iterations, {
    sampled_data = group_data[sample(nrow(group_data), replace = TRUE), ]
    bootstrap_mean_Ka = mean(sampled_data$Ka, na.rm = TRUE)
    bootstrap_mean_Ks = mean(sampled_data$Ks, na.rm = TRUE)
    if (bootstrap_mean_Ks > 0) {
      bootstrap_mean_Ka / bootstrap_mean_Ks
    } else {
      NA
    }
  })
  
  # Return results as a list
  list(
    group_name = gsub("\\.txt$", "", basename(file)),
    mean_Ka = mean_Ka,
    mean_Ks = mean_Ks,
    original_ratio = original_ratio,
    bootstrap_ratios = bootstrap_ratios
  )
}

# Process each group and store results
bootstrap_results = lapply(group_files, process_group)


# Save bootstrap ratios for each group into a combined table
bootstrap_ratios_table = do.call(rbind, lapply(bootstrap_results, function(res) {
  data.frame(
    Group = res$group_name,
    Ka = res$mean_Ka,
    Ks = res$mean_Ks,
    Ratio = c(res$original_ratio, res$bootstrap_ratios),
    Type = c("Original", rep("Bootstrap", bootstrap_iterations))
  )
}))

# Extract the relevant part from the Group column
bootstrap_ratios_table$Group = sub(".*_table_(.*?)_genes.*", "\\1", bootstrap_ratios_table$Group)

bootstrap_ratios_file="kaks_ratios_bootstrap1k_allgroups_intra.csv"
#write.csv(bootstrap_ratios_table, bootstrap_ratios_file, row.names = FALSE)


# Load ggplot2 library
library(ggplot2)

# Load the bootstrap ratios table
#bootstrap_ratios_table = read.csv(bootstrap_ratios_file)
bootstrap_ratios_table$Group = as.factor(bootstrap_ratios_table$Group)

bootstrap_ratios_table$Group = factor(bootstrap_ratios_table$Group, levels = c("highA","highD", 
                                                                               "highR", "iVg",
                                                                               "random4k"))

library(dplyr)
bootstrap_ratios_table$Group = recode(bootstrap_ratios_table$Group,
                                      "highA" = "High VA",
                                      "highD" = "High VNA",
                                      "highR" = "High VR",
                                      "iVg" = "iVg",
                                      "random4k" = "Random",)


ggplot(bootstrap_ratios_table, aes(x = Group, y = Ratio, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "grey30", size = 0.5) +  # Boxplot
  geom_point(data = subset(bootstrap_ratios_table, Type == "Original"), 
             aes(x = Group, y = Ratio), 
             color = "brown", size = 2, alpha = 0.8) +  # Add Original points
  labs(
    title = "Bootstrap Ratios of KaKs by Group",
    x = "Gene Groups",
    y = "Ka/Ks Ratio"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#82D08C","#FFAE55","#C5C6C7","#909090", "#41424C"))


# Pairwise test
library(rstatix)
bootstrap_ratios_table %>% 
  group_by(Group) %>%
  get_summary_stats(Ratio)

res.kruskal = bootstrap_ratios_table %>% kruskal_test(Ratio ~ Group)
res.kruskal

pwc = bootstrap_ratios_table %>% 
  dunn_test(Ratio ~ Group, p.adjust.method = "bonferroni") 
pwc



########################################################## 
########## similar but also saving bootstraped Ks  #######
process_group = function(file) {
  # Read the group file
  group_data = read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Ensure pN and pS are numeric; drop any bad rows
  group_data$Ka = as.numeric(group_data$Ka)
  group_data$Ks = as.numeric(group_data$Ks)
  group_data = subset(group_data, !is.na(Ka) & !is.na(Ks))
  
  # Avoid divide‐by‐zero
  group_data$Ks[group_data$Ks == 0] <- 0.001
  
  # Calculate original means & ratio
  mean_Ka = mean(group_data$Ka, na.rm = TRUE)
  mean_Ks = mean(group_data$Ks, na.rm = TRUE)
  original_ratio = mean_Ka / mean_Ks

  # Bootstrap means
  # we'll capture both in a 2×B matrix
  boot_mat = replicate(bootstrap_iterations, {
    samp = group_data[sample(nrow(group_data), replace = TRUE), ]
    c(ka = mean(samp$Ka, na.rm = TRUE),
      ks = mean(samp$Ks, na.rm = TRUE))
  })
  # extract the two rows
  bootstrap_Ka      = boot_mat["ka", ]
  bootstrap_Ks      = boot_mat["ks", ]
  bootstrap_ratios  = bootstrap_Ka / bootstrap_Ks
  
  # Return all of them
  list(
    group_name       = gsub("\\.txt$", "", basename(file)),
    mean_Ka          = mean_Ka,
    mean_Ks          = mean_Ks,
    original_ratio   = original_ratio,
    bootstrap_Ka     = bootstrap_Ka,
    bootstrap_Ks     = bootstrap_Ks,
    bootstrap_ratios = bootstrap_ratios
  )
}


# Process each group and store results
bootstrap_results = lapply(group_files, process_group)

bootstrap_table = do.call(rbind, lapply(bootstrap_results, function(res) {
  B = length(res$bootstrap_ratios)
  data.frame(
    Group     = res$group_name,
    Type      = c("Original",   rep("Bootstrap", B)),
    Ka        = c(res$mean_Ka, res$bootstrap_Ka),
    Ks        = c(res$mean_Ks, res$bootstrap_Ks),
    Ratio     = c(res$original_ratio, res$bootstrap_ratios),
    stringsAsFactors = FALSE
  )
}))

# Inspect
head(bootstrap_table)

# Extract the relevant part from the Group column
bootstrap_table$Group = sub(".*_table_(.*?)_genes$", "\\1",bootstrap_table$Group)


bootstrap_ratios_file="ks_bootstrap1k_allgroups_intra.csv"
#write.csv(bootstrap_ratios_table, bootstrap_ratios_file, row.names = FALSE)

# Load the bootstrap ratios table
#bootstrap_ratios_table = read.csv(bootstrap_ratios_file)
bootstrap_table$Group = as.factor(bootstrap_table$Group)

bootstrap_table$Group = factor(bootstrap_table$Group, levels = c("highA","highD", 
                                                                               "highR", "iVg",
                                                                               "random4k"))

library(dplyr)
bootstrap_table$Group = recode(bootstrap_table$Group,
                                      "highA" = "High VA",
                                      "highD" = "High VNA",
                                      "highR" = "High VR",
                                      "iVg" = "iVg",
                                      "random4k" = "Random",)

ggplot(bootstrap_table, aes(x = Group, y = Ks, fill = Group)) +
  # standard boxplot
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    colour = "black",
    alpha = 0.8
  ) +
  # highlight the original means
  geom_point(
    data = subset(bootstrap_table, Type == "Original"),
    aes(x = Group, y = Ks),
    color = "salmon4",
    size = 2.5,
    shape = 18
  ) +
  # your custom fill palette
  scale_fill_manual(values = c("#82D08C","#FFAE55","#C5C6C7","#909090","#41424C")) +
  #theme_minimal(base_size = 14) +
  labs(
    title = "Bootstrap pS in Intra F1",
    x     = NULL,
    y     = "pS"
  ) +
  theme(
    plot.title     = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x    = element_text(size = 12),
    axis.text.y    = element_text(size = 12),
    legend.position= "none"
  )


# Pairwise test
library(rstatix)
bootstrap_ratios_table %>% 
  group_by(Group) %>%
  get_summary_stats(Ratio)

res.kruskal = bootstrap_ratios_table %>% kruskal_test(Ratio ~ Group)
res.kruskal

pwc = bootstrap_ratios_table %>% 
  dunn_test(Ratio ~ Group, p.adjust.method = "bonferroni") 
pwc
