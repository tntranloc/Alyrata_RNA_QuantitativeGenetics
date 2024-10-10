# Load Necessary Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("JASPAR2018", "TFBSTools", "Biostrings", "rtracklayer", "randomForest", "ggplot2"))

library(JASPAR2018)
library(TFBSTools)
library(Biostrings)
library(rtracklayer)
library(dplyr)
library(randomForest)
library(ggplot2)

# Step 1: Load Transcriptome Data and Gene Annotation
transcriptome_data <- read.csv("/path/to/transcriptome.csv")  # Load transcriptome data (CSV format)
gene_annotation <- read.csv("/path/to/gene_annotation.csv")  # Load gene annotation (CSV format)

# Step 2: Define Gene Groups
# Assuming you have a list or file indicating which genes belong to which group
group_1_genes <- c("gene1", "gene2", "gene3")  # Example gene IDs for group 1
group_2_genes <- c("gene4", "gene5", "gene6")  # Example gene IDs for group 2

group_1_annotation <- gene_annotation %>% filter(gene_id %in% group_1_genes)
group_2_annotation <- gene_annotation %>% filter(gene_id %in% group_2_genes)

# Step 3: Get Transcription Factor Profiles from JASPAR
opts <- list(species = "Arabidopsis thaliana", collection = "CORE")
pfm_list <- getMatrixSet(JASPAR2018, opts)  # Get TF PFMs for Arabidopsis thaliana

# Step 4: Convert PFMs to PWMs and Scan for TF Binding Sites
pwm_list <- lapply(pfm_list, toPWM)

# Function to Scan for TF Binding Sites in a Given Group
scan_tf_binding_sites <- function(pwm_list, transcriptome_data) {
  binding_sites_counts <- data.frame(gene_id = character(), tf_name = character(), count = integer())

  for (pwm in pwm_list) {
    tf_name <- name(pwm)
    counts <- sapply(transcriptome_data$sequence, function(seq) {
      sites <- searchSeq(pwm, DNAString(seq), strand = "*", min.score = "80%")
      length(sites)
    })
    binding_sites_counts <- rbind(
      binding_sites_counts,
      data.frame(gene_id = transcriptome_data$gene_id, tf_name = tf_name, count = counts)
    )
  }
  return(binding_sites_counts)
}

# Scan TF Binding Sites for Each Gene Group
binding_sites_counts_group_1 <- scan_tf_binding_sites(pwm_list, transcriptome_data %>% filter(gene_id %in% group_1_genes))
binding_sites_counts_group_2 <- scan_tf_binding_sites(pwm_list, transcriptome_data %>% filter(gene_id %in% group_2_genes))

# Step 5: Save TF Binding Sites Data for Each Group
write.csv(binding_sites_counts_group_1, file = "/path/to/binding_sites_counts_group_1.csv", row.names = FALSE)
write.csv(binding_sites_counts_group_2, file = "/path/to/binding_sites_counts_group_2.csv", row.names = FALSE)

# Step 6: Extract Gene Features for Each Group
extract_gene_features <- function(annotation, transcriptome_data) {
  gene_features <- annotation %>%
    mutate(
      transcript_length = transcriptome_data$length[match(gene_id, transcriptome_data$gene_id)],
      gene_length = annotation$gene_length,
      exon_count = annotation$exon_count
    ) %>%
    select(gene_id, transcript_length, gene_length, exon_count)
  return(gene_features)
}

gene_features_group_1 <- extract_gene_features(group_1_annotation, transcriptome_data)
gene_features_group_2 <- extract_gene_features(group_2_annotation, transcriptome_data)

# Save Gene Features Data for Each Group
write.csv(gene_features_group_1, "/path/to/gene_features_group_1.csv", row.names = FALSE)
write.csv(gene_features_group_2, "/path/to/gene_features_group_2.csv", row.names = FALSE)

# Step 7: Random Forest Analysis
# Load the feature data for each group
gene_features_group_1 <- read.csv("/path/to/gene_features_group_1.csv")
gene_features_group_2 <- read.csv("/path/to/gene_features_group_2.csv")

# Assuming VA and VD are columns in your data representing additive and dominance variance
# For this example, create dummy VA and VD columns
gene_features_group_1$VA <- runif(nrow(gene_features_group_1), min = 0, max = 1)
gene_features_group_1$VD <- runif(nrow(gene_features_group_1), min = 0, max = 1)

gene_features_group_2$VA <- runif(nrow(gene_features_group_2), min = 0, max = 1)
gene_features_group_2$VD <- runif(nrow(gene_features_group_2), min = 0, max = 1)

# Train Random Forest Models for VA and VD for Each Group
set.seed(123)  # For reproducibility

# Group 1 Random Forest for VA
rf_va_group_1 <- randomForest(VA ~ transcript_length + gene_length + exon_count,
                              data = gene_features_group_1, ntree = 500, importance = TRUE)
print(rf_va_group_1)

# Group 1 Random Forest for VD
rf_vd_group_1 <- randomForest(VD ~ transcript_length + gene_length + exon_count,
                              data = gene_features_group_1, ntree = 500, importance = TRUE)
print(rf_vd_group_1)

# Group 2 Random Forest for VA
rf_va_group_2 <- randomForest(VA ~ transcript_length + gene_length + exon_count,
                              data = gene_features_group_2, ntree = 500, importance = TRUE)
print(rf_va_group_2)

# Group 2 Random Forest for VD
rf_vd_group_2 <- randomForest(VD ~ transcript_length + gene_length + exon_count,
                              data = gene_features_group_2, ntree = 500, importance = TRUE)
print(rf_vd_group_2)

# Save Random Forest Results
saveRDS(rf_va_group_1, "/path/to/rf_va_group_1.rds")
saveRDS(rf_vd_group_1, "/path/to/rf_vd_group_1.rds")
saveRDS(rf_va_group_2, "/path/to/rf_va_group_2.rds")
saveRDS(rf_vd_group_2, "/path/to/rf_vd_group_2.rds")

# Step 8: Visualize Feature Importance from Random Forest
# Function to Plot Feature Importance
plot_feature_importance <- function(rf_model, title) {
  importance_data <- data.frame(Feature = rownames(importance(rf_model)),
                                Importance = importance(rf_model)[, 1])
  ggplot(importance_data, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = title, x = "Feature", y = "Importance")
}

# Plot Feature Importance for Each Model
plot_feature_importance(rf_va_group_1, "Feature Importance for VA - Group 1")
plot_feature_importance(rf_vd_group_1, "Feature Importance for VD - Group 1")
plot_feature_importance(rf_va_group_2, "Feature Importance for VA - Group 2")
plot_feature_importance(rf_vd_group_2, "Feature Importance for VD - Group 2")
