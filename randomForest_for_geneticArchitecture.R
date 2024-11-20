#### Necessary Libraries ####
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("JASPAR2018", "TFBSTools", "Biostrings", "rtracklayer", "randomForest", "ggplot2"))

library(JASPAR2018)
library(TFBSTools)
library(Biostrings)
library(rtracklayer)
library(dplyr)
library(ranger)
library(ggplot2)

##### Find all used files as example in folder randomForest_for_geneticArchitecture #####


# Step 1: Load Transcriptome Data from FASTA file
sequences = readDNAStringSet("NT1_cds.fa") # my reference genome
gene_ids = gsub(".t1.v2.1", "", names(sequences)) # gene name was followed by these extra 
transcriptome_data = data.frame(gene_id = gene_ids, sequence = as.character(sequences), stringsAsFactors = FALSE) 

# Clean sequences by removing invalid characters and whitespace
transcriptome_data_clean = transcriptome_data %>%
  mutate(sequence = gsub("[^ATCGNatcgn]", "N", sequence)) %>%
  mutate(sequence = gsub("\\s+", "", sequence))


# Step 2: Define Geneset

genes = read.csv("/Users/nhutran/Documents/PhD/epidom/runs/genomicfactors_forRF/norm_counts_df_310_samples_intra_19k.csv", header = T)
rownames(genes) = genes$X
genes= genes[,-1]
#genes[1:10,1:10]
allgenes = colnames(genes)

# Check if all IDs in highA are present in transcriptome_data
missing_ids = setdiff(allgenes, transcriptome_data$gene_id)
if(length(missing_ids) > 0) print(missing_ids)

# Filter to keep only genes that are present in both datasets
common_genes = intersect(allgenes, transcriptome_data$gene_id)
genes_filtered = genes[, common_genes, drop = FALSE]

# Add missing genes with placeholder sequences if needed
for (gene in missing_ids) {
  transcriptome_data_clean = rbind(transcriptome_data_clean, data.frame(gene_id = gene, sequence = "NNNNNNNNNN"))
}

# Step 3: Load GFF and Extract Gene Annotations 
gff_file = "NT1_annotation_Novikova_edited_Leen.gff3" # my gff file
gff_data = import(gff_file, format = "gff")
gene_annotation = gff_data %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  select(gene_id = ID, seqnames, start, end) %>%
  mutate(gene_length = end - start + 1)

# Clean gene IDs
gene_annotation$gene_id = gsub(".v2.1", "", gene_annotation$gene_id)
head(gene_annotation) 
## GOT GENE LENGTH!


# Extract Gene Features function
extract_gene_features = function(annotation, transcriptome_data) {
  gene_features = annotation %>%
    mutate(transcript_length = nchar(transcriptome_data$sequence[match(gene_id, transcriptome_data$gene_id)])) %>%
    select(gene_id, transcript_length, gene_length)
  return(gene_features)
}

gene_features = extract_gene_features(gene_annotation, transcriptome_data_clean) 
head(gene_features)

# GOT GENE LENGTH and TRANSCRIPT LENGTH!!

# Count the number of CDS per gene, ensuring gene_id is character
cds_counts = gff_data %>%
  as.data.frame() %>%
  filter(type == "CDS") %>%
  mutate(gene_id = as.character(Parent)) %>%  # Convert Parent to character
  group_by(gene_id) %>%
  summarise(cds_count = n())

cds_counts$gene_id = gsub(".t1.v2.1", "", cds_counts$gene_id)  # Adjust gene_id format


# Merge CDS counts with gene annotations
gene_features = gene_features %>%
  left_join(cds_counts, by = "gene_id") %>%
  replace_na(list(cds_count = 0))

head(gene_features)

## NOTE that CDS in this gff is EXON

# GENE LENGTH, TRANSCRIPT LENGTH, and EXON COUNT


# Step 4: Get Transcription Factor Profiles from JASPAR
# TF motifs for Arabidopsis thaliana from the “CORE” collection, which includes well-characterized TFs.
opts = list(species = "Arabidopsis thaliana", collection = "CORE")
#  Fetches a set of PFMs matching your options, resulting in pfm_list, a list of frequency matrices representing the TF binding motifs
pfm_list = getMatrixSet(JASPAR2018, opts)
# Converts each PFM in pfm_list to a PWM, which is more suitable for scoring binding sites based on sequence
pwm_list = lapply(pfm_list, toPWM)

# Step 5: Define Function to Scan TF Binding Sites
# Define a modified function to accept a subset of PWMs and count binding sites to not jam the system
# upstream means how many kb upstream you scan for TF binding sites, here 1kb upstream

# this function counts how many sites a transcription factor binds to a gene
# out put has TFs as colnames, rownames as genes, and value is the number of binding sites 

scan_tf_binding_sites_upstream_count_batch = function(pwm_subset, gene_annotation, genome_seq, upstream_distance = 1000) {
  # Initialize an empty list to store counts of binding sites for each TF
  binding_sites_count = list()
  
  for (tf_name in names(pwm_subset)) {
    pwm = as.matrix(pwm_subset[[tf_name]]@profileMatrix)  # Extract PWM as a standard matrix
    tf_counts = integer(nrow(gene_annotation))  # Vector to store counts for each gene in this batch
    
    for (i in 1:nrow(gene_annotation)) {
      gene_id = gene_annotation$gene_id[i]
      start_pos = gene_annotation$start[i]
      seqname = gene_annotation$seqnames[i]
      
      if (!is.null(seqname) && seqname %in% names(genome_seq)) {
        # Define the upstream region
        upstream_start = max(1, start_pos - upstream_distance)
        upstream_end = start_pos - 1
        upstream_sequence = subseq(genome_seq[[seqname]], upstream_start, upstream_end)
        
        # Find all binding sites for the PWM in the upstream sequence
        sites = matchPWM(pwm, upstream_sequence, min.score = "80%")
        
        # Store the count of binding sites
        tf_counts[i] = length(sites)
      }
    }
    
    # Store counts for this TF in the list
    binding_sites_count[[tf_name]] = tf_counts
  }
  
  # Convert list of counts to a data frame and add gene IDs
  binding_sites_count_df = as.data.frame(binding_sites_count)
  binding_sites_count_df$gene_id = gene_annotation$gene_id
  binding_sites_count_df = binding_sites_count_df %>% select(gene_id, everything())
  
  return(binding_sites_count_df)
}


# Set batch size
batch_size = 60

# Split pwm_list into batches of items
pwm_batches = split(pwm_list, ceiling(seq_along(pwm_list) / batch_size))

# Confirm the number of batches
length(pwm_batches)

# Define a function that processes a single batch and saves the output
process_and_save_batch = function(batch_num) {
  cat("Processing batch", batch_num, "of", length(pwm_batches), "\n")
  
  # Select the current batch of PWMs
  pwm_batch = pwm_batches[[batch_num]]
  
  # Run the function on the current batch
  batch_result = scan_tf_binding_sites_upstream_count_batch(pwm_batch, gene_annotation, genome_seq)
  
  # Save the result immediately
  write.csv(batch_result, paste0("binding_sites_number_batch_", batch_num, ".csv"), row.names = FALSE)
  
  cat("Batch", batch_num, "saved successfully.\n")
}

# Manually run each batch as needed
process_and_save_batch(1)  # Run batch 1 and save 
## and etc until batch 8

# Combine everything!
combined = gene_features %>%
  left_join(binding_sites_number, by = "gene_id")

## save it!!
# without TF binding sites
write.csv(gene_features, "gene_features.csv", row.names = FALSE)
# with TF binding sites
write.csv(combined, "all_features.csv", row.names = FALSE)


## Prepare genetic variance data
genevar = read.csv("F1_ordered_hg_highest_gene_intra_19k.csv", header = T) # my df has colnames as genetic variance components and rownames as genes
dim(genevar) # 19642 10
dim(combined) # 29152 444

combined_sub = subset(combined, combined$gene_id %in% genevar$Gene)
genevar_sub = subset(genevar, genevar$Gene %in% combined_sub$gene_id)

## double check if it's 19641

genevar_sub = genevar_sub[match(combined_sub$gene_id, genevar_sub$Gene),]
identical(genevar_sub$Gene, combined_sub$gene_id) # TRUE

genevar_sub = genevar_sub[,-1]
colnames(genevar_sub)[1] = "gene_id"
identical(genevar_sub$gene_id, combined_sub$gene_id) # TRUE

combined_sub$VD = genevar_sub$hd

## save it! 
write.csv(combined_sub, "all_features_VD.csv", row.names = FALSE)


## Prepare inter data sets
df = read.csv("all_features_VD_intra.csv", header = T)

inter = read.csv("combined_F1_inter_19kgenes_animalmod.csv", header = F, sep = "\t")
inter = inter[,c(1,3)]
colnames(inter) = c("gene_id","VD_inter")
inter = subset(inter, inter$gene_id %in% df$gene_id)
inter = left_join(inter, df, by = "gene_id")
inter = inter[,-446]
#inter pop df
inter = inter[,-1]
#intra pop df
df = df[,-1]

### READY TO RUN RF

## Split train and test data
set.seed(42)
train_indices = sample(1:nrow(df), size = 0.7 * nrow(df))
train_data = df[train_indices,]
test_data = df[-train_indices,]

# run a simple basic rf model
names(df)
nfeatures= length(colnames(df)[1:443]) # everything except VD value
# train a default random forest model
intra_rf1 = ranger(
  VD ~ ., 
  data = train_data,
  mtry = floor(nfeatures / 3),
  respect.unordered.factors = "order",
  seed = 42
)
# get OOB RMSE
(default_rmse <- sqrt(ames_rf1$prediction.error))

## Find number of trees

# tuning grid
tuning_grid = expand.grid(
  trees = seq(10, 1000, by = 20),
  rmse  = NA
)
for(i in seq_len(nrow(tuning_grid))) {
  # Fit a random forest
  fit = ranger(
    formula = VD ~ ., 
    data = train_data, 
    num.trees = tuning_grid$trees[i],
    mtry = floor(nfeatures / 3),
    respect.unordered.factors = 'order',
    verbose = FALSE,
    seed = 42
  )
  
  # Extract OOB RMSE
  tuning_grid$rmse[i] <- sqrt(fit$prediction.error)
  
}
tuning_grid = na.omit(tuning_grid)
ggplot(tuning_grid, aes(trees, rmse)) +
  geom_line(size = 1) +
  ylab("OOB Error (RMSE)") +
  xlab("Number of trees")

# 200 trees would be good

## Find number of mtry

tuning_grid2 = expand.grid(
  trees = seq(10, 1000, by = 20),
  mtry  = floor(c(seq(2, 80, length.out = 5), 26)),
  rmse  = NA
)
for(i in seq_len(nrow(tuning_grid2))) {
  fit = ranger(
  formula    = VD ~ ., 
  data       = train_data, 
  num.trees  = tuning_grid2$trees[i],
  mtry       = tuning_grid2$mtry[i],
  respect.unordered.factors = 'order',
  verbose    = TRUE,
  seed       = 420
)
  
  tuning_grid2$rmse[i] = sqrt(fit$prediction.error)
  
}
tuning_grid2 = na.omit(tuning_grid2)
labels = tuning_grid2 %>%
  filter(trees == 500) %>%
  mutate(mtry = as.factor(mtry))

tuning_grid2 %>%
  mutate(mtry = as.factor(mtry)) %>%
  ggplot(aes(trees, rmse, color = mtry)) +
  geom_line(size = 1, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = labels, aes(trees, rmse, label = mtry), nudge_x = 50, show.legend = T) +
  ylab("OOB Error (RMSE)") +
  xlab("Number of trees")

# 20-25 trees is good


# create hyperparameter grid
hyper_grid <- expand.grid(
  mtry = floor(nfeatures * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  rmse = NA                                               
)
# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit <- ranger(
    formula         = VD ~ ., 
    data            = train_data, 
    num.trees       = nfeatures * 10,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = TRUE,
    seed            = 42,
    respect.unordered.factors = 'order',
  )
  # export OOB error 
  hyper_grid$rmse[i] = sqrt(fit$prediction.error)
}
# assess top 10 models
hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
  head(10)



## RUN different models 
rf_impurity = ranger(
  formula = VD ~ ., 
  data = train_data, 
  num.trees = 200,
  mtry = 21,
  min.node.size = 1,
  sample.fraction = .80,
  replace = FALSE,
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 42
)
rf_impurity_cor = ranger(
  formula = VD ~ ., 
  data = train_data, 
  num.trees = 200,
  mtry = 21,
  min.node.size = 1,
  sample.fraction = .80,
  replace = FALSE,
  importance = "impurity_corrected",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 42
)


rf_impurity_cor2 = ranger(
  formula = VD ~ ., 
  data = train_data, 
  num.trees = 500,
  mtry = 200,
  min.node.size = 1,
  sample.fraction = .80,
  replace = FALSE,
  importance = "impurity_corrected",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 42
)

rf_permutation = ranger(
  formula = VD ~ ., 
  data = train_data, 
  num.trees = 200,
  mtry = 21,
  min.node.size = 1,
  sample.fraction = .80,
  replace = FALSE,
  importance = "permutation",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 420
)

# View everything together
p1 = vip::vip(rf_impurity, num_features = 25, bar = FALSE)
p2 = vip::vip(rf_permutation, num_features = 25, bar = FALSE)
p3 = vip::vip(rf_impurity_cor, num_features = 25, bar = FALSE)
p4 = vip::vip(rf_impurity_cor2, num_features = 25, bar = FALSE)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)

# Test the model on test data
```{r}
test_predictions = predict(rf_impurity, data = test_data)$predictions
test_probablities = predict(rf_impurity, data = test_data, type = "response")$predictions
rmse = sqrt(mean((test_data$VD-test_predictions)^2))
baseline = sqrt(mean(mean((test_data$VD-test_predictions)^2)))
print(base)
print(rmse)
```

### GOOD LUCK! ####

