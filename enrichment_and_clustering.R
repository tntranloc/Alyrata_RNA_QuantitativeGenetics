## Goal: enrichment result of clustering of genes
## input: expression data, and gene universe that will be categorised into different groups

# Load required packages
library(topGO)
library(org.At.tair.db) # database of A.thaliana
library(stringr)
library(dplyr)
library(stats)
library(ggplot2)

# Set your data files
athaliana_orthologues <- "path/to/athaliana_orthologues.csv" # replace with your orthologue file
expression_data <- "path/to/expression_data.csv" # replace with your expression data

# Read in your orthologue data and transcriptome data
gene_universe <- read.csv(athaliana_orthologues, header = TRUE)
expression <- read.csv(expression_data, header = TRUE)

# Ensure your gene list for enrichment
ranked_genes <- gene_universe %>% 
  arrange(desc(VA)) # Ranking by VA (replace 'VA' by 'VNA' if needed)

# Prepare gene universe for topGO
gene_list <- factor(as.integer(ranked_genes$gene_id %in% ranked_genes$gene_id)) # Binary membership
names(gene_list) <- ranked_genes$gene_id

go_data <- new("topGOdata",
                description = "GO analysis",
                ontology = "BP", # Biological Process
                allGenes = gene_list,
                geneSel = function(p) p == 1,
                annot = annFUN.org, # Use organism annotation
                mapping = "org.At.tair.db")

# Fisher's exact test for enrichment
go_result <- runTest(go_data, 
                    algorithm = "classic", 
                    statistic = "fisher")

# Extract results
go_table <- GenTable(go_data, go_result, orderBy = "result", topNodes = 20)

# p-value threshold
significant_results <- go_table[go_table$P.value <= 0.0381, ]

# Save results
write.csv(significant_results, "GO_enrichment_results.csv")

# Compute pairwise Spearman correlation coefficients between gene expressions
spearman_matrix <- cor(expression, method = "spearman")

# Convert pairwise Spearman correlation to Euclidean distance matrix
D <- as.dist((1 - spearman_matrix) / 2)

# Hierarchical clustering
hclust_result <- hclust(D)

# Identify clusters by within-group sum of squares
wss <- sapply(1:300, function(k) {
  clusters <- cutree(hclust_result, k = k)
  sum(sapply(unique(clusters), function(cl) {
    cluster_data <- expression[clusters == cl, ]
    sum((cluster_data - rowMeans(cluster_data))^2)
  }))
})

# Plot within-group sum of squares to determine optimal number of clusters
plot(1:300, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")

# Prune to specific cluster numbers and analyze
cluster_23 <- cutree(hclust_result, k = 23)
cluster_50 <- cutree(hclust_result, k = 50)
cluster_100 <- cutree(hclust_result, k = 100)
cluster_200 <- cutree(hclust_result, k = 200)

# Save clustering results
write.csv(cluster_23, "cluster_23.csv")
write.csv(cluster_50, "cluster_50.csv")
write.csv(cluster_100, "cluster_100.csv")
write.csv(cluster_200, "cluster_200.csv")

# Correlate median VA, VNA, or transcript length with cluster size
cluster_data <- data.frame(cluster = cluster_23, VA = ranked_genes$VA)
cluster_medians <- cluster_data %>% 
  group_by(cluster) %>% 
  summarize(median_VA = median(VA), cluster_size = n())

correlation <- cor(cluster_medians$median_VA, cluster_medians$cluster_size)
print(correlation)
