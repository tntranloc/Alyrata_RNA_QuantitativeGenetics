library(MatrixEQTL)

# Load your data (replace with your file path)
load("genotype_matrix.RData")   # Should result in a matrix or data.frame
load("expression_matrix.RData")
load("snp_positions.RData")     # Data frame with columns: snp, chr, pos
load("gene_positions.RData")    # Data frame with columns: gene, chr, start, end

# Optional: Load covariates (e.g., PCs or batch effects)
# load("covariates.RData")

# Step 2. Prepare Matrix eQTL input objects
# Wrap genotype
snps = SlicedData$new()
snps$CreateFromMatrix(as.matrix(genotype_matrix))

# Wrap expression
gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(expression_matrix))

# If using covariates (can be empty)
cvrt = SlicedData$new()  # or: cvrt$CreateFromMatrix(as.matrix(covariates))

# Step 3. SNP and gene position formatting
# Format: SNP | CHR | POS
snp_pos = data.frame(snp = snp_positions$SNP_ID,
                     chr = snp_positions$CHR,
                     pos = snp_positions$POS)

# Format: Gene | CHR | START | END
gene_pos = data.frame(gene = gene_positions$Gene_ID,
                      chr = gene_positions$CHR,
                      left = gene_positions$START,
                      right = gene_positions$END)

# Step 4. Run Matrix eQTL
# Set model: linear or ANOVA (additive is most common)
useModel = modelLINEAR  # or modelANOVA

# Define output files
output_cis = "cis_eqtls.txt"
output_trans = "trans_eqtls.txt"

# Run!
eqtl_results = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_trans,
  pvOutputThreshold = 1e-5,               # trans threshold
  useModel = useModel,
  errorCovariance = numeric(),            # assume independence
  verbose = TRUE,
  output_file_name.cis = output_cis,
  pvOutputThreshold.cis = 1e-5,           # cis threshold
  cisDist = 5e6,                          # <5Mb = cis
  snpspos = snp_pos,
  genepos = gene_pos
)


# Output

#You will get two files:
#	cis_eqtls.txt → significant SNP–gene pairs where distance < 5 Mb
# trans_eqtls.txt → SNP–gene pairs where distance > 5 Mb

#Each file contains:
#	SNP ID, gene ID
#	beta, t-stat, p-value
#	(possibly adjusted p-value if you post-process)

# NOTE!
	#Ensure row/column names match exactly between matrices and position tables
	#Normalize expression data (e.g., log2(TPM + 1), or inverse-normal) for best results
	#SNPs and expression samples should be in the same order!


