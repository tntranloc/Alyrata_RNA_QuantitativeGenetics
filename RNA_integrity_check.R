### Check the RNA integrity by checking the covered transcript region
### Fei 2020223
### Irina modified Feb-2022
### This script runs for Athaliana mapping files with Athaliana annotation file
### Read the depth files according to the chr
### First, get the depth file from uniquely mapping files by samtools in bash
# for file in *bam; do
#    samtools view -bq20 $file > temp.bam 
#    samtools depth temp.bam > $file.dp.txt
# done
# rm temp.bam

### R part for RNA integrity
setwd("/media/bene/Genomic_data/Abdul/short_submergence_data/short_submerg_exp_2023/transcriptome/mapping/map_both_with_nemgenome/")

### Load A.alpina gene annotation
load("./arabis_nem_chr_changed.RData")  ## We will focus here on gene
chr <- unique(annotations$V1)
gene.anno <- subset(annotations, V3 == "gene")
exon.anno <- subset(annotations, V3 == "exon")

files <- list.files(pattern = "depth.bam.txt")
cat("Files to process:", files, "\n")

### A function to calculate the proportion of sites in the gene that is covered by more than 1/2 of the mean coverage in a gene 
calculateProportionCovered <- function(gene.anno, depth) {
  start <- as.numeric(gene.anno[1])
  stop <- as.numeric(gene.anno[2])
  sites <- which(depth[, 1] >= start & depth[, 1] <= stop)  # Get all the sequence depth site for one gene
  if (length(sites) == 0) {
    proportionCovered <- NA
  } else {
    gene.cover <- depth[sites, 2]  # Depth for one gene
    mean.cover <- mean(gene.cover)  # Mean coverage
    covered <- sum(gene.cover >= mean.cover / 2)  # Coverage > 1/2 mean coverage
    proportionCovered <- covered / length(gene.cover)  # The proportion of covered to the length of the transcribed region
  }
  return(proportionCovered)  # Return as a number
}

### Parallel calculation in R
library(parallel)
start.time <- Sys.time()
no_cores <- detectCores() - 1  # Using the total number of CPUs - 1
cl <- makeCluster(no_cores)

proportionCoveredList <- vector("list", length(files))  # Pre-allocate list for all samples

for (k in seq_along(files)) {  # Loop through each file
  cat("Processing file:", files[k], "\n")
  proportionCoveredSum <- c()
  
  for (i in 1:8) {  # Loop through chromosomes 1 to 8
    cmd <- paste("sed 's/Chrom_//g'", files[k], "| awk '{if($1==", i, "&& $3 >4) print $2, $3}'", "> temp1.txt")
    cat("Running command:", cmd, "\n")
    system(cmd)

    # Check if "temp1.txt" is not empty
    if (file.size("temp1.txt") > 0) {
      cover.region <- read.table("temp1.txt")
      cat("File 'temp1.txt' read successfully with", nrow(cover.region), "rows\n")

      # Ensure chr[i] exists in mRNA.anno
      if (chr[i] %in% mRNA.anno$V1) {
        chr.anno <- subset(mRNA.anno, V1 == chr[i], select = c(V4, V5))
        proportionCovered <- parApply(cl, chr.anno, 1, calculateProportionCovered, depth = cover.region)

        cat("Chromosome:", i, "Processed\n")
        proportionCoveredSum <- c(proportionCoveredSum, proportionCovered)
      } else {
        cat("Chromosome:", i, "not found in mRNA.anno\n")
      }
    } else {
      warning("File 'temp1.txt' is empty or does not exist.")
    }
  }
  
  if (length(proportionCoveredSum) > 0) {
    proportionCoveredList[[k]] <- proportionCoveredSum
    cat("Processed data for file:", files[k], "with", length(proportionCoveredSum), "entries\n")
  } else {
    warning("No data processed for file:", files[k])
    proportionCoveredList[[k]] <- NA  # Ensure the list entry is not empty
  }
}

# Check if proportionCoveredList has entries and names vector matches its length
cat("Expected length:", length(files), "\n")
cat("Actual length:", length(proportionCoveredList), "\n")
if (length(proportionCoveredList) != length(files)) {
  stop("The length of proportionCoveredList does not match the number of files processed.")
}

# Assign names to the list using full filenames without extension
names(proportionCoveredList) <- sub(".depth.bam.txt$", "", files)

save(proportionCoveredList, file = "proportionCoveredList.RData")

stopCluster(cl)  # Stop the parallel
cat("Total processing time:", Sys.time() - start.time, "\n")

### Plot density for all samples
pdf("temp.pdf")
if (length(proportionCoveredList) > 0) {
  cover <- proportionCoveredList[[1]]
  n <- which(!is.na(cover) & cover < 1)
  plot(density(cover[n]), lwd = 2, col = 1, xlab = "% covered transcript", cex.lab = 1.3, main = "", ylim = c(0, 5))
  for (i in 2:length(proportionCoveredList)) {
    cover <- proportionCoveredList[[i]]
    n <- which(!is.na(cover) & cover < 1)
    lines(density(cover[n]), lwd = 2, col = i)
  }
}
dev.off()

# Discard samples that are not good
# Showing the median proportion of covered region, to spot the samples that show a problem.
cat("Median proportion covered per file:\n")
print(unlist(lapply(proportionCoveredList, median, na.rm = TRUE)))
