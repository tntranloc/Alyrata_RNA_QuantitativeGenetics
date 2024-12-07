## Required: vcftools and bedtools 
## Input vcf file, bed file of location of genes of interest

# Bed file structure, must follow!
# chr, start, end, strand, gene id
#chr1    50      150     + Gene1
#chr1    180     350    - Gene2
#chr1    400     600    + Gene3

### First calculate pi for every SNP position
vcftools --vcf input.vcf --site-pi --out snp_pi_values

### Add "end" column to file
awk 'NR>1 {print $1"\t"$2"\t"($2+1)"\t"$3}' snp_pi_values.sites.pi > snp_pi_values.bed
# If this doesn't work, put end position the same as start position
awk 'NR>1 {print $1"\t"$2"\t"($2)"\t"$3}' snp_pi_values.sites.pi > snp_pi_values.bed

### Then simply intersect the genes of interest and SNP positions
bedtools intersect -a genes.bed -b snp_pi_values.bed -wa -wb > snps_in_genes.bed

# output should look like
#chr1    50      150     Gene1   chr1    100    101    0.0012
#chr1    180     350     Gene2   chr1    200    201    0.0008
#chr1    180     350     Gene2   chr1    300    301    0.0015


### DONE, good luck! ####
