#### METHOD 1 #####
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



#### METHOD 2 #####
## Required: Python 3.8, conda, and pixy 

## Installation
# Set up environment
mkdir ~/pixy_env
cd ~/pixy_env

# Create a virtual environment
python3 -m venv pixy_venv

# Activate the virtual environment
source pixy_venv/bin/activate

pip install pixy
pixy --version 

# Activate it anew by 
source ~/pixy_env/pixy_venv/bin/activate

# Then run simply by
pixy --stats pi fst dxy \
--vcf input.vcf.gz \
--populations popinfo.txt \
--bed_file genes.bed


  ## simply put in pop info 2 columns, sample name and population name
  ## if there is only one population, put one, no problem here

# output headers are	chromosome	window_pos_1	window_pos_2	avg_pi	no_sites	count_diffs	count_comparisons	count_missing

### DONE, good luck! ####
