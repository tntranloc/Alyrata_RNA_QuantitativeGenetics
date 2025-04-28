## Prepare bimbam file for TWAS
  # first column: gene id
  # second column: a nucleotide (A f.e.)
  # third column: a corresponding nucleotide (T f.e.)
  # fourth column onwards: sample names
  # save the data with write.table(data, "dataname.txt", sep = ",", quote = F, row.names = F, col.names = F)

  
# Make kinship
gemma -g bimbamfile.txt -p trait.txt -gk 2 -o kinshipname 
  # gk 2 is for making kinship with bimbam input

# Make kinship with PLINK (highly recommended)
# VCF was filtered to keep INFO/DP >= 10 && F_MISSING < 0.2 && MAF >= 0.05, as well as biallelic SNPs only
# bcftools view -v snps -m2 -M2 -i 'INFO/DP >= 10 && F_MISSING < 0.2 && MAF >= 0.05
plink --vcf gwas_plinkdata --make-bed --out gwas_plinkdata
plink --bfile gwas_plinkdata_clean --make-rel square 
  # the plink.rel is the kinship I'm going to use

# Run TWAS
# if kinship was generated with gemma (which I experienced issues and ambiguity so I don't recommend)
gemma -g bimbamfile.txt -p traits.txt -k kinshipname.sXX.txt -lmm 4 -o twas_outputname

# if kinship was calculated with plink
gemma -g bimbamfile.txt -p traits.txt -k plink.rel -lmm 4 -o twas_outputname

## NOTE
# SNP QC OPTIONS
# -miss     [num]           specify missingness threshold (default 0.05)
# -maf      [num]           specify minor allele frequency threshold (default 0.01)
# -hwe      [num]           specify HWE test p value threshold (default 0; no test)
# -r2       [num]           specify r-squared threshold (default 0.9999)
# -notsnp                   minor allele frequency cutoff is not used

# When I do normal GWAS, I leave things as default
# However, when I run TWAS, since these are not "real" SNPs, to avoid them being filtered out, 
# I have to add -miss 1 -maf 0 -notsnp flag so nothing is filtered out
