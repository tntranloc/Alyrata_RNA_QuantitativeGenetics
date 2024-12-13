## This is to achieve pN/pS ratios through bootstrapping 

  # step 1 - divide the vcf regions into synonymous vs non-synonymous (snpeff)
  # step 2 -  for each category, calculate mean pi --> bootstrap n times
  # step 3 - for each bootstrap time, divide pi nonsynonymous / pi synonymous (aka piN/piS) - hence n values of piN/piS
  
# Use snpeff to categorise a vcf into synonymous and nonsynonymous regions
  # check snpeff.sh for more details

# On a good cluster, the following steps take 30min-1hr

# Required vcftools, bcftools, htslib

# Input files and parameters
NON_SYNON_VCF="nonsynonymous.vcf.gz"
SYNON_VCF="synonymous.vcf.gz"
POP_FILE="populations.txt" # column 1 is sample name, column 2 is population it belongs to. one population only is no problem. 2 columns are tab separated
BOOTSTRAPS=400
SNPS=10000
OUTFILE="bootstrap_pnps_ratios.txt"

# Initialize the output file
echo "Bootstrap_ID MeanPi_Ratio" > $OUTFILE

# Perform bootstrapping
cd /working/dir/

for i in $(seq 1 $BOOTSTRAPS); do
    echo "Bootstrap iteration $i"

    # Subsample nonsynonymous SNPs
    bcftools view -H $NON_SYNON_VCF | shuf -n $SNPS > nonsyn_sampled_snps.txt
    bcftools view -R nonsyn_sampled_snps.txt $NON_SYNON_VCF -Oz -o nonsyn_sampled.vcf.gz
    tabix -p vcf nonsyn_sampled.vcf.gz

    # Subsample synonymous SNPs
    bcftools view -H $SYNON_VCF | shuf -n $SNPS > syn_sampled_snps.txt
    bcftools view -R syn_sampled_snps.txt $SYNON_VCF -Oz -o syn_sampled.vcf.gz
    tabix -p vcf syn_sampled.vcf.gz

    # Calculate π for nonsynonymous variants using vcftools
    vcftools --gzvcf nonsyn_sampled.vcf.gz --site-pi --out nonsyn_pi_$i
    NONSYN_PI=$(awk '{if ($3 != "nan") sum+=$3; count++} END {if (count > 0) print sum/count; else print "NA"}' nonsyn_pi_>

    # Calculate π for synonymous variants using vcftools
    vcftools --gzvcf syn_sampled.vcf.gz --site-pi --out syn_pi_$i
    SYN_PI=$(awk '{if ($3 != "nan") sum+=$3; count++} END {if (count > 0) print sum/count; else print "NA"}' syn_pi_$i.sit>

    # Calculate the ratio
    if [[ "$SYN_PI" != "NA" && "$SYN_PI" != "0" ]]; then
        MEAN_PI_RATIO=$(echo "$NONSYN_PI / $SYN_PI" | bc -l)
    else
        MEAN_PI_RATIO="NA"
    fi
    
     # Save the result
    echo "$i $MEAN_PI_RATIO" >> $OUTFILE

    # Clean up intermediate files
    rm nonsyn_sampled_snps.txt syn_sampled_snps.txt nonsyn_sampled.vcf.gz* syn_sampled.vcf.gz*
    rm nonsyn_pi_$i.pi syn_pi_$i.pi
done

echo "Bootstrap completed. Results saved in $OUTFILE."






