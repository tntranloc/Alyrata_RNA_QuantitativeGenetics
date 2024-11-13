### REQUIRED ###
# gatk 4.6.1.0
# Python/3.9.6
# Miniconda3/23.9.0
# GCCcore-11.2.0 
# Java 17.0.6
# picard.jar

#### To be ran on cluster ####
### last updated 12.11.2024 ###

#REF=reference_genome.fasta
#OUTDIR=/output/directory
#INPUTDIR=/input/directory
#SCRIPTDIR=/job/scripts/directory
#INPUT=input.vcf.gz
#OUTPUT=output.vcf.gz

### STEP 0: add headers for gatk ####  
##check your bam files to see if it's needed

# Loop through each BAM file in both input directories and create a separate script for each sample
for BAM in "$INPUTDIR"/*.bam; do
    #SAMPLE_NAME=$(basename "$BAM" .bam)
    # Extract only the base name without the suffix, f.e., for me, every sample name is followed by "Aligned.sortedByCoord.out.bam"
    SAMPLE_NAME=$(basename "$BAM" Aligned.sortedByCoord.out.bam)
    SCRIPT_FILE="$SCRIPT_DIR/${SAMPLE_NAME}_addheaders.sh"
    # Write the SLURM headers and GATK command to the script file
    cat <<EOT > "$SCRIPT_FILE"
java -Xmx4g -jar picard.jar AddOrReplaceReadGroups \
I="$BAM" \
O="$OUTDIR/${SAMPLE_NAME}.dedup_sort_withRG.bam" \
RGID="${SAMPLE_NAME}" \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit \
RGSM="${SAMPLE_NAME}" \
CREATE_INDEX=true

EOT
    # Make the script executable
    chmod +x "$SCRIPT_FILE"
done

### Excute all the script in the directory ###
for script in /job/scripts/directory/*.sh; do
sbatch $script
done


### STEP 1: Haplotype CAller 
# Loop through each BAM file in both input directories and create a separate script for each bam file
for BAM in "$INPUTDIR"/*.bam; do
    #SAMPLE_NAME=$(basename "$BAM" .bam)
    # Extract only the base name without the suffix ".dedup_sort_withRG.bam"
    SAMPLE_NAME=$(basename "$BAM" .dedup_sort_withRG.bam)
    SCRIPT_FILE="$SCRIPT_DIR/${SAMPLE_NAME}_haplotypecaller.sh"
    # Write the SLURM headers and GATK command to the script file
    cat <<EOT > "$SCRIPT_FILE"
  gatk --java-options "-XX:+UseSerialGC" HaplotypeCaller \\
    -R "$REF" \\
    -I "$BAM" \\
    -O "$OUTDIR/${SAMPLE_NAME}.g.vcf.gz" \\
    --min-pruning 1 \\
    --dont-use-soft-clipped-bases \\
    --disable-read-filter MappingQualityAvailableReadFilter \\ ## must have for RNA data
    -ERC GVCF
  EOT    
    # Make the scripts executable
  chmod +x "$SCRIPT_FILE"
done


### STEP 2: Combine all the GVCFs created in Step 1
gatk CombineGVCFs \
    -R "$REF" \
    $(for f in ${INPUTDIR}/*.g.vcf.gz; do echo "-V $f"; done) \
    -O ${OUTDIR}/combined.g.vcf.gz


### STEP 3: Calling the SNPs

gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R $REF \
    -V combined.g.vcf.gz \
    -O ${OUTDIR}/genotyped.vcf.gz


### STEP 4: Filter based on the needs ###
## To extract SNPs only (f.e. no INDELs)
gatk --java-options "-Xmx4g" SelectVariants \
    -R $REF \
    -V $INPUT \
    --select-type-to-include SNP \
    -O ${OUTDIR}/genotyped_311_chr1to8_renamed_SNPsOnly.vcf.gz 

## To filter max missingness, min and max depth 
bcftools view -i 'F_MISSING < 0.1 && FORMAT/DP > 10 && FORMAT/DP < 100 && MAF >= 0.05' $INPUT -o ${OUTDIR}/$OUTPUT -Oz

## To have statistics of the vcf file
bcftools stats input.vcf.gz > stats.txt

## To extract INFO DP for filtering
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' input.vcf.gz > info_dp_values.txt
# the in R:
info_dp = read.table("info_dp_values.txt", header = F, sep = "\t") # column 3 is DP value
ggplot(info_dp, aes(x = V3)) +
  geom_histogram(binwidth = 10, fill = "skyblue", color = "black") +
  labs(title = "Histogram of INFO DP (Limited x-axis)", x = "INFO DP", y = "Frequency") +
  xlim(0, 1000) +
  theme_minimal()

### BONUS ###
## To rename samples in vcf file
bcftools reheader -s ${INPUTDIR}/new_sample_names.txt -o ${OUTDIR}/$OUTPUT $INPUT 
tabix -p vcf ${OUTDIR}/$OUTPUT

## To subset chromosomes/scaffolds 
# for example, subset scaffold 1 to 8
bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8 $INPUT -Oz -o ${OUTDIR}/$OUTPUT
tabix -p vcf ${OUTDIR}/$OUTPUT

### PLEASE NOTE: don't demand too much memory, it's not necessary and will slow you down!!
## all scripts here used 4 cpus
## simple things like statistics check, add headers, filtering vcf using bcftools - ask for 4gb (max 8gb)
# haplotype caller: 8gb max
# combine vcf, genotype vcf, select variants: 24gb max


### on a good cluster, a whole process should take approx. 5 days, 1-2 days each step ###
### on a slow cluster, expect a month ##


    




