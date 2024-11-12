### REQUIRED
# gatk 4.6.1.0
# Python/3.9.6
# Miniconda3/23.9.0
# GCCcore-11.2.0 
# Java 17.0.6

#REF=reference_genome.fasta
#OUTDIR=/output/directory
#INPUTDIR=/input/directory

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

# Excute all the script in the directory
for script in /script/folder/*.sh; do
sbatch $script
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


### STEP 4: Filter based on the needs




