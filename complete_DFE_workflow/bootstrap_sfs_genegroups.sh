#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --account=ag-demeaux
#SBATCH --error=/scratch/ntran5/errors/sfs.%J.err
#SBATCH --output=/scratch/ntran5/errors/sfs.%J.out

export PATH=/projects/ag-demeaux/ltntran/mytools/easySFS:$PATH

module load lang/Python/3.11.3-GCCcore-12.3.0
module load lang/Miniconda3/23.9.0-0
conda activate /projects/ag-demeaux/ltntran/myconda_p313

module load bio/BCFtools/1.18-GCC-12.3.0

INPUT=/genegroups/input_NoRepeats_NoHet.VarFiltration.rmINDELS.minQ30.minDP10.biallelic.maxmiss20.recode.SNPsOnly.vcf.gz
    # the vcf must be already marked for repeats, removed fixed heterozygosity, removed INDELS, 
    # minQ 30, minDP 10, kept for biallelic sites only, allowed for maximum 20% missingness
    # fixed heterozygosity means where all individuals were called as heterozygotes and no homozygotes were present 
POP=/pi/AlyrataParents_17PL_populations.tsv
    # simply a pop info file of 2 columns, individual name - pop name, f.e. 
    # Alyr1   PL
    # Alyr2   PL
GENE_BED=/genegroups/NT1_annotation_gene_group.bed 
    # 4 columns, chromosome - start - end - geneID, f.e. 
    #scaffold_1      10155851        10158686        AL1G36520
    #scaffold_1      10261716        10262735        AL1G36820
BOOTSTRAP_REPS=200  # Number of bootstrap replicates

#mkdir -p bootstrap_output
#mkdir -p bootstrap_output/FinalResults

cd /genegroups

for i in $(seq 1 $BOOTSTRAP_REPS); do
    echo "Running bootstrap replicate $i"

    # 1. Sample genes with replacement
    shuf -r $GENE_BED | head -n $(wc -l < $GENE_BED) > bootstrap_PL_highA_0fold/sampled_genes.bed

    # 2. Extract SNPs from resampled genes
    bcftools view -R bootstrap_output/sampled_genes.bed $INPUT -Oz -o bootstrap_output/boot_${i}.vcf.gz --threads 8
    #tabix -p vcf bootstrap_output/boot_${i}.vcf.gz

    # 3. Run EasySFS on the resampled VCF
    python3 easySFS.py -i bootstrap_output/boot_${i}.vcf.gz \
                       -p $POP \
                       --proj 28 \
                       --prefix boot_${i} \
                       -o bootstrap_output/boot_${i} \
                       -a -f
    # 4. Move the final SFS file from dadi/ to bootstrap_sfs/
    mv bootstrap_output/boot_${i}/dadi/PL.sfs bootstrap_output/FinalResults/boot_${i}.sfs

done

echo "Bootstrap complete!"
