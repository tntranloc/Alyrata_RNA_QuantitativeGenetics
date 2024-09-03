## To be run on a cluster ##
############################

## STEP 1
########################## quality check ################################

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc1.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc1.%J.out

module load fastqc/0.11.9

cd /scratch/ltntran/test_transcriptome_workflow/raw_data

find ./ -type f -name "*.fastq.gz" | while read -r file; do
     fastqc "${file}" -o /scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc_initial -t 4
done


## STEP2
##########################  trimming with Trimmomatic ################################

# NOTE!!! make sure the adapter fa is in trimming directory or wherever you are

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/raw_data/trim.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/raw_data/trim.%J.out

module load trimmomatic/0.39
module load openjdk/1.8.0_60

RAWDIR=/scratch/ltntran/test_transcriptome_workflow/raw_data/
TRIMDIR=/scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/

export PATH=/projects/ag-demeaux/ltntran/mytools/:$PATH

cd "${TRIMDIR}"

find "${RAWDIR}" -name "*_R1_001.fastq.gz" | while read -r R1; do
    R2=$(echo $R1 | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
    base=$(basename "$R1" _R1_001.fastq.gz)
    java -jar /projects/ag-demeaux/ltntran/mytools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 -phred 33 "$R1" "$R2" \
"$base.trimmed.R1P.fastq.gz" "$base.trimmed.R1U.fastq.gz" "$base.trimmed.R2P.fastq.gz" "$base.trimmed.R2U.fastq.gz" \
ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 
done


########################## index reference genome ###########################################
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /ref_gen/ --genomeFastaFiles Alyrata.fasta --genomeSAindexNbases 8

########################## convert annotation gff file to gtf ###############################
gffread annotation.gff -T -o annotation.gtf


########################## quality check for trimmed product ################################

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc2.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc2.%J.out

module load fastqc/0.11.9

cd /scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/

find ./ -type f -name "*.fastq.gz" | while read -r file; do
     fastqc "${file}" -o /scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/fastqc_trimmed/
done


########################## mapping with STAR ################################

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/raw_data/map.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/raw_data/map.%J.out

module load star/2.7.8a
module load samtools/1.13

MAPDIR=/scratch/ltntran/test_transcriptome_workflow/mapping/
TRIMDIR=/scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/
REFDIR=/scratch/ltntran/test_transcriptome_workflow/ref_gen2

export PATH=/projects/ag-demeaux/ltntran/mytools/:$PATH

cd "${MAPDIR}"

find "${TRIMDIR}" -name "*.trimmed.R1P.fastq.gz" | while read -r R1; do
    R2=$(echo $R1 | sed 's/.trimmed.R1P.fastq.gz/.trimmed.R2P.fastq.gz/')
    base=$(basename "$R1" .trimmed.R1P.fastq.gz)
    STAR --runThreadN 2 --genomeDir "${REFDIR}" --sjdbGTFfile "${REFDIR}/Alyrata.gtf" --readFilesIn "$R1" "$R2" --readFilesCommand zcat --outFileNamePrefix "${base}" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
done





####################################################################################################################################
####################################################################################################################################
###################### running things in parallel by creating and submitting simultaneous bash scripts ##########################

###### STEP1 fastqc ######

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc_make_sh.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc_make_sh.%J.out

module load fastqc/0.11.9

cd /scratch/ltntran/test_transcriptome_workflow/raw_data

find ./ -type f -name "*.fastq.gz" | while read -r file; do
echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc1.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc1.%J.out
module use fastqc/0.11.9
cd /scratch/ltntran/test_transcriptome_workflow/raw_data/" > $file.submit.sh
echo "fastqc "${file}" -t 12 -o /scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc_initial/" >> $file.submit.sh
done

###### STEP2 trimming ######

#make sure to have the adapter fa file in trimming directory or wherever you are

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/trim.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/trim.%J.out

module load trimmomatic/0.39
module load openjdk/1.8.0_60

RAWDIR=/scratch/ltntran/test_transcriptome_workflow/raw_data/
TRIMDIR=/scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/
JOBDIR=/scratch/ltntran/test_transcriptome_workflow/jobs_parallel/trimming

export PATH=/projects/ag-demeaux/ltntran/mytools/:$PATH

find "${RAWDIR}" -name "*_1*.fq.gz" | while read -r R1; do
    R2=$(echo $R1 | sed 's/_1.fq.gz/_2.fq.gz/')
    base=$(basename "$R1" _1.fq.gz)
    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=480:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/trim.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/trim.%J.out

module load trimmomatic/0.39
module load openjdk/1.8.0_60
cd "${TRIMDIR}"
" > "${JOBDIR}/${base}.submit.sh"

echo "
java -jar /projects/ag-demeaux/ltntran/mytools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 -phred33 "$R1" "$R2" \
"$base.trimmed.1P.fq.gz" "$base.trimmed.1U.fq.gz" "$base.trimmed.2P.fq.gz" "$base.trimmed.2U.fq.gz" \
ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 
" >> "${JOBDIR}/${base}.submit.sh"
done





###### STEP2.5 optional: checking QC after trimming ######

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/fastqc_make_sh.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/fastqc_make_sh.%J.out

module load fastqc/0.11.9

cd /scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/

find /scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/ -type f -name "*.fastq.gz" | while read -r file; do
echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc2.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/raw_data/fastqc2.%J.out
module use fastqc/0.11.9
cd /scratch/ltntran/test_transcriptome_workflow/raw_data/trimming" > $file.submit.sh
echo "fastqc "${file}" -o /scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/fastqc_trimmed/" >> $file.submit.sh
done

###### STEP3 mapping ######


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/map_make_sh.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/map_make_sh.%J.out

module load star/2.7.8a 
module load samtools/1.13

MAPDIR=/scratch/ltntran/test_transcriptome_workflow/mapping/
TRIMDIR=/scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/
REFDIR=/scratch/ltntran/test_transcriptome_workflow/ref_gen
JOBDIR=/scratch/ltntran/test_transcriptome_workflow/jobs_parallel/mapping

export PATH=/projects/ag-demeaux/ltntran/mytools/:$PATH

cd "${MAPDIR}"

find "${TRIMDIR}" -name "*.trimmed.1P.fq.gz" | while read -r R1; do
    R2=$(echo $R1 | sed 's/.trimmed.1P.fq.gz/.trimmed.2P.fq.gz/')
    base=$(basename "$R1" .trimmed.1P.fq.gz)
    echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=480:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --account=UniKoeln
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntran5@uni-koeln.de
#SBATCH --error=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/map.%J.err
#SBATCH --output=/scratch/ltntran/test_transcriptome_workflow/error_out_stuff/map.%J.out
module load star/2.7.8a
module load samtools/1.13
MAPDIR=/scratch/ltntran/test_transcriptome_workflow/mapping/
TRIMDIR=/scratch/ltntran/test_transcriptome_workflow/raw_data/trimming/
REFDIR=/scratch/ltntran/test_transcriptome_workflow/ref_gen
cd "${MAPDIR}"
' > ${JOBDIR}/${base}.submit.sh
    echo "STAR --runThreadN 12 --genomeDir "${REFDIR}" --sjdbGTFfile "${REFDIR}/NT1_annotation_novikova_edited.gtf" --readFilesIn "$R1" "$R2" --readFilesCommand zcat --outFileNamePrefix "${base}" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts" >> ${JOBDIR}/${base}.submit.sh
    done



##to run all bash scripts in the folder (because I want to check if created scripts look okay, otherwise you can already submit on one go
#create a dir to move the scripts to, for organising sake
#make the file submit_all.sh

#!/bin/bash

cd  /where/the/sh_scripts/are/

for script in *.sh; do
    sbatch "$script"
done

#to run the script
chmod a+x submit_all.sh
bash submit_all.sh


## at this point according to my pipeline above I would have 320 tables for 320 samples
#STAR outputs read counts per gene into PREFIXReadsPerGene.out.tab file with 4 columns 
    #which correspond to different strandedness options: 
        #column 1: gene ID 
        #column 2: counts for unstranded RNA-seq column 
        #3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes) 
        #column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)


##### bring this to normalisation using DESeq2 and it's ready

##### perform PCA 

################################ GOOD LUCK ################################
