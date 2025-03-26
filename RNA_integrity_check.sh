## Required: Qualimap v2.3
## Optional: MultiQC

## Usage is very simple

## Check the details here : http://qualimap.conesalab.org/doc_html/command_line.html

qualimap rnaseq -bam sample_name.bam -gtf your_annotation.gtf -outdir /path/to/outdir/sample_name


###### On cluster ######
## Use this script to generate sh file for each sample, so that all jobs can be executed simultaneously 

############ Beginning of the "job generator" script  ################

#!/bin/bash

# Directories
BAMDIR="/path/to/bam_files"  # BAM file directory
JOBDIR="/path/to/save/job_scripts/"  # Output directory for .sh job files
REPORT_DIR="/path/to/save/qualimap/output"  # Directory for reports

GTF="/path/to/your_annotation.gtf"

#mkdir -p "$JOBDIR" "$FILTERED_DIR"

# Loop through BAM files and generate job scripts
find "$BAMDIR" -name "*Aligned.sortedByCoord.out.bam" | while read -r BAM; do
    base=$(basename "$BAM" Aligned.sortedByCoord.out.bam)
    
    JOB_FILE="${JOBDIR}/${base}_qualmap.submit.sh"
    
    echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --account=ag-demeaux
#SBATCH --error=/path/error_out_stuff/qual.%J.err
#SBATCH --output=/path/error_out_stuff/qual.%J.out

export PATH=/projects/ag-demeaux/ltntran/mytools/qualimap_v2.3:$PATH

cd /projects/ag-demeaux/ltntran/RNA_F2/

qualimap rnaseq -bam '"$BAM"' -gtf '"$GTF"' -outdir '"$REPORT_DIR/${base}"'
' > "$JOB_FILE"


############## End of the script ####################

## Submit all jobs by going to the JOBDIR and
for job in *.sh; do
sbatch $job
done

## Run multiqc to get a nice combined report and plots for all samples
multiqc ${REPORT_DIR} /path/multiqc_results


####### GOOD LUCK ##########
