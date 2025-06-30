#!/bin/bash
#SBATCH --job-name=mafft_align_array
#SBATCH --output=/your/path/mafft_%A_%a.out
#SBATCH --error=/your/path/mafft_%A_%a.err
#SBATCH --array=1-300
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --mem=4G

# Move to working directory 
cd /your/path

# Load necessary modules
module load lang/Miniconda3/23.9.0-0
module load lang/Python/3.11.3-GCCcore-12.3.0
module load lang/Perl/5.36.1-GCCcore-12.3.0

# Activate your conda environment
conda activate /your/conda/env

# Gather all input fasta files into an array
input_files=(/your_path_to_cds_pairs/*.fa)
total_files=${#input_files[@]}

# Calculate number of files per task (distributes roughly evenly among 300 tasks)
files_per_task=$(( (total_files + 299) / 300 ))

# Calculate this task's start and end indices
start_index=$(( ($SLURM_ARRAY_TASK_ID - 1) * files_per_task ))
end_index=$(( $SLURM_ARRAY_TASK_ID * files_per_task - 1 ))
if [ $end_index -ge $total_files ]; then
    end_index=$((total_files-1))
fi

# Process this chunk of files
for i in $(seq $start_index $end_index); do
    input_file="${input_files[$i]}"
    # Only process if the file actually exists (avoid out-of-bounds issues)
    if [[ -f "$input_file" ]]; then
        python step3_mafft_alignment_pal2nal.py "$input_file"
    fi
done
