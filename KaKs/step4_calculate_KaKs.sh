#!/bin/bash
#SBATCH --job-name=kaks_array
#SBATCH --output=/your_path/kaks_%A_%a.out
#SBATCH --error=/your_path/kaks_%A_%a.err
#SBATCH --array=1-300
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=2G

cd /your/path

module load lang/Miniconda3/23.9.0-0
module load lang/Python/3.11.3-GCCcore-12.3.0
module load lang/Perl/5.36.1-GCCcore-12.3.0

conda activate /path/to/your/conda_env

# List all codon alignments to process
input_files=(/path_to/aligned_pairs/*.codon.aln.fa)
total_files=${#input_files[@]}
files_per_task=$(( (total_files + 299) / 300 ))

start_index=$(( ($SLURM_ARRAY_TASK_ID - 1) * files_per_task ))
end_index=$(( $SLURM_ARRAY_TASK_ID * files_per_task - 1 ))
if [ $end_index -ge $total_files ]; then
    end_index=$((total_files-1))
fi

#mkdir -p kaks_results
for i in $(seq $start_index $end_index); do
    input_file="${input_files[$i]}"
    if [[ -f "$input_file" ]]; then
        base=$(basename "$input_file" .codon.aln.fa)
        KaKs_Calculator -i "$input_file" -o "/your/path/to/kaks_results/${base}.kaks.txt" -m YN
    fi
done
