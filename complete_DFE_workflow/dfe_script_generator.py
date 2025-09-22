# Define the template for the SLURM job script
# required python 3.9 and above
# recommended to use python poetry if having trouble installing dadi

job_script_template = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=150gb
#SBATCH --account=ag-demeaux
#SBATCH --error=/scratch/ntran5/errors/dfe.%J.err
#SBATCH --output=/scratch/ntran5/errors/dfe.%J.out

module load system/OpenSSL/1.1

source /projects/ag-demeaux/ltntran/poetry_env/bin/activate

cd /projects/ag-demeaux/ltntran/SFS_and_DFE/DFE_May2025

python dfe_April2025_VA.py --start {start} --end {end} --folder /projects/ag-demeaux/ltntran/SFS_and_DFE/DFE_May2025/INTRA/hig>
"""

# Define the range and step
start_range = 1
end_range = 200
step = 10

import os
os.chdir("/projects/ag-demeaux/ltntran/SFS_and_DFE/DFE_May2025/jobs_INTRA_highA")

# Loop to generate job scripts
for job_id, start in enumerate(range(start_range, end_range, step), start=1):
    end = start + step - 1
    job_script_content = job_script_template.format(job_id=job_id, start=start, end=end)

    # Write the job script to a file
    job_script_filename = f"job_script_{job_id}.sh"
    with open(job_script_filename, 'w') as job_script_file:
        job_script_file.write(job_script_content)

    print(f"Generated {job_script_filename}")
