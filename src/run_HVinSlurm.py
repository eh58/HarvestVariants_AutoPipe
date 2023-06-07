import os
import time
from datetime import datetime
import subprocess
import argparse


def get_n_ids_from_list(input_file, ns, ne, output_file):
    with open(input_file, 'r') as f:
        ids = f.readlines()

    with open(output_file, 'w') as f:
        for i in range(ns - 1, ne):
            f.write(ids[i])


def split_n_ids_to_n_files(input_file, output_folder):
    with open(input_file, 'r') as f:
        ids = f.readlines()

    filenames = []
    for id_line in ids:
        id_name = id_line.strip()
        output_file = os.path.join(output_folder, f"{id_name}.txt")
        with open(output_file, 'w') as f:
            f.write(id_line)
        filenames.append(output_file)

    return filenames

def check_slurm_jobs_complete(job_name):
    completed = False
    while not completed:
        squeue_output = subprocess.check_output(["squeue", "-u", os.environ['USER'], "-n", job_name]).decode('utf-8')
        job_ids = squeue_output.split("\n")[1:-1]
        if len(job_ids) == 0:
            completed = True
        else:
            print(f"Sbatch Job {job_name} for user {os.environ['USER']} still has {len(job_ids)} running")
            time.sleep(600)  # Wait for 1 minute before checking again
    return completed

def main(args):
    original_file = args.original_sra_run_id_list
    ns = args.nStart
    ne = args.nEnd

    # Count the number of lines (ids) in the original file and handle out of boundary error
    with open(original_file, 'r') as f:
        total_ids = len(f.readlines())

    if ns is 0:
        ns = 1
    if ne is 0:
        ne = total_ids

    n = ne - ns + 1

    if n > total_ids or ne > total_ids:
        print("The id index range is out of bound. The tool will automatically end with the last SRA Run ID of the list.")
        ne = total_ids
        n = ne - ns + 1

    nsTOne_file = original_file.replace(".txt", f"_{ns}to{ne}.txt")
    get_n_ids_from_list(original_file, ns, ne, nsTOne_file)

    # Create directory for intermediate id files
    output_base_folder = os.path.dirname(original_file)
    output_folder_name = os.path.splitext(os.path.basename(nsTOne_file))[0]
    output_folder = os.path.join(output_base_folder, output_folder_name)
    os.makedirs(output_folder, exist_ok=True)

    split_n_ids_to_n_files(nsTOne_file, output_folder)

    # Create directory for slurm tasks
    # Get the current date as a string in the format yyyy-mm-dd
    current_date = datetime.now().strftime("%Y-%m-%d")
    # Create the directory name
    dir_name = f"{current_date}_slurm"
    # Specify the parent directory where the new directory will be created
    parent_dir = "./Data/slurm"
    # Create the full path for the new directory
    slurm_dir = os.path.join(parent_dir, dir_name)
    # Create the directory if it doesn't exist
    os.makedirs(slurm_dir, exist_ok=True)
    os.makedirs(os.path.join(slurm_dir, "logs"), exist_ok=True)

    # Create directory for harvest variants outputs
    hv_base_dir = "./Data/hv"
    current_hv_dir = os.path.join(hv_base_dir, output_folder_name)
    os.makedirs(current_hv_dir, exist_ok=True)

    slurm_job_name = "hvSlurm"
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=hvSlurm
#SBATCH --account=tgen
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1024M
#SBATCH --cpus-per-task=16
#SBATCH --output={slurm_dir}/logs/output_%j.log
#SBATCH --export=ALL
#SBATCH --array=1-{n}

PROGRAM_PATH="../sars-cov-2-harvest-variants/src/HarvestVariants/sra2vcf.py"
SINGLE_ID_FILE=$(ls {output_folder}/*.txt | awk "NR==${{SLURM_ARRAY_TASK_ID}}")
REFERENCE_FILE="../sars-cov-2-harvest-variants/data/reference_file.fasta"

PYTHON_EXECUTABLE=/home/Users/eh58/miniconda3/envs/harvest_variants/bin/python

python $PROGRAM_PATH -r $REFERENCE_FILE -s $SINGLE_ID_FILE -o {current_hv_dir}
"""

    with open(f"{slurm_dir}/run_slurm.sh", "w") as f:
        f.write(slurm_script)

    sbatch_output = subprocess.run(["sbatch", f"{slurm_dir}/run_slurm.sh"], capture_output=True)
    print(sbatch_output.stdout.decode("utf-8"))

    if check_slurm_jobs_complete(slurm_job_name):
        print(f"All Slurm tasks have been completed and the outputs are saved to {current_hv_dir}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare and run a program with Slurm")
    parser.add_argument("original_sra_run_id_list", type=str, help="Path to the original sra id list file")
    parser.add_argument("--nStart", type=int, help="Starting index (inclusive) of IDs to process", default=0)
    parser.add_argument("--nEnd", type=int, help="Ending index (inclusive) of IDs to process", default=0)
    # parser.add_argument("-o", "--hv_output_dir", type=str, required=True, help="Output directory for the HarvestVariants results")

    args = parser.parse_args()
    main(args)
