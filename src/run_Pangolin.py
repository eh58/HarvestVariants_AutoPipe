import os
import argparse

# Set the directory containing the input FASTA files from parser
parser = argparse.ArgumentParser(description="Run pangolin for all the SRA IDs from given directory")
parser.add_argument("task_subdir", type = str, help = "Specify the subdirectory of Harvest Variants for this task")
args = parser.parse_args()

# Get the input directory from harvest variants for pangolin
input_dir = f"Data/hv/{args.task_subdir}/consensus_genomes"

# Set the output directory for the pango results
output_path = "Data/pango"
output_dir = os.path.join(output_path, args.task_subdir)
os.makedirs(output_dir, exist_ok=True)

# Loop through all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta"):
        # Build the full path to the input file
        input_file = os.path.join(input_dir, filename)

        # Build the output filename by replacing ".fasta" with ".csv"
        output_file = filename.replace(".fasta", ".csv")

        # Run pangolin with the input and output paths
        os.system("pangolin {} -o {} --outfile {}".format(input_file, output_dir, output_file))
