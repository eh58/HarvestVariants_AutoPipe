import os
import argparse

# Set the directory containing the input FASTA files from parser
parser = argparse.ArgumentParser(description="Run pangolin for all the SRA IDs from given directory")
parser.add_argument("input_dir", type = str, help = "Specify the subdirectory of Harvest Variants for this task")
args = parser.parse_args()

# Set the output directory for the pango results
output_dir = "Data/pango/"

# Loop through all files in the input directory
for filename in os.listdir(args.input_dir):
    if filename.endswith(".fasta"):
        # Build the full path to the input file
        input_file = os.path.join(args.input_dir, filename)

        # Build the output filename by replacing ".fasta" with ".csv"
        output_file = filename.replace(".fasta", ".csv")

        # Run pangolin with the input and output paths
        os.system("pangolin {} -o {} --outfile {}".format(input_file, output_dir, output_file))
