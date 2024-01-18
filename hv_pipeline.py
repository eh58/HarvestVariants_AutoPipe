import argparse
import os
import subprocess
import time


# Create directory structure if it doesn't already exist
data_dir = "Data"
hv_dir = os.path.join(data_dir, "hv")
pango_dir = os.path.join(data_dir, "pango")
suggester_dir = os.path.join(data_dir, "suggester")
xml_dir = os.path.join(data_dir, "xml")
slurm_dir = os.path.join(data_dir, "slurm")

os.makedirs(hv_dir, exist_ok=True)
os.makedirs(pango_dir, exist_ok=True)
os.makedirs(suggester_dir, exist_ok=True)
os.makedirs(xml_dir, exist_ok=True)


start_time1 = time.time()
# Program 1: Run downloader.py and save original xml outputs to ./Data/xml
parser1 = argparse.ArgumentParser()
parser1.add_argument("--email", required=False)
parser1.add_argument("--api_key", required=False)
parser1.add_argument("--pdat", required=True)
# parser1.add_argument("--output", default="./Data/xml/")
parser1.add_argument("--suggester_jar", required=True)
args1 = parser1.parse_known_args()[0]

start_date, end_date = args1.pdat.split(":")
start_date, end_date = start_date.replace("/", ""), end_date.replace("/", "")
output_xml_filename = f"{start_date}_{end_date}_original.xml"
output_xml_file_path = os.path.join("./Data/xml/", output_xml_filename)

cmd1 = ["python", "src/downloader.py",
        "--email", args1.email,
        "--api_key", args1.api_key,
        "--pdat", args1.pdat]
subprocess.run(cmd1, check=True)

time.sleep(5)
elapsed_time1 = round(time.time() - start_time1, 2)
print("Downloader completed in", elapsed_time1, "seconds")


start_time2 = time.time()
# Program 2: Run xml_pp and save the pretty print xml outputs to ./Data/xml
input_pp_file_path = output_xml_file_path
output_pp_file_path = input_pp_file_path.replace("original", "pp")

cmd2 = ["xml_pp", input_pp_file_path]
with open(output_pp_file_path, "w") as outfile:
    subprocess.run(cmd2, check=True, stdout=outfile)

print(f"The original sra xml file has been converted to pretty print xml file in {output_pp_file_path}")
elapsed_time2 = round(time.time() - start_time2, 2)
print("XML file conversion completed in", elapsed_time2, "seconds")
time.sleep(5)


start_time3 = time.time()
# Program 3: Run run_suggester.py: trim the unused warning messages to avoid errors,
# and save the sra run id list txt files to ./Data/suggester
parser3 = argparse.ArgumentParser()
# parser3.add_argument("--output_dir", default="./Data/suggester/")
args3, unknown = parser3.parse_known_args()

cmd3 = ["python", "src/run_Suggester.py",
        args1.suggester_jar,
        output_pp_file_path]
subprocess.run(cmd3, check=True)

print(f"The ordered SRA IDs have been saved into ./Data/suggester/{start_date}_{end_date}.txt")
elapsed_time3 = round(time.time() - start_time3, 2)
print("run_Suggester completed in", elapsed_time3, "seconds")
time.sleep(5)


start_time4 = time.time()
# Program 4: Run hv in slurm and save the outputs to ./Data/hv/{current_hv_dir}
parser4 = argparse.ArgumentParser()
parser4.add_argument("--nStart", type=int, help="Starting index (inclusive) of IDs to process", default=None)
parser4.add_argument("--nEnd", type=int, help="Ending index (inclusive) of IDs to process", default=None)
args4, unknown = parser4.parse_known_args()

hv_program_path = "src/run_HVinSlurm.py"
input_sra_id_file_name = output_xml_filename.replace("_original.xml", ".txt")
input_sra_id_path = f"./Data/suggester/{input_sra_id_file_name}"
print(f"Test check: the input sra id file is {input_sra_id_path}")

cmd4 = ["python", hv_program_path,
        input_sra_id_path,
        "--nStart", str(args4.nStart),
        "--nEnd", str(args4.nEnd)]
subprocess.run(cmd4, check=True)
elapsed_time4 = round(time.time() - start_time4, 2)
print("run_harvest_variants in Slurm mode is completed in", elapsed_time4, "seconds")
time.sleep(5)


start_time5 = time.time()
# Program 5: Run Pangolin Command line tool and save the outputs to ./Data/pango
pangolin_program_path = "src/run_Pangolin.py"
task_subdir = f"{start_date}_{end_date}_{args4.nStart}to{args4.nEnd}"
# input_fasta_file_dir = f"./Data/hv/{task_subdir}/consensus_genomes/"

cmd5 = ["python", pangolin_program_path,
        task_subdir]
subprocess.run(cmd5, check=True)

print("Pangolin has been executed and SRA ID csv files have been saved into ./Data/pango")
elapsed_time5 = round(time.time() - start_time5, 2)
print("Pangolin completed in", elapsed_time5, "seconds")
time.sleep(5)


# Save processing times to a txt file
with open("processing_times.txt", "w") as file:
    file.write(f"Program 1: {elapsed_time1} seconds\n")
    file.write(f"Program 2: {elapsed_time2} seconds\n")
    file.write(f"Program 3: {elapsed_time3} seconds\n")
    file.write(f"Program 4: {elapsed_time4} seconds\n")
    file.write(f"Program 5: {elapsed_time5} seconds\n")

print("Processing times saved to processing_times.txt")