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

os.makedirs(hv_dir, exist_ok=True)
os.makedirs(pango_dir, exist_ok=True)
os.makedirs(suggester_dir, exist_ok=True)
os.makedirs(xml_dir, exist_ok=True)

# Program 1
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

cmd1 = ["python", "auto_downloader/downloader.py",
        "--email", args1.email,
        "--api_key", args1.api_key,
        "--pdat", args1.pdat]
subprocess.run(cmd1, check=True)

time.sleep(5)

# Program 2
input_pp_file_path = output_xml_file_path
output_pp_file_path = input_pp_file_path.replace("original", "pp")

cmd2 = ["xml_pp", input_pp_file_path]
with open(output_pp_file_path, "w") as outfile:
    subprocess.run(cmd2, check=True, stdout=outfile)

print(f"The original sra xml file has been converted to pretty print xml file in {output_pp_file_path}")
time.sleep(5)

# Program 3
parser3 = argparse.ArgumentParser()
# parser3.add_argument("--output_dir", default="./Data/suggester/")
args3, unknown = parser3.parse_known_args()

cmd3 = ["python", "run_Suggester/run_Suggester.py",
        args1.suggester_jar,
        output_pp_file_path]
subprocess.run(cmd3, check=True)

print(f"The ordered SRA IDs haven been saved into {output_pp_file_path}")
time.sleep(5)
