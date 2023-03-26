import os
import time
import random
import argparse
from Bio import Entrez
from xml.etree.ElementTree import Element, SubElement, Comment, tostring, fromstring

# Argument parsing
parser = argparse.ArgumentParser(description="Download SRA XML data.")
parser.add_argument("--email", required=False, help="Your email address.")
parser.add_argument("--api_key", required=False, help="Your NCBI API key.")
parser.add_argument("--pdat", required=True, help="Publication date for the query.")
parser.add_argument("--output_dir", default="./", help="Output directory for the XML file.")
args = parser.parse_args()

# Optional Entrez Authentication
Entrez.email = args.email
Entrez.api_key = args.api_key

# Database and query settings
db = "sra"
search_query = f'txid2697049[Organism] AND "illumina"[Platform] {args.pdat}[PDAT]'
batch_size = 10000

# Search record count and query key, and print out the total record number
search_handle = Entrez.esearch(db=db, term=search_query, usehistory="y")
search_results = Entrez.read(search_handle)
search_handle.close()

num_records = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

print(f"Total records to download: {num_records}")

# Calculate the number of batches
num_batches = (num_records + batch_size - 1) // batch_size

# Initialize the concatenated XML data
root = Element("EXPERIMENT_PACKAGE_SET")

# Function to handle connection failure issue
def fetch_records(start, end, retries=10):
    for _ in range(retries):
        try:
            fetch_handle = Entrez.efetch(
                db=db,
                rettype="full",
                retmode="xml",
                retstart=start,
                retmax=end - start,
                webenv=webenv,
                query_key=query_key,
            )
            fetch_data = fetch_handle.read().decode("utf-8")
            fetch_handle.close()
            return fetch_data
        except Exception as e:
            print(f"Error occurred while fetching records {start + 1} to {end}: {e}")
            if _ < retries - 1:
                wait_time = 5 + random.uniform(0, 5)
                print(f"Retrying in {wait_time:.2f} seconds...")
                time.sleep(wait_time)
            else:
                print(f"Failed to fetch records {start + 1} to {end} after {retries} attempts.")
                return None

# Download the XML files in batches and concatenate with connection failure function
for i in range(num_batches):
    start = i * batch_size
    end = min((i + 1) * batch_size, num_records)

    print(f"Downloading records {start + 1} to {end}...")

    fetch_data = fetch_records(start, end)
    if fetch_data is None:
        print(f"Skipping records {start + 1} to {end} due to fetch error.")
        continue

    # Parse the fetched XML data
    fetched_tree = fromstring(fetch_data)

    # Merge the fetched XML data with the root element
    for experiment_package in fetched_tree.findall("EXPERIMENT_PACKAGE"):
        root.append(experiment_package)

    # Wait 1 second between requests to avoid overwhelming the server
    time.sleep(1)

# Extract the start and end dates from the pdat argument and generate the output filename
start_date, end_date = args.pdat.split(":")
start_date, end_date = start_date.replace("/", ""), end_date.replace("/", "")
output_filename = f"{start_date}_{end_date}_original.xml"
output_file_path = os.path.join(args.output_dir, output_filename)

# Save the concatenated XML data to the specified output file
with open(output_file_path, "w", encoding="utf-8") as xml_file:
    xml_file.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
    xml_file.write(tostring(root, encoding="unicode"))

print(f"Download complete. Concatenated XML data saved to {output_file_path}.")
