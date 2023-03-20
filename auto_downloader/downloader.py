from Bio import Entrez
from tqdm import tqdm
import os
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(description="Download SRA records from NCBI")
parser.add_argument("pdat", type=str, help="Publication date range (e.g., 2022/05/20:2022/06/03)")
parser.add_argument("file_location", type=str, help="Location to save the downloaded file")
# Optional argument parser
parser.add_argument("--email", type=str, help="Email address for Entrez (optional)")
parser.add_argument("--api_key", type=str, help="API key for Entrez (optional)")

# Parse the arguments
args = parser.parse_args()
pdat = args.pdat
file_location = args.file_location

# Set the email and API key for Entrez
Entrez.email = args.email
Entrez.api_key = args.api_key

# Set the search parameters
search_query = f"txid2697049[Organism] AND illumina[Platform] AND {pdat}[PDAT]"

# Set the batch size
batch_size = 1000

# Get the total number of records for the search query
handle = Entrez.esearch(db="sra", term=search_query, retmax=0, retmode="xml")
records = Entrez.read(handle)
total_count = int(records["Count"])
print(f"Total number of records for search query: {total_count}")
handle.close()

# Download the search results in batches
for batch_start in tqdm(range(0, total_count, batch_size)):
    # Set the batch end
    batch_end = min(total_count, batch_start + batch_size)

    # Download the batch
    batch_handle = Entrez.esearch(db="sra", term=search_query, retstart=batch_start, retmax=batch_size, idtype="acc", retmode="xml")
    batch_records = Entrez.read(batch_handle)
    batch_handle.close()

    # Get the list of IDs for the batch
    id_list = batch_records["IdList"]

    # Download the full records for the batch
    batch_handle = Entrez.efetch(db="sra", id=id_list, rettype="full", retmode="xml")
    batch_xml = batch_handle.read()
    batch_handle.close()

    # Save the batch to file
    with open(file_location, "wb") as f:
        f.write(batch_xml)

    tqdm.write(f"Downloaded batch {batch_start + 1} to {batch_end} of {total_count}")

print("All batches downloaded and saved to {file_location}")
