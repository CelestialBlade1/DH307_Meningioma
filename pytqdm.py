import csv
import re
import subprocess
from tqdm import tqdm

# Function to get SRR number for a given SAMN number
def get_srr_number(samn_number):
    # Run esearch command
    command = f"esearch -db biosample -query {samn_number} | elink -target sra | efetch -format docsum | xtract -pattern Runs -element Run@acc"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout.strip()

# Read data from the file
with open('metadata.txt', 'r') as file:
    data = file.read()

# Initialize CSV data
csv_data = []

# Define CSV headers
headers = ["SAMN", "Cohort", "Pathology", "Tissue", "Sex", "ID", "Organism", "Age", "SRR"]

# Compile regular expressions for matching patterns
sample_pattern = re.compile(r'BioSample: (\w+);')
cohort_pattern = re.compile(r'/cohort="([^"]+)"')
pathology_pattern = re.compile(r'/pathology="([^"]+)"')
tissue_pattern = re.compile(r'/tissue="([^"]+)"')
sex_pattern = re.compile(r'/sex="([^"]+)"')
id_pattern = re.compile(r'ID: (\d+)')
organism_pattern = re.compile(r'Organism: (.+)')
age_pattern = re.compile(r'/age="([^"]+)"')

# Split the data into individual records
records = data.split('\n\n')

# Create a tqdm progress bar
progress_bar = tqdm(records, desc="Processing records", unit="record")

# Extract information for each record
for record in progress_bar:
    match_sample = sample_pattern.search(record)
    match_cohort = cohort_pattern.search(record)
    match_pathology = pathology_pattern.search(record)
    match_tissue = tissue_pattern.search(record)
    match_sex = sex_pattern.search(record)
    match_id = id_pattern.search(record)
    match_organism = organism_pattern.search(record)
    match_age = age_pattern.search(record)

    if match_sample and match_cohort and match_pathology and match_tissue and match_sex and match_id and match_organism and match_age:
        samn_number = match_sample.group(1)
        srr_number = get_srr_number(samn_number)

        csv_data.append([
            samn_number,
            match_cohort.group(1),
            match_pathology.group(1),
            match_tissue.group(1),
            match_sex.group(1),
            match_id.group(1),
            match_organism.group(1),
            match_age.group(1),
            srr_number,
        ])

# Close the tqdm progress bar
progress_bar.close()

# Write to CSV file
with open('output__.csv', 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(headers)
    csv_writer.writerows(csv_data)

print("CSV file created successfully.")

