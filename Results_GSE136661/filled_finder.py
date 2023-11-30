# Given a CSV that looks like:
# SAMN,Cohort,Pathology,Tissue,Sex,ID,Organism,Age,SRR
# SAMN12661878,Discovery,WHO I,Meningioma,female,12661878,Homo sapiens,35,SRR10042693
# SAMN12661871,Discovery,WHO II,Meningioma,female,12661871,Homo sapiens,69,SRR10042694
# ...,...,...,etc
# 

import csv

# Function to read SRR values from the txt file
def read_srr_from_file(file_path):
    with open(file_path, 'r') as file:
        return set(line[2:].strip() for line in file)

# Function to filter rows based on SRR values
def filter_csv(input_csv_path, output_csv_path, srr_values):
    with open(input_csv_path, 'r') as input_file, open(output_csv_path, 'w', newline='') as output_file:
        csv_reader = csv.reader(input_file)
        csv_writer = csv.writer(output_file)

        # Write the header to the new CSV file
        header = next(csv_reader)
        csv_writer.writerow(header)

        # Filter rows based on SRR values
        for row in csv_reader:
            if row[-1] in srr_values:
                csv_writer.writerow(row)

# Path to the old CSV file
old_csv_path = 'GSE136661.csv'

# Path to the txt file with SRR values
srr_txt_path = 'GSE136661.txt'

# Path to the new CSV file
new_csv_path = 'GSE136661r.csv'

# Read SRR values from the txt file
srr_values = read_srr_from_file(srr_txt_path)

# Filter and create a new CSV file
filter_csv(old_csv_path, new_csv_path, srr_values)

print("New CSV file created successfully.")
