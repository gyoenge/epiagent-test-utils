import os
import subprocess
from epiagent.preprocessing import construct_cell_by_ccre_matrix
import pandas as pd

# Define paths
CCRE_FILE_PATH = "../cCRE.bed"
INPUT_DIR = "./fragment/"
OUTPUT_INTERSECT_DIR = "./temp/output_intersect/"
METADATA_PATH = './temp/tissue_positions_list.csv'
OUTPUT_DIR = "./raw_h5ad/"

"""
Fragment Overlap with cCREs
"""

# Ensure the output directory exists
os.makedirs(OUTPUT_INTERSECT_DIR, exist_ok=True)

# Iterate through .bed files in the input directory
for fragments_file in os.listdir(INPUT_DIR):
    if fragments_file.endswith(".bed"):
        basename = os.path.splitext(fragments_file)[0]
        output_file = os.path.join(OUTPUT_INTERSECT_DIR, f"{basename}.bed")

        # Logging
        print(f"Processing {fragments_file}...")

        # Construct bedtools command
        command = (
            f"bedtools intersect -a {os.path.join(INPUT_DIR, fragments_file)} "
            f"-b {CCRE_FILE_PATH} -wa -wb > {output_file}"
        )

        # Execute the command
        subprocess.run(command, shell=True, check=True)

print("Overlap calculation between fragments and cCREs completed.")

"""
Create AnnData from Intersect Results and cCRE Definitions
"""

# Process all intersect files
os.makedirs(OUTPUT_DIR, exist_ok=True)

intersect_files = os.listdir(OUTPUT_INTERSECT_DIR)
for intersect_file in intersect_files:
    if intersect_file.endswith(".bed"):
        output_file_path = os.path.join(OUTPUT_INTERSECT_DIR, intersect_file)
        output_filename = intersect_file.replace('.bed', '.h5ad')
        prefix = output_filename.split('.')[0]  # Remove file suffix
        final_output_path = os.path.join(OUTPUT_DIR, output_filename)

        # Logging
        print(f"Creating AnnData for {intersect_file}...")

        # Construct AnnData
        adata = construct_cell_by_ccre_matrix(output_file_path, CCRE_FILE_PATH)

        # Add metadata
        print("Adding metadata...")
        metadata = pd.read_csv(METADATA_PATH, index_col=0)
        print(f"| For {metadata.shape} metadata entries")
        
        adata.obs = adata.obs.join(metadata, how='left')

        # Save AnnData
        adata.write(final_output_path)
        print(f"Generated {final_output_path}")

print("AnnData creation from intersect results completed.")

