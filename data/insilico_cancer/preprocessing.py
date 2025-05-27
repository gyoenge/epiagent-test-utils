import os
import subprocess
from epiagent.preprocessing import construct_cell_by_ccre_matrix
import pandas as pd

# Define paths
CCRE_FILE_PATH = "../cCRE.bed"
INPUT_DIR = "./fragment/"
OUTPUT_INTERSECT_DIR = "./temp/output_intersect/"
METADATA_PATH = './temp/metadata.tsv'
OUTPUT_DIR = "./raw_h5ad/"

"""
Fragment Overlap with cCREs
"""

# # Ensure the output directory exists
# os.makedirs(OUTPUT_INTERSECT_DIR, exist_ok=True)

# # Iterate through .bed files in the input directory
# for fragments_file in os.listdir(INPUT_DIR):
#     if fragments_file.endswith(".bed"):
#         basename = os.path.splitext(fragments_file)[0]
#         output_file = os.path.join(OUTPUT_INTERSECT_DIR, f"{basename}.bed")

#         # Logging
#         print(f"Processing {fragments_file}...")

#         # Construct bedtools command
#         command = (
#             f"bedtools intersect -a {os.path.join(INPUT_DIR, fragments_file)} "
#             f"-b {CCRE_FILE_PATH} -wa -wb > {output_file}"
#         )

#         # Execute the command
#         subprocess.run(command, shell=True, check=True)

# print("Overlap calculation between fragments and cCREs completed.")

"""
Create AnnData from Intersect Results and cCRE Definitions
"""

# Process all intersect files
os.makedirs(OUTPUT_DIR, exist_ok=True)

# intersect_files = os.listdir(OUTPUT_INTERSECT_DIR)
intersect_files = [
    'C3L-00908-T1_CPT0086350004_snATAC_ccRCC.bed', 
    'C3L-00917-T1_CPT0023690004_snATAC_ccRCC.bed', 
    'C3L-01287-T1_CPT0079410004_snATAC_ccRCC.bed', 
    'C3L-01302-T1_CPT0063630004_snATAC_ccRCC.bed', 
    'C3L-01313-T1_CPT0086820004_snATAC_ccRCC.bed', 
    # 'C3N-00242-N_CPT0014470002_snATAC_ccRCC.bed', 
    # 'C3N-00242-T1_CPT0014450005_snATAC_ccRCC.bed', 
    # 'C3N-00317-T1_CPT0012280004_snATAC_ccRCC.bed', 
    # 'C3N-00437-T1_CPT0012550012_snATAC_ccRCC.bed', 
    # 'C3N-00495-T1_CPT0078510004_snATAC_ccRCC.bed', 
    # 'C3N-00733-T1_CPT0025880013_snATAC_ccRCC.bed', 
    # 'C3N-01200-N_CPT0075170013_snATAC_ccRCC.bed', 
    # 'C3N-01200-T1_CPT0075130004_snATAC_ccRCC.bed', 
    # 'C3N-01213-T1_CPT0075720013_snATAC_ccRCC.bed', 
    # 'C3L-00004-T1_CPT0001540013_snATAC_ccRCC.bed', 
    # 'C3L-00010-T1_CPT0001220012_snATAC_ccRCC.bed', 
    # 'C3L-00026-T1_CPT0001500003_snATAC_ccRCC.bed', 
    # 'C3L-00079-N_CPT0001270002_snATAC_ccRCC.bed', 
    # 'C3L-00079-T1_CPT0001260013_snATAC_ccRCC.bed', 
    # 'C3L-00088-N_CPT0000890002_snATAC_ccRCC.bed', 
    # 'C3L-00088-T1_CPT0000870003_snATAC_ccRCC.bed', 
    # 'C3L-00088-T2_CPT0000880001_snATAC_ccRCC.bed', 
    # 'C3L-00096-T1_CPT0001180011_snATAC_ccRCC.bed', 
    # 'C3L-00416-T2_CPT0010100001_snATAC_ccRCC.bed', 
    # 'C3L-00448-T1_CPT0010160004_snATAC_ccRCC.bed', 
    # 'C3L-00583-T1_CPT0019130004_snATAC_ccRCC.bed', 
    # 'C3L-00610-T1_CPT0025110004_snATAC_ccRCC.bed', 
    # 'C3L-00790-T1_CPT0065690004_snATAC_ccRCC.bed'
]
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
        metadata = pd.read_csv(METADATA_PATH, sep='\t', index_col=1)
        print(f"| For {metadata.shape} metadata entries")
        metadata = metadata[metadata['GEO.sample'].astype(str).str.contains(prefix)] # Filter metadata for current sample
        # metadata = metadata[metadata['Cancer'].astype(str)=='ccRCC'] # Filter metadata for ccRCC
        # metadata = metadata[~metadata.index.duplicated(keep='first')] # Remove duplicates based on index
        print(f"| Selected {metadata.shape}")
        
        # print(f"| {metadata['Merged_barcode'].unique().shape[0]}")
        # print(f"| {metadata['Merged_barcode'].unique()}")
        # print(f"| {metadata.index.unique().shape[0]}")
        # print(f"| {metadata.index.unique()}")

        adata.obs = adata.obs.join(metadata, how='left')
        adata.obs.index = prefix + '_' + adata.obs.index
        adata.obs = adata.obs.astype(str) # Ensure obs is of type str

        # Save AnnData
        adata.write(final_output_path)
        print(f"Generated {final_output_path}")

print("AnnData creation from intersect results completed.")

