import scanpy as sc
import numpy as np
import os
import re

TOTAL_DATA_PATH = "./raw_h5ad/merged_hasmeta_cells.h5ad"
SAVE_DIR = "./raw_h5ad"
DOWNSAMPLE_SET = [
    "C3L-00088-T2_CPT0000880001_snATAC_ccRCC",
    "C3L-00088-N_CPT0000890002_snATAC_ccRCC",
    "C3L-00088-T1_CPT0000870003_snATAC_ccRCC",
]
# DOWNSAMPLE_SET = [
#     # blue
#     "C3L-00088-T2_CPT0000880001_snATAC_ccRCC",
#     "C3L-00088-N_CPT0000890002_snATAC_ccRCC",
#     "C3L-00088-T1_CPT0000870003_snATAC_ccRCC",
#     # orange
#     "C3N-00242-T1_CPT0014450005_snATAC_ccRCC",
#     "C3N-00242-N_CPT0014470002_snATAC_ccRCC",
#     # green
#     "C3L-00079-T1_CPT0001260013_snATAC_ccRCC",
#     "C3L-00079-N_CPT0001270002_snATAC_ccRCC",
#     # purple
#     "C3N-01200-T1_CPT0075130004_snATAC_ccRCC",
#     "C3N-01200-N_CPT0075170013_snATAC_ccRCC",
# ]
# DOWNSAMPLE_SET = [
#     # paper
#     "C3N-00242-T1_CPT0014450005_snATAC_ccRCC",
#     "C3N-01200-T1_CPT0075130004_snATAC_ccRCC",
#     "C3L-00004-T1_CPT0001540013_snATAC_ccRCC",
#     "C3L-00610-T1_CPT0025110004_snATAC_ccRCC",
# ]

os.makedirs(SAVE_DIR, exist_ok=True)

adata = sc.read_h5ad(TOTAL_DATA_PATH)

pattern = '|'.join(map(re.escape, DOWNSAMPLE_SET))
mask = adata.obs['GEO.sample'].str.contains(pattern, na=False)

filtered_adata = adata[mask].copy()
n_sample = filtered_adata.shape[0]

output_path = os.path.join(SAVE_DIR, f"downsampled_stable_{n_sample}.h5ad")
# output_path = os.path.join(SAVE_DIR, f"downsampled_paper_{n_sample}.h5ad")
filtered_adata.write(output_path)

print(f"✅ Saved downsampled data: {output_path}")