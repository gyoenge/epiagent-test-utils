import scanpy as sc 
import pandas as pd

DATA_PATH = "./raw_h5ad/merged_hasmeta_cells.h5ad"
METADATA_PATH = "./temp/normal_metadata.tsv"
OUTPUT_PATH = "./raw_h5ad/normal_withmeta.h5ad"

adata = sc.read_h5ad(DATA_PATH)
adata.obs.index = adata.obs['Merged_barcode'].tolist()

print("Adding metadata...")
metadata = pd.read_csv(METADATA_PATH, sep='\t', index_col=0)
print(f"| For {metadata.shape} metadata entries")
metadata = metadata[metadata['cell_type.normal'].astype(str) == 'Proximal Tubule'] # Filter metadata for 'Proximal Tubule' cell type
print(f"| Selected {metadata.shape}")

adata = adata[adata.obs_names.isin(metadata.index)].copy()
metadata = metadata.loc[adata.obs_names]
adata.obs = adata.obs.join(metadata, how='right')
adata.obs = adata.obs.astype(str) 
print(f"| Final adata shape: {adata.shape}")

# Save AnnData
adata.write(OUTPUT_PATH)
print(f"Generated {OUTPUT_PATH}")

