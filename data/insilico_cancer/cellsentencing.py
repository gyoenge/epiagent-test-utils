from epiagent.preprocessing import global_TFIDF
from epiagent.tokenization import tokenization
import numpy as np
import scanpy as sc
import os

CCRE_DOCFREQ_PATH = '../cCRE_document_frequency.npy'
DATA_PATH = "./raw_h5ad/downsampled_train_20000.h5ad"
# DATA_PATH = "./raw_h5ad/downsampled_test_10000.h5ad"
SAVE_DIR = "./processed_h5ad" 

adata = sc.read_h5ad(DATA_PATH)
cCRE_document_frequency = np.load(CCRE_DOCFREQ_PATH)

# Apply TF-IDF
print("Applying TF-IDF...")
global_TFIDF(adata, cCRE_document_frequency)

# Tokenize the data
print("Tokenizing the data...")
tokenization(adata)

# Save the processed AnnData
os.makedirs(SAVE_DIR, exist_ok=True)
processed_output_path = os.path.join(SAVE_DIR, os.path.basename(DATA_PATH).replace('.h5ad', '_cellsentenced.h5ad'))
adata.write(processed_output_path)
print(f"Processed data saved at {processed_output_path}")
