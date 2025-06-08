import scanpy as sc
import numpy as np
import os

TOTAL_DATA_PATH = "./raw_h5ad/merged_hasmeta_cells.h5ad"
SAVE_DIR = "./raw_h5ad"
DOWNSAMPLE_TRAIN = 20000
DOWNSAMPLE_TEST = 10000

os.makedirs(SAVE_DIR, exist_ok=True)

adata = sc.read_h5ad(TOTAL_DATA_PATH)
n_cells = adata.n_obs


if n_cells < DOWNSAMPLE_TRAIN + DOWNSAMPLE_TEST:
    raise ValueError(f"셀 수가 부족합니다. 총 셀 수: {n_cells}, 필요한 수: {DOWNSAMPLE_TRAIN + DOWNSAMPLE_TEST}")

all_indices = np.arange(n_cells)
train_indices = np.random.choice(all_indices, size=DOWNSAMPLE_TRAIN, replace=False)
remaining_indices = np.setdiff1d(all_indices, train_indices)
test_indices = np.random.choice(remaining_indices, size=DOWNSAMPLE_TEST, replace=False)

train_adata = adata[train_indices].copy()
test_adata = adata[test_indices].copy()

train_adata.write(os.path.join(SAVE_DIR, f"downsampled_train_{DOWNSAMPLE_TRAIN}.h5ad"))
test_adata.write(os.path.join(SAVE_DIR, f"downsampled_test_{DOWNSAMPLE_TEST}.h5ad"))

print("✅ 저장 완료:")
print(f" - train: downsampled_train_{DOWNSAMPLE_TRAIN}.h5ad")
print(f" - test : downsampled_test_{DOWNSAMPLE_TEST}.h5ad")