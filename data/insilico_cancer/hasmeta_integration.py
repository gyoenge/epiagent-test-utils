import glob, os
import scanpy as sc
import anndata as ad

SAVE_DIR = "./raw_h5ad"
os.makedirs(SAVE_DIR, exist_ok=True)

h5ad_files = glob.glob(os.path.join("./raw_h5ad", "*.h5ad"))

filtered_adatas = []
total_cells = 0
meta_cells = 0

for i, h5ad_path in enumerate(h5ad_files):
    print(f"Processing {i+1} : {h5ad_path}")
    h5ad_adata = sc.read_h5ad(h5ad_path)

    # 전체 셀 수
    total = h5ad_adata.shape[0]

    # Merged_barcode 있는 셀 필터
    valid_mask = ~h5ad_adata.obs['Merged_barcode'].isna() & (h5ad_adata.obs['Merged_barcode'].astype(str).str.lower() != 'nan')
    filtered = h5ad_adata[valid_mask].copy()

    meta = filtered.shape[0]
    print(f"| Total : {total} cells")
    print(f"| HasMeta : {meta} cells")

    # 누적 합계
    total_cells += total
    meta_cells += meta

    if meta > 0:
        filtered_adatas.append(filtered)

# 병합 및 저장
if filtered_adatas:
    print("Merging filtered cells into one AnnData...")
    merged_adata = ad.concat(filtered_adatas, axis=0)
    merged_adata.write(os.path.join(SAVE_DIR, "merged_hasmeta_cells.h5ad"))
    print("✅ Saved to: merged_hasmeta_cells.h5ad")
else:
    print("⚠️ No metadata-containing cells found. Nothing to save.")

# 통계 출력
print(f"Total cells across all datasets: {total_cells}")
print(f"Total cells with metadata across all datasets: {meta_cells}")
print(f"Total cells with metadata ratio: {meta_cells / total_cells:.2%}" if total_cells > 0 else "No cells found.")
