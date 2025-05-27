#!/bin/bash

# 1. 다운로드할 URL
url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE240822&format=file&file=GSE240822%5FGBM%5FccRCC%5FATAC%5Fmetadata%5FCPTAC%5Fsamples%2Etsv%2Egz"

# 2. 저장할 임시 파일 경로
downloaded="temp/GSE240822_GBM_ccRCC_ATAC_metadata_CPTAC_samples.tsv.gz"
final_output="temp/metadata.tsv"

# 3. temp 디렉토리 생성
mkdir -p temp

echo "📥 Downloading metadata file..."
wget -O "$downloaded" "$url"

echo "📂 Unzipping into $final_output ..."
gunzip -c "$downloaded" > "$final_output"

echo "🧹 Removing original .gz file..."
rm -f "$downloaded"

echo "✅ Done. Extracted file:"
ls -lh "$final_output"
