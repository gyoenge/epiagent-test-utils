#!/bin/bash

# mkdir -p bed

for tsv in fragment/*.tsv; do
  name=$(basename "$tsv" .tsv)
  out="fragment/${name}.bed"

  echo "Converting: $tsv → $out"

  # cut -f1-5 "$tsv" > "$out"
  # 주석(#) 제거 + 1~5열만 추출
  grep -v '^#' "$tsv" | cut -f1-5 > "$out"
done