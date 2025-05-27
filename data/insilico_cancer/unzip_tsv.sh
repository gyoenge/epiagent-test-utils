#!/bin/bash

# 압축해제한 .tsv 파일을 저장할 디렉토리
output_dir="./fragment"

# fragment 디렉토리 없으면 생성
mkdir -p "$output_dir"

# 모든 .tsv.gz 파일 찾아서 반복
find . -type f -name "*.tsv.gz" | while read -r gzfile; do
  # 원본 파일명 추출 (경로 제거 + .gz 제거)
  base_name=$(basename "$gzfile" .gz)

  # 상위 상위 폴더 이름 추출
  grandparent_dir=$(basename "$(dirname "$(dirname "$gzfile")")")
  # 파일 이름 앞에 prefix 붙이기
  final_name="${grandparent_dir}.tsv"

  # 압축 해제 메시지 
  echo "Unzipping: $gzfile → $output_dir/$final_name"

  # 압축 해제하여 fragment 디렉토리에 저장
  gunzip -c "$gzfile" > "$output_dir/$final_name"
done
