#!/bin/bash

# 현재 디렉토리에서 모든 .tar.gz 파일을 찾음
for file in *.tar.gz; do
  echo "Extracting: $file"

  # tar.gz 압축 해제
  tar -xvzf "$file"

  # 원본 .tar.gz 파일은 필요 시 삭제
  rm "$file"
done
