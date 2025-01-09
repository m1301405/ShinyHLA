#!/bin/bash

# 檢查是否提供輸入檔案
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <input_file> [output_file]"
    exit 1
fi

# 讀取輸入檔案和輸出檔案
INPUT_FILE=$1
OUTPUT_FILE=${2:-hla_coverage.csv}  # 如果未指定輸出檔案，默認為 hla_coverage.csv

# 檢查輸入檔案是否存在
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# 清空或初始化輸出檔案，並寫入表頭
HEADER=$(head -n 1 "$INPUT_FILE")
echo -e "HLA\t$HEADER" > "$OUTPUT_FILE"

# HLA-A
awk -v OFS='\t' '
    BEGIN {HLA="HLA-A"}
    NR > 1 && $1 == "chr6" && $2 > 29941260 && $3 < 29949572 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# HLA-B
awk -v OFS='\t' '
    BEGIN {HLA="HLA-B"}
    NR > 1 && $1 == "chr6" && $2 > 31353872 && $3 < 31367067 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# HLA-C
awk -v OFS='\t' '
    BEGIN {HLA="HLA-C"}
    NR > 1 && $1 == "chr6" && $2 > 31268749 && $3 < 31272130 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# HLA-DPA1
awk -v OFS='\t' '
    BEGIN {HLA="HLA-DPA1"}
    NR > 1 && $1 == "chr6" && $2 > 33064569 && $3 < 33080775 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# HLA-DPB1
awk -v OFS='\t' '
    BEGIN {HLA="HLA-DPB1"}
    NR > 1 && $1 == "chr6" && $2 > 33075936 && $3 < 33089696 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# HLA-DQA1
awk -v OFS='\t' '
    BEGIN {HLA="HLA-DQA1"}
    NR > 1 && $1 == "chr6" && $2 > 32628179 && $3 < 32647062 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# HLA-DQB1
awk -v OFS='\t' '
    BEGIN {HLA="HLA-DQB1"}
    NR > 1 && $1 == "chr6" && $2 > 32659467 && $3 < 32668383 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# HLA-DRB1
awk -v OFS='\t' '
    BEGIN {HLA="HLA-DRB1"}
    NR > 1 && $1 == "chr6" && $2 > 32577902 && $3 < 32589848 {
        print HLA, $0
    }
' "$INPUT_FILE" >> "$OUTPUT_FILE"

# 完成
echo "Results have been saved to $OUTPUT_FILE"
