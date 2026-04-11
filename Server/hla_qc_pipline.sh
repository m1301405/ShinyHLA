#!/bin/bash

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <input_file> [output_file]"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=${2:-hla_coverage.csv}

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

awk -v OFS='\t' '
BEGIN {
    genes["HLA-A"]    = "29941260 29949572"
    genes["HLA-B"]    = "31353872 31367067"
    genes["HLA-C"]    = "31268749 31272130"
    genes["HLA-DPA1"] = "33064569 33080775"
    genes["HLA-DPB1"] = "33075936 33089696"
    genes["HLA-DQA1"] = "32628179 32647062"
    genes["HLA-DQB1"] = "32659467 32668383"
    genes["HLA-DRB1"] = "32577902 32589848"
}

NR == 1 {
    $1 = "Chr"
    print "HLA", $0
    next
}

{
    chr   = $1
    start = $2
    end   = $3

    if (chr != "chr6" && chr != "6") next

    for (gene in genes) {
        split(genes[gene], pos, " ")
        gstart = pos[1]
        gend   = pos[2]

        if (start < gend && end > gstart) {
            print gene, $0
        }
    }
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Results have been saved to $OUTPUT_FILE"