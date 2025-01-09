from Bio import SeqIO
import argparse

def extract_sequences(input_list_file, fasta_file, output_file):
    
    with open(input_list_file, "r") as bam_file:
        sequences_from_bam = bam_file.read().splitlines()

    
    with open(output_file, "w") as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.id.split()[0]  
            if header in sequences_from_bam:
                SeqIO.write(record, out_fasta, "fasta")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Extract sequences from a list file and match them to a reference FASTA.")
    parser.add_argument("input_list", help="Path to the input list file (.txt)")
    parser.add_argument("input_fasta", help="Path to the input FASTA file")
    parser.add_argument("output_fasta", help="Path to the output FASTA file")
    args = parser.parse_args()

    
    extract_sequences(args.input_list, args.input_fasta, args.output_fasta)
