import pysam
from Bio import SeqIO
import argparse

def process_bam_and_fasta(bam_file, reference_fasta, output_fasta):
    
    hla_identifiers = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            hla_identifier = read.reference_name  
            hla_identifiers.add(hla_identifier)

    
    with open(output_fasta, "w") as out_fasta:
        for record in SeqIO.parse(reference_fasta, "fasta"):
            header = record.id.split()[0]
            if header in hla_identifiers:  
                SeqIO.write(record, out_fasta, "fasta")  

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Extract HLA identifiers from BAM file and match them to reference FASTA.")
    parser.add_argument("input_bam", help="Path to the input BAM file")
    parser.add_argument("input_fasta", help="Path to the input reference FASTA file")
    parser.add_argument("output_fasta", help="Path to the output FASTA file")
    args = parser.parse_args()

    
    process_bam_and_fasta(args.input_bam, args.input_fasta, args.output_fasta)

