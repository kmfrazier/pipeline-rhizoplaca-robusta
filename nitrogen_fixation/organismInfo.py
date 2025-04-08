#!/usr/bin/env python3
import os
import csv
from Bio import SeqIO

def parse_fasta_headers(fasta_dir):
    id_to_org = {}
    for root, _, files in os.walk(fasta_dir):
        for file in files:
            if file.endswith(".fasta") or file.endswith(".fa"):
                file_path = os.path.join(root, file)
                for record in SeqIO.parse(file_path, "fasta"):
                    parts = record.description.split("OS=")
                    if len(parts) > 1:
                        org = parts[1].split(" ")[0]
                    else:
                        org = "Unknown"
                    id_to_org[record.id] = org
    return id_to_org

def process_blastx_files(blast_dir, id_to_org, output_file):
    with open(output_file, "w", newline="") as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["Query_ID", "Subject_ID", "Percent_Identity", "Alignment_Length",
                         "Mismatches", "Gap_Openings", "Q_Start", "Q_End", "S_Start", "S_End",
                         "Evalue", "Bitscore", "Organism"])
        
        for file in os.listdir(blast_dir):
            if file.endswith(".txt"):
                file_path = os.path.join(blast_dir, file)
                with open(file_path, "r") as f:
                    for line in f:
                        fields = line.strip().split("\t")
                        if len(fields) < 12:
                            continue
                        subject_id = fields[1]
                        organism = id_to_org.get(subject_id, "Unknown")
                        writer.writerow(fields + [organism])

def main():
    fasta_path = "/grphome/grp_lichenscapstone/combined_db_uniprot"
    blastx_dir = "/grphome/grp_lichenscapstone/blastx3-5/blastx_nitrogenfixation"
    output_csv = os.path.join(blastx_dir, "annotated_blastx_results.csv")
    
    id_to_org = parse_fasta_headers(fasta_path)
    process_blastx_files(blastx_dir, id_to_org, output_csv)
    print("Finished annotating BLASTX results. Output written to:", output_csv)

if __name__ == "__main__":
    main()

