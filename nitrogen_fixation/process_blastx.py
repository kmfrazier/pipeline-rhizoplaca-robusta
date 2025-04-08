#!/usr/bin/env python3
import glob
import csv
import os

# Pattern to match BLASTX output files; adjust if needed.
blast_files = glob.glob("*_blastx.txt")

# Data structure to store aggregated results per query
query_summary = {}

# Process each file found
for filename in blast_files:
    print(f"Processing file: {filename}")
    with open(filename, 'r') as fh:
        for line in fh:
            # Skip comment or empty lines
            if line.startswith("#") or not line.strip():
                continue
            # Split line into columns (assuming tab-delimited)
            fields = line.strip().split("\t")
            if len(fields) < 12:
                # If the line does not have enough fields, skip it
                continue

            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = fields[:12]
            try:
                pident = float(pident)
                length = int(length)
                mismatch = int(mismatch)
                gapopen = int(gapopen)
                qstart = int(qstart)
                qend = int(qend)
                sstart = int(sstart)
                send = int(send)
                evalue = float(evalue)
                bitscore = float(bitscore)
            except ValueError:
                # If conversion fails, skip this line
                continue

            # Aggregate information per query sequence id
            if qseqid not in query_summary:
                query_summary[qseqid] = {
                    "hits": 0,
                    "min_evalue": evalue,
                    "max_bitscore": bitscore,
                    "files": set()
                }
            query_summary[qseqid]["hits"] += 1
            query_summary[qseqid]["min_evalue"] = min(query_summary[qseqid]["min_evalue"], evalue)
            query_summary[qseqid]["max_bitscore"] = max(query_summary[qseqid]["max_bitscore"], bitscore)
            query_summary[qseqid]["files"].add(filename)

# Write out a summary CSV file
output_csv = "blastx_summary.csv"
with open(output_csv, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    # Write header row
    writer.writerow(["Query_ID", "Hit_Count", "Min_Evalue", "Max_Bitscore", "Source_Files"])
    for qseqid, summary in query_summary.items():
        writer.writerow([
            qseqid,
            summary["hits"],
            summary["min_evalue"],
            summary["max_bitscore"],
            ";".join(sorted(summary["files"]))
        ])

print(f"Summary written to {output_csv}")

