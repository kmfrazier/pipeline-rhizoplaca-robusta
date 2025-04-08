#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, default="blastx_summary.csv")
    parser.add_argument("--output", type=str, default="blastx_summary_graph.png")
    parser.add_argument("--threshold", type=float, default=40.0)
    args = parser.parse_args()

    # Use a built-in style that is available
    plt.style.use("ggplot")

    # Load the CSV file
    df = pd.read_csv(args.input)
    threshold = args.threshold

    # Split the data into "good" hits and others
    good_hits = df[df["Max_Bitscore"] >= threshold]
    other_hits = df[df["Max_Bitscore"] < threshold]

    # Create the scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(good_hits["Hit_Count"], good_hits["Max_Bitscore"],
                color="forestgreen", alpha=0.8, s=50, edgecolor="black", label="Good hits")
    plt.scatter(other_hits["Hit_Count"], other_hits["Max_Bitscore"],
                color="lightgray", alpha=0.8, s=50, edgecolor="black", label="Other hits")
    plt.axhline(y=threshold, color="red", linestyle="--", linewidth=2, label=f"Threshold (Bitscore = {threshold})")

    plt.xlabel("Number of Hits per Query", fontsize=12)
    plt.ylabel("Maximum Bitscore", fontsize=12)
    plt.title("Summary of BLASTX Hits for Nitrogen Fixation Genes", fontsize=14)
    plt.legend(loc="best")
    # Adjust layout to reserve space at the bottom for the caption.
    plt.tight_layout(rect=[0, 0.15, 1, 1])

    # Define the caption text to be embedded in the figure.
    caption = (
        "Figure 3. Summary of BLASTX Hits for Nitrogen Fixation Genes. This scatter plot displays, for each query sequence, "
        "the number of BLASTX hits (x-axis) versus the maximum bitscore (y-axis). Queries with a maximum bitscore ≥ 40 "
        "(considered 'good' hits) are shown in forest green, while those with lower scores are in light gray. A horizontal "
        "dashed red line at a bitscore of 40 marks the threshold for statistically significant alignments, highlighting the "
        "variation in alignment quality across the dataset."
    )
    # Add the caption at the bottom center of the figure.
    plt.figtext(0.5, 0.02, caption, wrap=True, horizontalalignment="center", fontsize=10,
                bbox=dict(facecolor='white', alpha=0.5))

    # Create output directory and save the figure
    output_dir = os.path.join("blastx3-5", "blastx nitrogen fixation")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, args.output)
    plt.savefig(output_path, dpi=300)
    print("Graph saved to", output_path)
    plt.show()

if __name__ == "__main__":
    main()

