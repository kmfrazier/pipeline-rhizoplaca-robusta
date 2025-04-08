#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import os

# Define input and output file paths
input_path = "/grphome/grp_lichenscapstone/blastx3-5/blastx_nitrogenfixation/annotated_blastx_results.csv"
output_path = "/grphome/grp_lichenscapstone/blastx3-5/blastx_nitrogenfixation/blastx_hits_per_organism.png"

# Load the annotated BLASTX results CSV file
df = pd.read_csv(input_path)

# Count the number of hits per organism
organism_counts = df["Organism"].value_counts().sort_values(ascending=False)

# Create the bar graph
plt.figure(figsize=(12, 6))
bars = plt.bar(organism_counts.index, organism_counts.values, color="skyblue", edgecolor="black", label="Hit Count")
plt.xlabel("Organism", fontsize=12)
plt.ylabel("Number of Hits", fontsize=12)
plt.title("Number of BLASTX Hits per Organism Involved in Nitrogen Fixation", fontsize=14)
plt.xticks(rotation=90, ha="right")
plt.legend(loc="upper right")

# Reserve space at the bottom for a caption and embed the caption into the figure
caption = ("Figure 3. Bar graph showing the number of BLASTX hits per organism involved in nitrogen fixation. "
           "Organisms with higher hit counts indicate a greater contribution to the nitrogen fixation process, "
           "revealing the distribution of functional genes among symbiotic partners.")
plt.figtext(0.5, 0.01, caption, wrap=True, horizontalalignment="center", fontsize=10,
            bbox=dict(facecolor="white", alpha=0.5, edgecolor="none"))

plt.tight_layout(rect=[0, 0.07, 1, 1])  # Adjust layout to leave room for the caption

# Save the figure at high resolution
plt.savefig(output_path, dpi=300)
print("Bar graph saved to:", output_path)
plt.show()

