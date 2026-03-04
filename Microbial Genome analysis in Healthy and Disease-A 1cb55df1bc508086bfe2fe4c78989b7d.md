# Microbial Genome analysis in Healthy and Disease-Affected  Human Oral Environment:

The steps  and their associated codes :

This study analyzed 30 metagenome-assembled genomes (MAGs) obtained from sequencing dental plaque samples of 30 individuals. Metadata associated with each MAG included age, sex, body mass index, study group (*healthy, mucositis, peri-implantitis*), and smoking status (*non-smoker, smoker, ex-smoker*).

# 1- quality control:

CheckM (v. 1.2.3.): we have folder (”mag”)contains 30 mags

```bash
(base) $ conda activate checkm
(checkm) $ checkm -h
(checkm) $ mkdir checkm_output
(checkm) $ bunzip2 mags/*
(checkm) $ checkm  taxonomy_wf domain Bacteria mags checkm_output -t 4
# t4 to make it faster 

```

 for the analysis and making plots out of the CheckM : by using the output files of CheckM 

```python
# checkm
# QC of the MAGs

import os
import pandas as pd 

os.chdir('path/checkm_output/storage')

# setting display option to show the whole DataFrame
pd.set_option('display.max_rows', None)  # Remove row limit
pd.set_option('display.max_columns', None)  # Remove column limit
pd.set_option('display.width', None)  # Allow unlimited width
pd.set_option('display.max_colwidth', None)  # Show full column content

# reading bin_stats.analyze.tsv
check_analyze = pd.read_csv("bin_stats.analyze.tsv", sep="\t", header=None)

# string -> real dictionary -> df -> concatenating two dfs
check_dicts = [eval(x) for x in check_analyze[1]]  
check_expanded = pd.DataFrame(check_dicts) 
check_final = pd.concat([check_analyze[0], check_expanded], axis=1)
check_final.rename(columns={0: "MAG_ID"}, inplace=True)  # final stat df

# reading bin_stats_ext.tsv
check_ext = pd.read_csv("bin_stats_ext.tsv", sep="\t", header=None, dtype=str)
check_ext_dict = [eval(x.replace("null", "None")) for x in check_ext[1]]
check_ext_expanded = pd.DataFrame(check_ext_dict)
check_ext_final = pd.concat([check_ext[0], check_ext_expanded], axis=1)
check_ext_final.rename(columns={0: 'MAG_ID'}, inplace=True)  # final extended df

# reading marker_gene_stats.tsv
check_marker = pd.read_csv("marker_gene_stats.tsv", sep="\t", header=None)
check_marker.iloc[:, 1] = check_marker.iloc[:, 1].apply(eval)
rows = []
for idx, row in check_marker.iterrows():
    mag_id = row[0]  # the first column: MAG id
    node_data = row[1]  # the second column: dictionary for nodes and genes
    
    # iterate through each node in the dictionary
    for node, gene_data in node_data.items():
        # for each gene in the node
        for gene, positions in gene_data.items():
            # for each position list associated with the gene
            for pos in positions:
                # ensure pos is a list and has at least one element
                start = pos[0] if isinstance(pos, list) and len(pos) > 0 else None
                end = pos[1] if isinstance(pos, list) and len(pos) > 1 else None
                rows.append([mag_id, node, gene, start, end])

check_marker_final = pd.DataFrame(rows, columns=["mag_id", "node", "gene_ID", "start", "end"])  # final marker df
```

```python
# completeness & contamination joint df

size_complete = pd.merge(
    check_ext_final[['MAG_ID', 'Completeness', 'Contamination']], 
    check_final[['MAG_ID', 'GC', 'GC std', 'Genome size']], 
    on='MAG_ID',  # Merge on the MAG_ID column
    how='outer'  # Keeps all MAGs, even if missing in one DF
)

# checking if the merge is right and corresponds to the correct MAGS forom 2 different dfs
merged_check = size_complete.merge(
    check_ext_final[['MAG_ID', 'Completeness', 'Contamination']], on="MAG_ID", suffixes=('', '_original_ext')
).merge(
    check_final[['MAG_ID', 'GC', 'GC std', 'Genome size']], on="MAG_ID", suffixes=('', '_original_final')
)

print((merged_check['Completeness'] == merged_check['Completeness_original_ext']).all())  # Should be True
print((merged_check['Contamination'] == merged_check['Contamination_original_ext']).all())  # Should be True
print((merged_check['GC'] == merged_check['GC_original_final']).all())  # Should be True

print(size_complete)
```

```python
# merging size_complete df and metadata

merged_meta = size_complete.merge(
    metadata[['magID', 'study_group', 'smoking_state']], 
    left_on='MAG_ID', right_on='magID', how='left'
)

# dropping the extra 'magID' column after merge
merged_meta.drop(columns='magID', inplace=True)
columns_order = ['MAG_ID','study_group', 'smoking_state', 'Genome size'] + [col for col in merged_meta.columns if col not in ['MAG_ID', 'study_group', 'smoking_state', 'Genome size']]
merged_meta = merged_meta[columns_order]
```

```python
# plotting completeness of MAGs against genome size
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# fitting regression model using statsmodels to get predicted values
X = sm.add_constant(merged_meta['Completeness'])  # Add constant for intercept
y = merged_meta['Genome size']
model = sm.OLS(y, X).fit()
predicted_values = model.predict(X)

# calculating residuals (difference between actual and predicted values)
residuals = y - predicted_values

# calculating the standard deviation of the residuals
residual_std = residuals.std()
# setting threshold as 2 standard deviations
threshold = 1.96 * residual_std

sns.set_palette(["#FFFACD", "#FAFAD2", "#FFEFD5"])

# creating scatter plot with regression line
plt.figure(figsize=(10, 6))
sns.regplot(x="Completeness", y="Genome size", data=merged_meta,
            scatter_kws={'alpha':0.6, 'color': 'green'},
            line_kws={'color': '#FFC107'})

# adding MAG names to points that are far from the regression line
for i in range(merged_meta.shape[0]):
    if abs(residuals.iloc[i]) > threshold:
        plt.text(merged_meta['Completeness'].iloc[i], merged_meta['Genome size'].iloc[i], 
                 merged_meta['MAG_ID'].iloc[i], fontsize=7, ha='center', va='bottom')

# Titles and labels
plt.title("Completeness vs Genome Size")
plt.xlabel("Completeness (%)")
plt.ylabel("Genome Size (bp)")

plt.show()
plt.close()
```

![image.png](image.png)

The rest of the graph are in the supplementary files: 

# 2- Taxonomic Analysis

Taxonomic classification of MAGs was carried out using Phylophlan metagenomic (v. 3.0.36). MAGs were assigned to species-level genome bins (SGBs) based on the computed average Mash distance between each MAG and reference genomes. MAGs with average Mash distance below 0.05 to a reference genome were assigned to the corresponding SGB.

```bash
# Create output directory
mkdir -p phylophlan_output

# Activate the PhyloPhlAn environment which save before as "ppa"
conda activate ppa

# Define variables for directories
INPUT_DIR="mags"
OUTPUT_DIR="phylophlan_output/ppa_m"
DATABASE_DIR="ppa_db"

# Run PhyloPhlAn metagenomic analysis
phylophlan_metagenomic \
    -i "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    --nproc 4 \
    -n 1 \  # Setting -n to 1 to get the closest match
    -d CMG2425 \  # Specify the database
    --database_folder "$DATABASE_DIR" \
    --verbose

# Move into the output directory
cd "$OUTPUT_DIR"

```

After running the code we have got that all our mags contain one species which is Prevotella_koreensis (P.koreensis).

# 3- Genome Annotation:

Prokka (v.1.14.6) was employed to annotate the genomic features of the contigs within bins. For each MAG, multiple output files were generated, detailing genomic features such as coding sequences (CDS), non-coding transcripts, repetitive elements, transfer RNAs (tRNA), ribosomal RNAs (rRNA), and predicted proteins.

```bash
#!/bin/bash
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate prokka
#prokka -h
mkdir prokka_output
#run prokka several times using a for loop
for f in mags/*; do
mag=$(basename $f .fna)
#echo $mag
#echo $f
 mkdir prokka_output/${mag}
prokka mags/${mag}.fna
\
 --outdir prokka_output/${mag}
\
 --prefix ${mag}
\
 --compliant
done
#this should take ~30 mins
conda deactivate
```

Analysis of Prokka Output: Gene, CDS, and Hypothetical Protein Counts: 

```bash
#for each folder of the prokka outputs we run this code 
for file in *.tsv; do
    echo "File: $file"
    echo "Genes: $(grep -c 'gene' "$file")"
    echo "CDS: $(grep -c 'CDS' "$file")"
    echo "Hypothetical proteins: $(grep -c 'hypothetical protein' "$file")"
    echo "------------------------"
done
#To extract protein-coding genes (excluding hypothetical proteins), the following command was used awk -F'\t' '$2 == "CDS" && $7 != "hypothetical protein" {print $4}' M1025245651.tsv | sort | uniq
 $ awk -F'\t' '$2 == "CDS" && $7 != "hypothetical protein" {print $4}' magID.tsv | sort | uniq
#The gene list is included in the supplementary files.
```

![image.png](image%201.png)

In order to visualize the previous results: 

```python
import matplotlib.pyplot as plt

# Data from the table
samples = [
    "M1025245651", "M1213021579", "M1365618114", "M1486539787", "M1671123810",
    "M1709480479", "M1825178389", "M1927030333", "M1048336460", "M1247523836",
    "M1402687461", "M1548705285", "M1702705864", "M1762225173", "M1831926750",
    "M1968198586", "M1129502254", "M1314714359", "M1419890899", "M1607287943",
    "M1703261764", "M1787348583", "M1846526184", "M1140039743", "M1320943765",
    "M1481385194", "M1650726873", "M1708806623", "M1801750263", "M1846526184",
    "M1968198586"
]

genes = [1604, 1638, 2301, 1438, 2098, 1884, 1657, 2085, 2034, 1617, 
         2463, 2088, 2180, 1278, 1637, 1917, 1611, 1931, 2098, 2003, 
         1919, 2281, 2099, 2083, 1882, 2140, 2165, 1527, 1834, 2099, 1917]

hypothetical_proteins = [923, 953, 1330, 834, 1155, 1025, 898, 1118, 1088, 908, 
                         1500, 1095, 1200, 844, 914, 1039, 882, 1031, 1099, 1057, 
                         1068, 1310, 1107, 1111, 1010, 1181, 1201, 870, 1025, 1107, 1039]

# Create a scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(genes, hypothetical_proteins, color='g', alpha=0.6, edgecolors='k')

# Highlight the outlier (M1402687461)
plt.scatter(2463, 1500, color='r', edgecolors='k', s=100, label="M1402687461")

# Labels and title
plt.xlabel("Total Genes")
plt.ylabel("Hypothetical Proteins")
plt.title("Scatter Plot of Genes vs. Hypothetical Proteins")
plt.legend()
plt.grid(True)
# Show the plot
plt.show()

```

![image.png](image%202.png)

# 4-Fisher Test and Further analysis:

Fisher test was applied to get the of the gene presence across the different samples:  

```python
import pandas as pd
from scipy.stats import fisher_exact

# Load the data
df = pd.read_csv("/content/sample_data/gene_presence_absence.csv")

# Define groups
healthy_samples = ["M1548705285", "M1702705864", "M1831926750"]
disease_samples = [
    "M1025245651", "M1213021579", "M1365618114", "M1486539787", "M1671123810",
    "M1709480479", "M1825178389", "M1927030333", "M1048336460", "M1247523836",
    "M1402687461", "M1762225173", "M1968198586", "M1129502254", "M1314714359",
    "M1419890899", "M1607287943", "M1703261764", "M1787348583", "M1846526184",
    "M1140039743", "M1320943765", "M1481385194", "M1650726873", "M1708806623",
    "M1801750263", "M1876746341"
]

# Extract gene presence/absence
samples = healthy_samples + disease_samples
presence_absence_df = df[["Gene"] + samples].copy()

# Convert to binary presence/absence
for sample in samples:
    presence_absence_df[sample] = presence_absence_df[sample].notnull().astype(int)

# Function for Fisher's exact test
def fisher_test(row):
    healthy_present = row[healthy_samples].sum()
    healthy_absent = len(healthy_samples) - healthy_present
    disease_present = row[disease_samples].sum()
    disease_absent = len(disease_samples) - disease_present
    table = [[healthy_present, healthy_absent], [disease_present, disease_absent]]
    _, p_value = fisher_exact(table, alternative='two-sided')
    return p_value

# Apply Fisher’s test
presence_absence_df["p_value"] = presence_absence_df.apply(fisher_test, axis=1)

# Get significant genes (p < 0.05)
significant_genes = presence_absence_df[presence_absence_df["p_value"] < 0.05]

# Save results
significant_genes[["Gene", "p_value"]].to_csv("significant_genes.csv", index=False)
```

 # the output file is in supplementary files

To check  the samples’ states  that  are enrich with the significant_genes

```python
import pandas as pd

# Load the data
df = pd.read_csv("/content/sample_data/gene_presence_absence.csv")

# Define groups
healthy_samples = ["M1548705285", "M1702705864", "M1831926750"]
disease_samples = [
    "M1025245651", "M1213021579", "M1365618114", "M1486539787", "M1671123810",
    "M1709480479", "M1825178389", "M1927030333", "M1048336460", "M1247523836",
    "M1402687461", "M1762225173", "M1968198586", "M1129502254", "M1314714359",
    "M1419890899", "M1607287943", "M1703261764", "M1787348583", "M1846526184",
    "M1140039743", "M1320943765", "M1481385194", "M1650726873", "M1708806623",
    "M1801750263", "M1876746341"
]

# List of significant genes
significant_gene_list = [
    "mnmE", "group_527", "group_750", "group_1897", "pqqE", "group_1906", 
    "group_2625", "group_5275", "group_1383", "group_5199", "group_1870", 
    "group_1890", "group_2686", "group_57", "group_2627", "group_2629", "group_2656"
]

# Extract presence/absence for significant genes
samples = healthy_samples + disease_samples
presence_absence_df = df[["Gene"] + samples].copy()

# Convert to binary presence/absence
for sample in samples:
    presence_absence_df[sample] = presence_absence_df[sample].notnull().astype(int)

# Filter significant genes
sig_gene_df = presence_absence_df[presence_absence_df["Gene"].isin(significant_gene_list)]

# Count presence in each group
sig_gene_df["Healthy_Present"] = sig_gene_df[healthy_samples].sum(axis=1)
sig_gene_df["Disease_Present"] = sig_gene_df[disease_samples].sum(axis=1)

# Select relevant columns
sig_gene_counts = sig_gene_df[["Gene", "Healthy_Present", "Disease_Present"]]

# Save results
sig_gene_counts.to_csv("gene_presence_counts.csv", index=False)

# Display results
print(sig_gene_counts.sort_values(by="Disease_Present", ascending=False))
import seaborn as sns
import matplotlib.pyplot as plt

# Set the plot style
sns.set(style="whitegrid")

# Prepare data for plotting
sig_gene_counts_long = sig_gene_counts.melt(id_vars="Gene", value_vars=["Healthy_Present", "Disease_Present"], 
                                            var_name="Group", value_name="Gene_Presence")

# Plot the barplot with light green color
plt.figure(figsize=(12, 8))
sns.barplot(x="Gene", y="Gene_Presence", hue="Group", data=sig_gene_counts_long, 
            palette=["lightgreen", "green"], ci=None)

# Rotate the x-axis labels for better readability
plt.xticks(rotation=90)

# Add labels and title
plt.xlabel("Gene")
plt.ylabel("Gene Presence Count")
plt.title("Gene Presence Comparison Between Healthy and Disease Samples")

# Show the plot
plt.tight_layout()
plt.show()
```

![image.png](image%203.png)

# 5- Pan-genome Analysis:

Roary (v3.13.0) was utilized in pan-genome characterization to assess the shared gene content and genetic diversity across MAGs. GFF files were used as input, from which protein sequences were extracted for the analysis to ensure more stable alignments. 

```bash
#!/bin/bash

# === Step 1: Create Conda Environment for Plotting ===
conda create -y -n roary_plots \
  r-ggplot2 \
  python=3.10.12 \
  seaborn=0.11.0 \
  matplotlib=3.8.2 \
  pandas=2.1.1 \
  numpy=1.26.0 \
  biopython=1.81

# === Step 2: Run Roary ===
# Assumes a conda environment named "roary" already exists
conda deactivate
conda activate roary

# Create output folder if not exists
mkdir -p roary_output

# Run Roary with specified parameters
roary prokka_output/*/*.gff -f roary_output -i 95 -cd 90 -p 4

# Check output
ls roary_output/

# === Step 3: Download Plotting Scripts ===
conda deactivate
conda activate roary_plots

# Download R script
curl -L https://raw.githubusercontent.com/sanger-pathogens/Roary/refs/heads/master/bin/create_pan_genome_plots.R -o create_pan_genome_plots.R

# Download Python script
curl -L https://raw.githubusercontent.com/sanger-pathogens/Roary/refs/heads/master/contrib/roary_plots/roary_plots.py -o roary_plots.py

# Move scripts to output folder
mv create_pan_genome_plots.R roary_output/
mv roary_plots.py roary_output/

# === Step 4: Generate Plots ===
cd roary_output/

# Make the R script executable and run it
chmod +x create_pan_genome_plots.R
Rscript create_pan_genome_plots.R

# Run the Python plot script
python roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv

echo " Roary analysis and plotting complete. Check the 'roary_output/' folder for results."

```

![image.png](image%204.png)

#The rest of the graphs could be found in the supplementary files

# 6- phylogenetic Analysis:

```bash
# Activate the conda environment containing Roary and FastTree
conda activate roary

# Check Roary's help menu to verify installation and see options
roary -h

# Run Roary to create a pan-genome analysis from multiple Prokka annotations
# -f: Output directory
# -cd: Minimum percentage identity for blastp (90%)
# -p: Number of CPU threads to use (4)
# -e: Create multiFASTA alignment of core genes
# -n: Disable fast core gene alignment with PRANK (slower but more accurate)
roary prokka_output/*/*.gff \
    -f roary_output_w_aln \
    -cd 90 \
    -p 4 \
    -e -n
# Note: This core genome analysis typically takes ~20 minutes with 4 threads

# Build a phylogenetic tree from the core gene alignment using FastTreeMP
# -pseudo: Handle pseudocounts for better support values
# -spr: SPR moves (4)
# -mlacc: ML accuracy (2)
# -slownni: Enable thorough NNI search
# -fastest: Combine with -slownni for optimized speed/accuracy
# -no2nd: Skip second round of NNIs
# -mlnni: ML NNI moves (4)
# -gtr: Use GTR model
# -nt: Nucleotide input
FastTreeMP -pseudo -spr 4 -mlacc 2 -slownni -fastest -no2nd -mlnni 4 -gtr -nt \
    -out roary_output_w_aln/core_gene_phylogeny.nwk \
    roary_output_w_aln/core_gene_alignment.aln
# Note: Tree construction typically takes ~10 minutes

# Deactivate the conda environment
conda deactivate

# The resulting phylogenetic tree (Newick format) can be visualized with:
# - FigTree (GUI)
# - iTOL (web-based)
# - R packages like ggtree or ape
# - Python libraries like ETE3 or DendroPy
```

![image.png](image%205.png)