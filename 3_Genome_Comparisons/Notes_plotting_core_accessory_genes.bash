############################
# Look at positions of core and accessory genes
# Feb 2025
############################

cd 

#----------------------------
# Core genes - in all
awk 'NR==1 {print; next} 
{
  present=1;
  for(i=2;i<=NF;i++){
    if($i==0){present=0}
  }
  if(present==1){print}
}' Orthogroups.GeneCount.tsv > core_orthogroups.tsv

#----------------------------
# Accessory genes - in at least one but not all (7 less than all)
awk 'NR==1 {print; next}
{
  count=0;
  for(i=2;i<=NF;i++){
    if($i>0){count++}
  }
  if(count>0 && count<NF-7){print}
}' Orthogroups.GeneCount.tsv > accessory_orthogroups.tsv

#----------------------------
# extract core gene IDs
cut -f1 core_orthogroups.tsv | tail -n +2 > core_ids.txt
grep -Ff core_ids.txt Orthogroups.tsv > core_genes.tsv

#----------------------------
# Extract core gene positions
grep -Ff core_gene_ids.txt all_genes.bed > core_genes.bed


#----------------------------
# make density windows
bedtools makewindows -g genome.fai -w 100000 > windows.bed
bedtools intersect -a windows.bed -b core_genes.bed -c > core_density.bed
# and use IGV Visualization of bed files


#----------------------------
# Make plot in python

import pandas as pd
import matplotlib.pyplot as plt

bed = pd.read_csv("core_genes.bed", sep="\t", header=None)
bed.columns = ["chr","start","end","gene"]

plt.figure(figsize=(12,4))

for i, chrom in enumerate(bed["chr"].unique()):
    subset = bed[bed["chr"]==chrom]
    plt.scatter(subset["start"], [i]*len(subset), s=5)

plt.yticks(range(len(bed["chr"].unique())), bed["chr"].unique())
plt.xlabel("Genomic position")
plt.ylabel("Chromosome")
plt.title("Core gene distribution")
plt.show()



# Or add to old plot in R


