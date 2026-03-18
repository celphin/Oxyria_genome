#########################
# Running HyPhy programs
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# Steps
# Run RELAX and absrel
# Explore results

######################
# Narval2
tmux attach-session -t total

# Load modules
module load  StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 hyphy/2.5.49

##########################
# Copy trees from Arctic folder back up to main to rerun with non-Arctic

tmux attach-session -t backup1
ls *_tree.txt | wc -l
#14602
ls *_tree.txt_pal2nal.fasta | wc -l
# 3993

cp OG*_tree.txt ..
cp OG*_tree.txt_pal2nal.fasta ..

cd ..
ls *_tree.txt | wc -l
ls *_tree.txt_pal2nal.fasta | wc -l

# 23145
# 23145


############################
# Check gene  names for closest non-Arctic relatives

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS
grep ">" Arabis_alpina.f*
# >AALP_AAs64319U000100
# >AALP_AAs49271U000100
# >AALP_AAs60671U000100
# >AALP_AAs48102U000100
# >AALP_AAs60109U000100
# >AALP_AAs52579U000100
# >AALP_AAs45409U000100
# >AALP_AAs64185U000100
# >AALP_AAs61717U000100
# >AALP_AAs47510U000100
# >AALP_AAs68684U000100
# >AALP_AAs73943U000100
# >AALP_AAs42227U000100

grep ">" Prunus_persica.f*
# >LOC18768756
# >LOC18767021
# >LOC18767694
# >LOC18766660
# >LOC18767366
# >LOC18768373
# >LOC18768715
# >LOC18766300
# >LOC18767169
# >LOC18767794
# >LOC18767676
# >LOC18767280
# >LOC18766388
# >LOC18768341
# >LOC18768949
# >LOC18767264
# >LOC18768104
# >LOC18767743
# >LOC18766280
# >LOC18766638
# >LOC18766399
# >LOC18766758
# >LOC18767356
# >LOC18768784
# >LOC109950649

grep ">" Rheum_nobile_H0.f*
# >RnoG0007793.1
# >RnoG0007792.1
# >RnoG0007791.1
# >RnoG0007790.1
# >RnoG0007789.1
# >RnoG0007788.1
# >RnoG0007787.1
# >RnoG0007786.1
# >RnoG0024551.1
# >RnoG0024550.1
# >RnoG0024549.1
# >RnoG0024548.1
# >RnoG0024547.1
# >RnoG0024546.1

grep ">" Arabidopsis_lyrata.f*
# >LOC110224666
# >LOC9326461
# >LOC9326463
# >LOC9326464
# >LOC9326465
# >LOC9326466
# >LOC9326467
# >LOC9329189
# >LOC9326468
# >LOC9326469
# >LOC9329190
# >LOC9326470
# >LOC9329192
# >LOC9326471
# >LOC9329193
# >LOC9329194
# >LOC9326472
# >LOC9326473
# >LOC9326474
# >LOC9326475
# >LOC9329195
# >LOC9326476

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS

grep ">" Arabis_alpina.f* | sed 's/>//' > foreground_genes.txt
grep ">" Prunus_persica.f* | sed 's/>//' >> foreground_genes.txt
grep ">" Rheum_nobile_H0.f* | sed 's/>//' >> foreground_genes.txt
grep ">" Arabidopsis_lyrata.f* | sed 's/>//' >> foreground_genes.txt


mv foreground_genes.txt /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

sort -u foreground_genes.txt > foreground_genes_unique.txt

#------------------------------------
# Add {Foreground} label to Arctic branches

# Narval1
tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

for f in *tree.txt
do
    echo "Processing $f"

    awk '
    NR==FNR {fg[$1]=1; next}

    {
        line=$0
        out=""
        pos=1

        while (match(substr(line,pos), /[^:(),]+:/)) {
            start = pos + RSTART - 1
            end = start + RLENGTH - 2
            gene = substr(line, start, RLENGTH-1)

            out = out substr(line, pos, RSTART-1)

            if (gene in fg)
                out = out gene "{Foreground}:"
            else
                out = out gene ":"

            pos = start + RLENGTH
        }

        out = out substr(line,pos)
        print out
    }
    ' foreground_genes_unique.txt "$f" > "${f}_HyPhy_nonArctic.txt"

done

ls *_HyPhy_nonArctic.txt | wc -l
# 23145

################################
# Remove new trees without Foreground branches
grep -L "{Foreground}" *_HyPhy_nonArctic.txt | xargs rm

ls *_HyPhy_nonArctic.txt | wc -l
# 18779

##################################
# loop through orthogroups that include nonArctic spp 
# run as slurm array

# get file list of non empty files
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/

mkdir -p logs

for f in *_tree.txt; do
    msa="${f}_pal2nal.fasta"
    tree="${f}_HyPhy_nonArctic.txt"

    # Require both files
    [[ -f "$msa" && -f "$tree" ]] || continue

    # Require at least 2 sequences
    seq_count=$(grep -c "^>" "$msa")
    [[ $seq_count -ge 2 ]] || continue

    # Skip zero-length sequences (ignoring gaps)
    if awk '
        /^>/ {
            if (seq_len == 0 && NR > 1) exit 1
            seq_len = 0
            next
        }
        { gsub("-", ""); seq_len += length($0) }
        END { if (seq_len == 0) exit 1 }
    ' "$msa"
    then
        echo "$f"
    fi

done > filtered_tree_list_nonArctic.txt

# check count
wc -l filtered_tree_list_nonArctic.txt

# 16401 filtered_tree_list.txt (Arctic)
# 16176 filtered_tree_list_nonArctic.txt

#----------------------
# use nano to import text
cat << EOF > absrel_array_nonArctic.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=logs/absrel_nonArc_%A_%a.out
#SBATCH --error=logs/absrel_nonArc_%A_%a.err

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/nonArctic_trees

CHUNK=20

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

TOTAL=$(wc -l < filtered_tree_list_nonArctic.txt)

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

sed -n "${START},${END}p" filtered_tree_list_nonArctic.txt | while read f
do
    echo "Processing $f"

    msa="${f}_pal2nal.fasta"
    tree="${f}_HyPhy_nonArctic.txt"
    output="${f}_unique_nonArctic.nxh"

    if [[ -f "$output" ]]; then
        echo "Output exists for $f — skipping"
        continue
    fi

    hyphy remove-duplicates.bf \
        --msa "$msa" \
        --tree "$tree" \
        --output "$output"

    hyphy absrel \
        --alignment "$output" \
        --branches FOREGROUND

    echo "Finished $f"
done

EOF

chmod +x absrel_array_nonArctic.sh
dos2unix absrel_array_nonArctic.sh

cp Arctic_trees/remove-duplicates.bf .

sbatch --array=1-809%100 absrel_array_nonArctic.sh
# Submitted batch job 57490810

###############################
# Checking results
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

ls *_unique_nonArctic.nxh |wc -l
# 16176
ls *_unique_nonArctic.nxh.ABSREL.json | wc -l 
# 16044

##########################
