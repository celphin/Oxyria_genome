#########################
# Running HyPhy programs
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# Steps
# Run RELAX, BUSTED and absrel on guidance aligned data

######################
# Narval2
tmux new-session -s total
tmux attach-session -t total

# Load modules
module load  StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 hyphy/2.5.49

##########################

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
tmux new-session -s total1
tmux attach-session -t total1

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

mv *_HyPhy_nonArctic.txt ./Arctic_trees

################################
# Rerun HyPhy with guidance prank alignments

# get file list of non empty files

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

mkdir -p logs_nonArctic

for msa in *_guidance_pal2nal.fasta; do

    # get OG prefix
    og=${msa%_guidance_pal2nal.fasta}

    tree="${og}_tree.txt_HyPhy_nonArctic.txt"

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
        echo "$og"
    fi

done > filtered_guidance_tree_list_nonArctic.txt

wc -l filtered_guidance_tree_list_nonArctic.txt
# 10250 filtered_guidance_tree_list_nonArctic.txt

##############################
# Test ABSREL, RELAX and BUSTED

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

msa="OG0000042_guidance_pal2nal.fasta"
tree="OG0000042_tree.txt_HyPhy_nonArctic.txt"
output="OG0000042_guidance_unique_nonArctic.nxh"

hyphy remove-duplicates.bf \
	--msa "$msa" \
	--tree "$tree" \
	--output "$output"

hyphy absrel \
	--alignment "$output" \
	--branches FOREGROUND

# Must be more than one branch for these tests?
hyphy relax \
	--alignment "$output" \
	--test FOREGROUND

hyphy busted \
	--alignment "$output"
	--branches FOREGROUND

#-------------------------
# use nano to import text
cat << 'EOF' > absrel_guidance_nonArctic_array.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=logs_nonArctic/nonArctic_HyPhy_%A_%a.out
#SBATCH --error=logs_nonArctic/nonArctic_HyPhy_%A_%a.err

set +e

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

CHUNK=20

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

TOTAL=$(wc -l < filtered_guidance_tree_list_nonArctic.txt)

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

sed -n "${START},${END}p" filtered_guidance_tree_list_nonArctic.txt | while read f
do
    echo "Processing $f"

    msa="${f}_guidance_pal2nal.fasta"
    tree="${f}_tree.txt_HyPhy_nonArctic.txt"
    output="${f}_guidance_unique_nonArctic.nxh"
    final="${f}_guidance_unique_nonArctic.nxh.ABSREL.json"
    clean_tree="${f}_tree.txt_HyPhy_nonArctic_relax.txt"

    if [[ -f "$final" ]]; then
        echo "Output exists for $f — skipping"
        continue
    fi

    # Remove duplicate sequences
    hyphy remove-duplicates.bf \
        --msa "$msa" \
        --tree "$tree" \
        --output "$output"

    # Run ABSREL
    hyphy absrel \
        --alignment "$output" \
        --branches FOREGROUND

    # Run BUSTED
    hyphy busted \
        --alignment "$output" \
        --branches Leaves

    # Count number of foreground branches
    FG_COUNT=$(grep -o "Foreground" "$tree" | wc -l)

    if [[ "$FG_COUNT" -ge 2 ]]; then
        echo "Running RELAX (FG=$FG_COUNT)"

		# RELAX
		sed -E 's/([,(])([0-9]+):/\1N_\2:/g' "$tree" > "$clean_tree"

		# Remove duplicate sequences
		hyphy remove-duplicates.bf \
			--msa "$msa" \
			--tree "$clean_tree" \
			--output "$output"

		hyphy relax \
			--alignment "$output" \
			--test FOREGROUND

		STATUS=$?

		if [ $STATUS -ne 0 ]; then
			echo "RELAX failed for $f (exit code $STATUS)" >> relax_failures.log
			continue
		fi

    else
        echo "Skipping RELAX (FG=$FG_COUNT)"
    fi


    echo "Finished $f"
done

EOF

chmod +x absrel_guidance_nonArctic_array.sh
dos2unix absrel_guidance_nonArctic_array.sh

sbatch --array=1-566%100 absrel_guidance_nonArctic_array.sh


#--------------------
# Check
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *_guidance_unique_nonArctic.nxh | wc -l
ls *_guidance_unique_nonArctic.nxh.ABSREL.json  | wc -l
ls *_guidance_unique_nonArctic.nxh.BUSTED.json  | wc -l
ls *_guidance_unique_nonArctic.nxh.RELAX.json  | wc -l

# 9238
# 8926
# 9182
# 7648

# 10250 filtered_guidance_tree_list_nonArctic.txt

more logs_nonArctic/absrel_guidance_57983652_1.out

#########################################