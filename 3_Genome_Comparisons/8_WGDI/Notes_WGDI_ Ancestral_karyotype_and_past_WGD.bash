#####################################
# WGDI - Ancestral karyotype adn past whole genome duplications
# Oct 2024
# https://github.com/SunPengChuan/wgdi
##############################

# On Cedar
tmux new-session -s wgdi
tmux attach-session -t wgdi

cd /home/celphin/scratch/Oxyria/wgdi


# Install
module load StdEnv/2020  python/3.9.6 scipy-stack/2021a 

pip install wgdi
# Successfully installed wgdi-0.6.5

# or can try
#git clone https://github.com/SunPengChuan/wgdi.git
#cd wgdi
#python setup.py install

#-------------------------
# dependancies

module load StdEnv/2020  python/3.9.6 scipy-stack/2021a ipykernel/2021a gcc/9.3.0 mafft/7.471 paml/4.9j muscle/3.8.1551 trimal/1.4 fasttree/2.1.11  diamond/2.1.6


# other programs not as modules
# https://www.bork.embl.de/pal2nal/#Download
wget https://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -xzvf pal2nal.v14.tar.gz
cd pal2nal.v14

# http://www.iqtree.org/#download
wget https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz
tar -xzvf iqtree-2.3.6-Linux-intel.tar.gz
cd iqtree-2.3.6-Linux-intel

# https://github.com/simonwhelan/Divvier
git clone https://github.com/simonwhelan/Divvier.git
cd Divvier
make

# check PATH of modules
$PATH 

# set paths to dependancies
wgdi -conf help > conf.ini

nano conf.ini
[ini]
mafft_path = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/mafft/7.471/bin/mafft
pal2nal_path = /home/celphin/scratch/Oxyria/wgdi/pal2nal.v14/pal2nal.pl
yn00_path = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/paml/4.9j/bin/yn00
muscle_path = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/muscle/3.8.1551/bin/muscle
iqtree_path = /home/celphin/scratch/Oxyria/wgdi/iqtree-2.3.6-Linux-intel/bin/iqtree2
trimal_path = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trimal/1.4/bin/trimal
fasttree_path = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/fasttree/2.1.11/bin/FastTree
divvier_path = /home/celphin/scratch/Oxyria/wgdi/Divvier/divvier

# finalize
wgdi -conf conf.ini

###################################
# Usage
# https://wgdi.readthedocs.io/en/latest/

# initial setup - make total file
wgdi -conf ? > total.conf

# gff1 , lens1 , genome1_name and gff2, lens2, genome2_name 
# represent the files of species 1 and 2 respectively

# get the configuration file in the current directory and modify to run
wgdi -* ? > *.conf

#----------------------------
# Make lens and gff3 required files
# https://github.com/ABDULLAHGUJJAR/WGDI-Tool-Installation-Tutorial

mkdir prep_files; cd prep_files
wget https://raw.githubusercontent.com/ABDULLAHGUJJAR/WGDI-Tool-Installation-Tutorial/refs/heads/main/01.getgff.py
wget https://raw.githubusercontent.com/ABDULLAHGUJJAR/WGDI-Tool-Installation-Tutorial/refs/heads/main/02.gff_lens.py
wget https://raw.githubusercontent.com/ABDULLAHGUJJAR/WGDI-Tool-Installation-Tutorial/refs/heads/main/03.seq_newname.py
wget https://raw.githubusercontent.com/ABDULLAHGUJJAR/WGDI-Tool-Installation-Tutorial/refs/heads/main/rundiamond.py

chmod +755 *

# copy over the gff3 and protein files from GeneSpace to Cedar by globus
# Beluga - /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/

#------------------------
# make CDS files from reference and gff3
cd /home/celphin/scratch/Oxyria/wgdi/prep_files

tmux new-session -s wgdi
tmux attach-session -t wgdi

# formatting
module load StdEnv/2020  python/3.9.6 scipy-stack/2021a ipykernel/2021a gcc/9.3.0 mafft/7.471 paml/4.9j muscle/3.8.1551 trimal/1.4 fasttree/2.1.11  diamond/2.1.6

python

import sys
import pandas as pd
data = pd.read_csv("./genomes/Malus_sylvestris/Malus_sylvestris.gff", sep="\t", comment='#', header=None)
data = data[data[2] == 'mRNA']
data = data.loc[:, [0, 8, 3, 4, 6]]
data[8] = data[8].str.split(':|=|;|;',expand=True)[1]
# data.drop_duplicates(subset=[8], keep='first', inplace=True)
data[0] = data[0].str.replace('PAV_r1.0chr','')
data.to_csv(sys.argv[2], sep="\t", header=None, index=False)





#------------
nano 01.getgff.py

import sys
import pandas as pd
data = pd.read_csv(sys.argv[1], sep="\t", comment='#', header=None)
data = data[data[2] == 'mRNA']
data = data.loc[:, [0, 8, 3, 4, 6]]
data[8] = data[8].str.split(':|=|;|;',expand=True)[1]
# data.drop_duplicates(subset=[8], keep='first', inplace=True)
data[0] = data[0].str.replace('PAV_r1.0chr','')
data.to_csv(sys.argv[2], sep="\t", header=None, index=False)
# add to all 3 comment='#',
#-----------
python 03.seq_newname.py

import sys
import pandas as pd
from Bio import SeqIO
data = pd.read_csv(sys.argv[1], sep="\t", comment='#', header=None, index_col=6)
id_dict = data[1].to_dict()
print(data.head())
seqs = []
n = 0
for seq_record in SeqIO.parse(sys.argv[2], "fa"):
        if seq_record.id in id_dict:
                seq_record.id = id_dict[seq_record.id]
                n += 1
        else:
                continue
        seqs.append(seq_record)
SeqIO.write(seqs, sys.argv[3], "fasta")
print(n)

#----------------------

spp_hap1=Malus_sylvestris
python 01.getgff.py ./genomes/${spp_hap1}/${spp_hap1}.gff ./${spp_hap1}.bed
python 02.gff_lens.py ./bed/${spp_hap1}.bed ${spp_hap1} ./genomes/${spp_hap1}/${spp_hap1}.gff ${spp_hap1}.lens
#python 03.seq_newname.py ${spp_hap1}.gff ${spp_hap1}.cds ${spp_hap1}.cds.fasta 
python 03.seq_newname.py ./genomes/${spp_hap1}/${spp_hap1}.gff ./peptide/${spp_hap1}.fa ${spp_hap1}.pep.fasta
# pep file is empty

spp_hap2=Dryas_octopetala
python 01.getgff.py ./genomes/${spp_hap2}/${spp_hap2}.gff3 ./${spp_hap2}.bed
python 02.gff_lens.py ${spp_hap2}.bed ${spp_hap2} ./genomes/${spp_hap2}/${spp_hap2}.gff3 ${spp_hap2}.lens
#python 03.seq_newname.py ${spp_hap1}.gff ${spp_hap1}.cds ${spp_hap1}.cds.fasta 
python 03.seq_newname.py ./genomes/${spp_hap2}/${spp_hap2}.gff3 ./peptide/${spp_hap2}.fa ${spp_hap2}.pep.fasta

# amybe tyr to use the files from GENESPACE more

# blast
# https://github.com/bbuchfink/diamond/wiki
#diamond makedb --in ${spp_hap1}.pep.fasta  -d db_${spp_hap1}
python rundiamond.py ${spp_hap1}.pep.fasta ${spp_hap2}.pep.fasta db_${spp_hap1}_${spp_hap2} ${spp_hap1}_${spp_hap2}.blast
# Error opening file db_Malus_sylvestris_Dryas_octopetaladia_db: No such file or directory


#----------------------------
# make a dotplot
# Try Rosaceae compare - Dryas and Prunus or Malus
# Try Polygonaceae compare - Oxyria and Rhuem
# https://wgdi.readthedocs.io/en/latest/dotplot.html
wgdi -d help >> total.conf

# edit file for my spp
[dotplot]
   blast = ../blast/pyu_cor.blast 
   gff1 =  ../pyu.gff
   gff2 =  ../../cor/cor.gff
   lens1 = ../pyu.lens
   lens2 = ../../cor/cor.lens
   genome1_name =  pyu
   genome2_name =  cor
   multiple  = 2
   score = 100
   evalue = 1e-5
   repeat_number = 2
   position = order
   blast_reverse = false
   ancestor_left = none
   ancestor_top = none
   markersize = 1
   figsize = 10,10
   savefig = pyu_cor.dotplot.order.local.pdf

# to run
wgdi -d total.conf

#-------------------------------
# https://wgdi.readthedocs.io/en/latest/collinearity.html
wgdi -icl ? >> total.conf

# edit file
[collinearity]
gff1 = gff1 file
gff2 = gff2 file
lens1 = lens1 file
lens2 = lens2 file
blast = blast file
blast_reverse = false
multiple  = 1
process = 8
evalue = 1e-5
score = 100
grading = 50,40,25
mg = 40,40
pvalue = 0.2
repeat_number = 10
positon = order
savefile = collinearity file

# to run
wgdi -icl total.conf

#---------------------------
# https://wgdi.readthedocs.io/en/latest/ks.html

wgdi -ks ? >> total.conf

[ks]
cds_file = cds file
pep_file = pep file
align software = muscle
pairs_file = gene  pairs file
ks_file = ks result

wgdi -ks total.conf

#------------------------------





