########################################################################################
#Installation of BUSCO version 5.1.2 and according databases to Alliance Canada clusters
########################################################################################
#Install Busco
cd ~

# Setup
#module purge
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap

virtualenv ~/busco_env

source ~/busco_env/bin/activate

pip install biopython pandas busco==5.1.2 --no-index

# to close
deactivate

# to activate the environment
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate

# check installed bbmap
bbduk.sh --version

java -ea -Xmx3235m -Xms3235m -cp /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bbmap/38.86/current/ jgi.BBDuk --version
Picked up JAVA_TOOL_OPTIONS: -Xmx2g
BBMap version 38.86
For help, please run the shellscript with no parameters, or look in /docs/.

cd ~

# download data needed
mkdir BUSCO_downloads; cd BUSCO_downloads

mkdir lineages; cd lineages
wget https://busco-data.ezlab.org/v5/data/lineages/viridiplantae_odb10.2020-09-10.tar.gz
wget https://busco-data.ezlab.org/v5/data/lineages/eudicots_odb10.2020-09-10.tar.gz
wget https://busco-data.ezlab.org/v5/data/lineages/embryophyta_odb10.2020-09-10.tar.gz
wget https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz
wget https://busco-data.ezlab.org/v5/data/lineages/stramenopiles_odb10.2024-01-08.tar.gz
cd ..

mkdir information; cd information
wget https://busco-data.ezlab.org/v5/data/information/lineages_list.2021-12-14.txt.tar.gz
cd ..
mkdir placement_files; cd placement_files
wget https://busco-data.ezlab.org/v5/data/placement_files/list_of_reference_markers.eukaryota_odb10.2019-12-16.txt.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/mapping_taxid-lineage.eukaryota_odb10.2019-12-16.txt.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/supermatrix.aln.eukaryota_odb10.2019-12-16.faa.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/tree.eukaryota_odb10.2019-12-16.nwk.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/tree_metadata.eukaryota_odb10.2019-12-16.txt.tar.gz
cd..

cd lineages
tar -xvf ./viridiplantae_odb10.2020-09-10.tar.gz
tar -xvf ./eudicots_odb10.2020-09-10.tar.gz
tar -xvf ./embryophyta_odb10.2020-09-10.tar.gz
tar -xvf ./eukaryota_odb10.2024-01-08.tar.gz
tar -xvf ./stramenopiles_odb10.2024-01-08.tar.gz
cd ..
cd information
tar -xvf ./lineages_list.2021-12-14.txt.tar.gz
cd ..
cd placement_files
tar -xvf ./list_of_reference_markers.eukaryota_odb10.2019-12-16.txt.tar.gz
tar -xvf ./mapping_taxid-lineage.eukaryota_odb10.2019-12-16.txt.tar.gz
tar -xvf ./mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz
tar -xvf ./supermatrix.aln.eukaryota_odb10.2019-12-16.faa.tar.gz
tar -xvf ./tree.eukaryota_odb10.2019-12-16.nwk.tar.gz
tar -xvf ./tree_metadata.eukaryota_odb10.2019-12-16.txt.tar.gz

cd ..
