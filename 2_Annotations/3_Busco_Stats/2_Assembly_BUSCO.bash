##############################################################################
#Running BUSCO version 5.1.2 on the following assemblies:
    #Oxyria digyna (our version)
    #Oxyria digyna (NCBI version)
    #Polygonum aviculare
    #Fagopyrum escelentum (H1, H2)
    #Fagopyrum tataricum (H1, H2, NCBI version)
    #Rheum Nobile
    #Rheum tangaticum
    #Rumex hastatulus
##############################################################################
tmux new-session -s BUSCO
salloc -c10 --time 9:55:00 --mem 120000m --account def-rieseber
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate
conda deactivate

#Oxyria digyna (our version):
busco --offline --in  Oxyria_digyna.fasta  \
--out  BUSCO_Oxyria_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Oxyria_NCBI:
busco --offline --in  Oxyria_ragtag.fasta  \
--out  BUSCO_Oxyria_NCBI_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Polygonum aviculare:
busco --offline --in Polavi_Main.fasta  \
--out  BUSCO_Polavi_Main_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum escelentum H1
busco --offline --in F_escelentum_H1.fasta  \
--out  BUSCO_F_escelentum_H1_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum escelentum H2
busco --offline --in F_escelentum_H2.fasta  \
--out  BUSCO_F_escelentum_H2_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum tataricum H1
busco --offline --in F_tataricum_H1.fasta  \
--out  BUSCO_F_tataricum_H1_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum tataricum H2
busco --offline --in F_tataricum_H2.fasta  \
--out  BUSCO_F_tataricum_H2_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum tataricum (NCBI)
busco --offline --in Fagopyrum.fasta  \
--out  BUSCO_Fagopyrum_Main_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Rheum nobile:
busco --offline --in R_nobile.fasta  \
--out  BUSCO_R_nobile_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Rheum tangaticum:
busco --offline --in R_tangaticum.fasta  \
--out  BUSCO_R_tangaticum_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Rumex hastatulus
busco --offline --in R_hastatulus.fasta  \
--out  BUSCO_R_hastatulus_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
