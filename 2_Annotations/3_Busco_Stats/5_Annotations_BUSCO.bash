############################################################
#Running BUSCO version 5.1.2 on the following annotations:
    #Oxyria digyna (our version: AED=0.6, AED=1.0)
    #Oxyria digyna (NCBI: AED=0.6, AED=1.0)
    #Polygonum aviculare
    #Fagopyrum escelentum (H1, H2)
    #Fagopyrum tataricum (H1, H2, NCBI AED=0.6, 1.0)
    #Rheum nobile
    #Rheum tangaticum
############################################################
tmux new-session -s BUSCO
salloc -c10 --time 9:55:00 --mem 120000m --account def-rieseber
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate
conda deactivate

#Oxyria digyna AED=0.6
busco --offline --in Oxyria.AED_0.6.genes.fasta  \
--out  BUSCO_Oxyria_AED0.6_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Oxyria digyna AED=1.0
busco --offline --in Oxyria.AED_1.genes.fasta  \
--out  BUSCO_Oxyria_AED1_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/


#Oxyria NCBI AED=0.6:
busco --offline --in Oxyria_NCBI.AED_0.6.genes.fasta  \
--out  BUSCO_Oxyria_NCBI_AED0.6_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/


#Oxyria NCBI AED=1.0
busco --offline --in Oxyria.AED_1.genes.fasta  \
--out  BUSCO_Oxyria_NCBI_AED1_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Polygonum aviculare AED=0.6
busco --offline --in Pol-avi.AED_0.6.genes.fasta  \
--out  BUSCO_Polavi_AED0.6_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Polygonum aviculare AED1.0
busco --offline --in Pol-avi.AED_1.genes.fasta  \
--out  BUSCO_Polavi_AED1_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum escelentum (H1)
busco --offline --in F_escelentum_H1.genes.cds.fasta  \
--out  BUSCO_F_escelentum_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum escelentum (H2)
busco --offline --in F_escelentum_H2.genes.cds.fasta  \
--out  BUSCO_F_escelentum_H2_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/


#Fagopyrum tataricum (H1)
busco --offline --in F_tataricum_H1.genes.cds.fasta  \
--out  BUSCO_F_tataricum__H1_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum tataricum (H2)
busco --offline --in F_tataricum_H2.genes.cds.fasta  \
--out  BUSCO_F_tataricum_H2_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum tataricum NCBI AED0.6
busco --offline --in Fagopyrum.AED_0.6.genes.fasta  \
--out  BUSCO_Fagopyrum_AED0.6_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Fagopyrum tataricum NCBI AED1.0
busco --offline --in Fagopyrum.AED_1.genes.fasta  \
--out  BUSCO_Fagopyrum_AED1_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/


#Rheum nobile
busco --offline --in R_nobile.genes.cds.fasta  \
--out  BUSCO_R_nobile_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

#Rheum tangaticum
busco --offline --in R_tangaticum.genes.cds.fasta  \
--out  BUSCO_R_tangaticum_Annotation_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

############################################################################
tmux new-session -s BUSCO
salloc -c10 --time 9:55:00 --mem 120000m --account def-rieseber
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate

#Protein mode:
#Oxyria NCBI AED=0.6:
busco --offline --in Oxyria_proteins_AED0.6.fasta  \
--out  BUSCO_Oxyria_MS_CE_AED0.6_Proteins_eudicots  --lineage_dataset eudicots_odb10 --mode protein --cpu 20 \
--download_path ~/BUSCO_downloads/

#Same output with augustus 

busco --offline --in Fagopyrum_proteins_AED0.6.fasta   \
--out  BUSCO_Fagopyrum_MS_CE_AED0.6_Proteins_eudicots  --lineage_dataset eudicots_odb10 --mode protein --cpu 20 \
--download_path ~/BUSCO_downloads/

busco --offline --in Polavi_proteins_AED0.6.fasta  \
--out  BUSCO_Polavi_AED0.6_Proteins_eudicots  --lineage_dataset eudicots_odb10 --mode protein --cpu 20 \
--download_path ~/BUSCO_downloads/

busco --offline --in Oxyria_NCBI_proteins_AED0.6.fasta  \
--out  BUSCO_Oxyria_NCBI_AED0.6_Proteins_eudicots  --lineage_dataset eudicots_odb10 --mode protein --cpu 20 \
--download_path ~/BUSCO_downloads/




busco --offline --in F_escelentum_H1.proteins.fasta  \
--out  BUSCO_F_escelentum_AED0.6_Proteins_eudicots  --lineage_dataset eudicots_odb10 --mode protein --cpu 20 \
--download_path ~/BUSCO_downloads/

