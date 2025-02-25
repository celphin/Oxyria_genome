###################################################################################
#Converting gff3 files to cds by mapping to reference
#Needs: gff3_to_cds.sh
    #To run: sbatch gff3_to_cds.sh annotation.gff3 assembly.fasta outputannotation.fasta
#Running on:
    #Oxyria digyna (our version AED=0.6, 1.0)
    #Rhuem nobile
    #Polygonun aviculare
    #Fagopyrum tataricum (AED-0.6, 1.0)
    #Oxyria digyna (NCBI version: AED=0.6, 1.0)
###################################################################################

sbatch gff3_to_cds.sh FINAL_Oxyria.AED_0.6.sorted.gff3   Oxyria_Main.fasta Oxyria.AED_0.6.genes.fasta
sbatch gff3_to_cds.sh FINAL_Oxyria.AED_1.sorted.gff3   Oxyria_Main.fasta Oxyria.AED_1.genes.fasta
sbatch gff3_to_cds.sh R_nobile.gff3   R_nobile.fasta R_nobile.genes.cds.fasta
sbatch gff3_to_cds.sh FINAL_Polavi.AED_0.6.sorted.gff3  Polavi_Main.fasta Pol-avi.AED_0.6.genes.fasta
sbatch gff3_to_cds.sh FINAL_Polavi.AED_1.sorted.gff3  Polavi_Main.fasta Pol-avi.AED_1.genes.fasta
sbatch gff3_to_cds.sh FINAL_Fagopyrum.AED_0.6.sorted.gff3   Fagopyrum_Main.fasta Fagopyrum.AED_0.6.genes.fasta
sbatch gff3_to_cds.sh FINAL_Fagopyrum.AED_1.sorted.gff3   Fagopyrum_Main.fasta Fagopyrum.AED_1.genes.fasta
sbatch gff3_to_cds.sh FINAL_Oxyria_NCBI.AED_0.6.sorted.gff3   Oxyria_NCBI_Main.fasta Oxyria_NCBI.AED_0.6.genes.fasta
sbatch gff3_to_cds.sh FINAL_Oxyria_NCBI.AED_1.sorted.gff3   Oxyria_NCBI_Main.fasta Oxyria.AED_1.genes.fasta

