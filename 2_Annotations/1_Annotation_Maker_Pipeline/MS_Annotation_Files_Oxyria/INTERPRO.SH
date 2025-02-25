sed 's/*//g' CE10.AED_0.8_protein.fasta > CE10.AED_0.8_protein.interpro.fasta
/DATA/home/jmlazaro/Software/interproscan-5.57-90.0/interproscan.sh --formats GFF3 --goterms --cpu 144 --input CE10.AED_0.8_protein.fasta --applications Pfam,SUPERFAMILY 
