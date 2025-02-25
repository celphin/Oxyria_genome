###########################################################
# Oct 2024
# following this: https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md
# GO analysis section - plotting and explore

###########################################################
cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/
mkdir ermine.results

cp *genesets.ermine.results ermine.results
cd ermine.results

#-------------------------
# Look at the files for each set combined
# Subset those with a p-value < 0.05

grep "!" *_Arctic_specific_genesets.ermine.results  | awk -F '\t' '$7 <0.5 { print }'| awk  -F $'\t' '{print $2, $3, $7, $1}' |sort  > Arctic_specific_GO.txt

#-------------------------
grep "!" *_totalfam_all_contracted_genesets.ermine.results  | awk -F '\t' '$7 <0.5 { print }'| awk  -F $'\t' '{print $2, $3, $7, $1}' |sort  > totalfam_all_contracted_GO.txt

#-------------------------
grep "!" *_totalfam_all_expanded_genesets.ermine.results  | awk -F '\t' '$7 <0.5 { print }'| awk  -F $'\t' '{print $2, $3, $7, $1}' |sort  > totalfam_all_expanded_GO.txt

#-------------------------
grep "!" *_totalfam_rapidly_contracted_genesets.ermine.results  | awk -F '\t' '$7 <0.5 { print }'| awk  -F $'\t' '{print $2, $3, $7, $1}' |sort  > totalfam_rapidly_contracted_GO.txt

#-------------------------
grep "!" *_totalfam_rapidly_expanded_genesets.ermine.results  | awk -F '\t' '$7 <0.5 { print }'| awk  -F $'\t' '{print $2, $3, $7, $1}' |sort  > totalfam_rapidly_expanded_GO.txt


#----------------------------
# Extract just GO terms and count duplicates

cat *_Arctic_specific_genesets.ermine.results | awk -F '\t' '$7 <0.5 { print }' | awk  -F $'\t' '{print $3, $2}' | sort | uniq -c | sort | awk '$1 > 1'  | wc -l 
#85
cat *_totalfam_all_contracted_genesets.ermine.results  | awk -F '\t' '$7 <0.5 { print }' |awk  -F $'\t' '{print $3, $2}' | sort | uniq -c | sort | awk '$1 > 3' | wc -l 
#31 in all 4
cat *_totalfam_all_expanded_genesets.ermine.results | awk -F '\t' '$7 <0.5 { print }' |awk  -F $'\t' '{print $3, $2}' | sort | uniq -c | sort | awk '$1 > 3' | wc -l
#15 in all 4
cat *_totalfam_rapidly_contracted_genesets.ermine.results | awk -F '\t' '$7 <0.5 { print }' |awk  -F $'\t' '{print $3, $2}' | sort | uniq -c | sort | awk '$1 > 1' | wc -l
#9 in 2 or more
cat *_totalfam_rapidly_expanded_genesets.ermine.results | awk -F '\t' '$7 <0.5 { print }' |awk  -F $'\t' '{print $3, $2}' | sort | uniq -c | sort | awk '$1 > 1' | wc -l
#15 in 2 or more




###################################
# summarize with : http://revigo.irb.hr/

#-------------------------
cat *_Arctic_specific_genesets.ermine.results | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' |\
sort  | uniq -c | awk  -F " " '{print $2, $1}' > Revigio_Arctic_specific.txt

#-------------------------
cat *_totalfam_all_contracted_genesets.ermine.results | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' |\
sort  | uniq -c | awk  -F " " '{print $2, $1}' > Revigio_totalfam_all_contracted.txt

#-------------------------
cat *_totalfam_all_expanded_genesets.ermine.results | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' |\
sort  | uniq -c | awk  -F " " '{print $2, $1}' > Revigio_totalfam_all_expanded.txt

#-------------------------
cat *_totalfam_rapidly_contracted_genesets.ermine.results | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' |\
sort  | uniq -c | awk  -F " " '{print $2, $1}' > Revigio_totalfam_rapidly_contracted.txt

#-------------------------
cat *_totalfam_rapidly_expanded_genesets.ermine.results | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' |\
sort  | uniq -c | awk  -F " " '{print $2, $1}' > Revigio_totalfam_rapidly_expanded.txt


##################################
#-------------------------
cat *_Arctic_specific_genesets.ermine.results | awk -F '\t' '$7 <0.05 { print $3 , $7}'  > Revigio_Arctic_specific_pval.txt

#-------------------------
cat *_totalfam_all_contracted_genesets.ermine.results | awk -F '\t' '$7 <0.05   { print $3 , $7}'  > Revigio_totalfam_all_contracted_pval.txt

#-------------------------
cat *_totalfam_all_expanded_genesets.ermine.results | awk -F '\t' '$7 <0.05   { print $3 , $7}'  > Revigio_totalfam_all_expanded_pval.txt

#-------------------------
cat *_totalfam_rapidly_contracted_genesets.ermine.results | awk -F '\t' '$7 <0.05   { print $3 , $7}'  > Revigio_totalfam_rapidly_contracted_pval.txt

#-------------------------
cat *_totalfam_rapidly_expanded_genesets.ermine.results | awk -F '\t' '$7 <0.05   { print $3 , $7}'  > Revigio_totalfam_rapidly_expanded_pval.txt


#-------------------------
cat *_totalfam_all_contracted_genesets.ermine.1results | awk -F '\t' '$7 <0.05  { print $3 , $7}'  > Revigio_totalfam_all_contracted_pval1.txt

#-------------------------
cat *_totalfam_all_expanded_genesets.ermine.1results | awk -F '\t' '$7 <0.05   { print $3 , $7}'  > Revigio_totalfam_all_expanded_pval1.txt

#-------------------------
cat *_totalfam_rapidly_contracted_genesets.ermine.1results | awk -F '\t' '$7 <0.05   { print $3 , $7}'  > Revigio_totalfam_rapidly_contracted_pval1.txt

#-------------------------
cat *_totalfam_rapidly_expanded_genesets.ermine.1results | awk -F '\t' '$7 <0.05   { print $3 , $7}'  > Revigio_totalfam_rapidly_expanded_pval1.txt

###################################
# Look at all and what species they are from
cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/ermine.results

grep "!" *totalfam*_genesets.ermine.results  | awk -F '\t' '$7 <1 { print }'| awk  -F $'\t' '{print $1, $2, $3, $7}' |sort  > List_all_GO.txt

#---------------------
# Search for specific terms
grep light List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! photosynthesis, light harvesting GO:0009765 0.56966575
Cochgro_totalfam_all_expanded_genesets.ermine.results:! blue light signaling pathway GO:0009785 0.13742534
Cochgro_totalfam_all_expanded_genesets.ermine.results:! cellular response to blue light GO:0071483 0.13742534
Cochgro_totalfam_all_expanded_genesets.ermine.results:! cellular response to light stimulus GO:0071482 0.08111565
Cochgro_totalfam_all_expanded_genesets.ermine.results:! photosynthesis, light harvesting GO:0009765 0.06565836
Cochgro_totalfam_all_expanded_genesets.ermine.results:! response to blue light GO:0009637 0.3472892
Cochgro_totalfam_all_expanded_genesets.ermine.results:! response to light stimulus GO:0009416 0.66391989
Cochgro_totalfam_all_expanded_genesets.ermine.results:! response to red or far red light GO:0009639 0.6887685
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! cellular response to light stimulus GO:0071482 0.40207568
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! response to blue light GO:0009637 0.44443905
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! response to light stimulus GO:0009416 0.47410388
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! response to red or far red light GO:0009639 0.20698667
Oxydig_totalfam_all_expanded_genesets.ermine.results:! photosynthesis, light harvesting GO:0009765 0.63790537
Oxydig_totalfam_all_expanded_genesets.ermine.results:! response to light stimulus GO:0009416 0.51315873

grep heat List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! response to heat GO:0009408 0.39232953
Cochgro_totalfam_all_expanded_genesets.ermine.results:! heat shock protein binding GO:0031072 0.48195864
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! heat shock protein binding GO:0031072 0.41913938
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! response to heat GO:0009408 0.28787441
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! heat shock protein binding GO:0031072 0.23851475

grep temperature List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! response to temperature stimulus GO:0009266 0.61643874
Cochgro_totalfam_all_expanded_genesets.ermine.results:! response to temperature stimulus GO:0009266 0.8689807
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! response to temperature stimulus GO:0009266 0.51480601

grep epigenetic List_all_GO.txt
Cochgro_totalfam_all_expanded_genesets.ermine.results:! epigenetic regulation of gene expression GO:0040029 1.623e-03
Cochgro_totalfam_all_expanded_genesets.ermine.results:! negative regulation of gene expression, epigenetic GO:0045814 1.207e-03
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! epigenetic regulation of gene expression GO:0040029 5.129e-03
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! negative regulation of gene expression, epigenetic GO:0045814 1.043e-03
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! epigenetic regulation of gene expression GO:0040029 0.33555383
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! negative regulation of gene expression, epigenetic GO:0045814 0.26759899
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! epigenetic regulation of gene expression GO:0040029 1.166e-04
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! negative regulation of gene expression, epigenetic GO:0045814 2.926e-03

grep tropism List_all_GO.txt
Cochgro_totalfam_all_expanded_genesets.ermine.results:! gravitropism GO:0009630 0.6887685
Cochgro_totalfam_all_expanded_genesets.ermine.results:! negative gravitropism GO:0009959 0.49502009
Cochgro_totalfam_all_expanded_genesets.ermine.results:! tropism GO:0009606 0.94844852

grep glyco List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! alkylbase DNA N-glycosylase activity GO:0003905 0.03396649
Cochgro_totalfam_all_contracted_genesets.ermine.results:! DNA-3-methyladenine glycosylase activity GO:0008725 0.02611688
Cochgro_totalfam_all_contracted_genesets.ermine.results:! DNA-3-methylbase glycosylase activity GO:0043733 0.02611688
Cochgro_totalfam_all_contracted_genesets.ermine.results:! DNA N-glycosylase activity GO:0019104 0.14537867
Cochgro_totalfam_all_contracted_genesets.ermine.results:! glycolipid biosynthetic process GO:0009247 0.72843748
Cochgro_totalfam_all_contracted_genesets.ermine.results:! glycolipid metabolic process GO:0006664 0.75797268
Cochgro_totalfam_all_contracted_genesets.ermine.results:! glycolytic process GO:0006096 0.38307034
Cochgro_totalfam_all_contracted_genesets.ermine.results:! hydrolase activity, hydrolyzing N-glycosyl compounds GO:0016799 0.19954174
Cochgro_totalfam_all_contracted_genesets.ermine.results:! macromolecule glycosylation GO:0043413 0.01215827
Cochgro_totalfam_all_contracted_genesets.ermine.results:! protein glycosylation GO:0006486 0.01215827
Cochgro_totalfam_all_contracted_genesets.ermine.results:! protein N-linked glycosylation GO:0006487 0.41518221
Cochgro_totalfam_all_contracted_genesets.ermine.results:! protein O-linked glycosylation GO:0006493 0.01912664
Cochgro_totalfam_all_expanded_genesets.ermine.results:! catalytic activity, acting on a glycoprotein GO:0140103 0.27937733
Cochgro_totalfam_all_expanded_genesets.ermine.results:! glycogen (starch) synthase activity GO:0004373 0.18832431
Cochgro_totalfam_all_expanded_genesets.ermine.results:! glycolipid biosynthetic process GO:0009247 0.95427503
Cochgro_totalfam_all_expanded_genesets.ermine.results:! glycolipid metabolic process GO:0006664 0.96820768
Cochgro_totalfam_all_expanded_genesets.ermine.results:! glycolytic process GO:0006096 0.09236191
Cochgro_totalfam_all_expanded_genesets.ermine.results:! lipid glycosylation GO:0030259 0.03562739
Cochgro_totalfam_all_expanded_genesets.ermine.results:! macromolecule glycosylation GO:0043413 0.93013376
Cochgro_totalfam_all_expanded_genesets.ermine.results:! protein glycosylation GO:0006486 0.93013376
Cochgro_totalfam_all_expanded_genesets.ermine.results:! protein N-linked glycosylation GO:0006487 0.28989824
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! glycolytic process GO:0006096 0.59282772
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! macromolecule glycosylation GO:0043413 0.44532437
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! protein glycosylation GO:0006486 0.44532437
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! catalytic activity, acting on a glycoprotein GO:0140103 0.08511246
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! glycolipid biosynthetic process GO:0009247 0.37798032
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! glycolipid metabolic process GO:0006664 0.43635196
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! glycolytic process GO:0006096 0.78836572
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! glycosyl compound metabolic process GO:1901657 0.10553209
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! macromolecule glycosylation GO:0043413 0.70158386
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! protein glycosylation GO:0006486 0.70158386
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! macromolecule glycosylation GO:0043413 0.85215413
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! protein glycosylation GO:0006486 0.85215413
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! protein N-linked glycosylation GO:0006487 0.16066531
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! glycolipid biosynthetic process GO:0009247 0.31550507
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! glycolipid metabolic process GO:0006664 0.33006423
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! macromolecule glycosylation GO:0043413 0.03860673
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! protein glycosylation GO:0006486 0.03860673
Oxydig_totalfam_all_contracted_genesets.ermine.results:! catalytic activity, acting on a glycoprotein GO:0140103 0.0726913
Oxydig_totalfam_all_contracted_genesets.ermine.results:! glycogen (starch) synthase activity GO:0004373 0.03948822
Oxydig_totalfam_all_contracted_genesets.ermine.results:! glycolytic process GO:0006096 6.446e-04
Oxydig_totalfam_all_contracted_genesets.ermine.results:! protein N-linked glycosylation GO:0006487 0.0116136
Oxydig_totalfam_all_expanded_genesets.ermine.results:! glycogen (starch) synthase activity GO:0004373 0.57107678
Oxydig_totalfam_all_expanded_genesets.ermine.results:! glycolipid biosynthetic process GO:0009247 0.34187861
Oxydig_totalfam_all_expanded_genesets.ermine.results:! glycolipid metabolic process GO:0006664 0.40189286
Oxydig_totalfam_all_expanded_genesets.ermine.results:! glycolytic process GO:0006096 0.76794664

grep oxidoreductase List_all_GO.txt
125
Cochgro_totalfam_all_contracted_genesets.ermine.results:! hydroquinone:oxygen oxidoreductase activity GO:0052716 6.954e-03
Cochgro_totalfam_all_contracted_genesets.ermine.results:! intramolecular oxidoreductase activity GO:0016860 0.72843748
Cochgro_totalfam_all_contracted_genesets.ermine.results:! intramolecular oxidoreductase activity, transposing S-S bonds GO:0016864 0.29164842
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, disulfide as acceptor GO:0016671 9.016e-03
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors GO:0016667 0.30023511
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors GO:0016679 0.0326917
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors, oxygen as acceptor GO:0016682 0.01714934
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on NAD(P)H GO:0016651 0.61756172
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor GO:0050664 0.04260119
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen GO:0016709 1.486e-05
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen GO:0016701 0.78429988
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen GO:0016702 0.70677431
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on superoxide radicals as acceptor GO:0016721 0.20535263
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors GO:0016903 0.44378768
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor GO:0016620 0.30023511
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors GO:0016627 0.51142375
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, NAD or NADP as acceptor GO:0016628 0.69530404
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, oxygen as acceptor GO:0016634 0.26397651
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, with a flavin as acceptor GO:0052890 0.4371767
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH group of donors GO:0016645 0.0326917
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH group of donors, NAD or NADP as acceptor GO:0016646 0.29164842
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor GO:0016616 0.1742375
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-OH group of donors, oxygen as acceptor GO:0016899 6.954e-03
Cochgro_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase complex GO:1990204 0.64646417
Cochgro_totalfam_all_expanded_genesets.ermine.results:! disulfide oxidoreductase activity GO:0015036 0.40816356
Cochgro_totalfam_all_expanded_genesets.ermine.results:! intramolecular oxidoreductase activity GO:0016860 0.30535625
Cochgro_totalfam_all_expanded_genesets.ermine.results:! intramolecular oxidoreductase activity, transposing S-S bonds GO:0016864 3.471e-04
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors GO:0016667 0.83517235
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors GO:0016679 0.48195864
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on NAD(P)H GO:0016651 0.31961803
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on NAD(P)H, heme protein as acceptor GO:0016653 1.469e-03
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor GO:0016655 0.8689807
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen GO:0016709 0.83607859
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water GO:0016717 0.21712634
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen GO:0016701 0.03599685
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of one atom of oxygen (internal monooxygenases or internal mixed function oxidases) GO:0016703 0.13742534
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen GO:0016702 0.05754599
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on superoxide radicals as acceptor GO:0016721 2.36e-04
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors GO:0016903 0.37833506
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor GO:0016620 0.14110842
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors GO:0016627 0.09852562
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, NAD or NADP as acceptor GO:0016628 0.2327817
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, oxygen as acceptor GO:0016634 0.01466469
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH2 group of donors GO:0016638 0.36828606
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH2 group of donors, oxygen as acceptor GO:0016641 0.10986683
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor GO:0016616 0.80689442
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-OH group of donors, oxygen as acceptor GO:0016899 0.44802441
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase complex GO:1990204 0.98877538
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen GO:0016709 0.42729377
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, disulfide as acceptor GO:0016624 0.01178669
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors GO:0016903 0.0897386
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor GO:0016620 0.59282772
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors GO:0016627 0.64036492
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, NAD or NADP as acceptor GO:0016628 0.31028195
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! disulfide oxidoreductase activity GO:0015036 5.278e-07
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! glutathione disulfide oxidoreductase activity GO:0015038 8.267e-05
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! hydroquinone:oxygen oxidoreductase activity GO:0052716 0.28787441
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, disulfide as acceptor GO:0016671 0.31843903
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors GO:0016667 2.343e-09
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, NAD(P) as acceptor GO:0016668 6.958e-03
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, oxygen as acceptor GO:0016670 0.08293735
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, quinone or similar compound as acceptor GO:0016672 9.996e-06
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors GO:0016679 0.57937354
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors, oxygen as acceptor GO:0016682 0.46825665
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen GO:0016709 2.487e-05
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen GO:0016701 0.47662986
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen GO:0016702 0.26516508
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on superoxide radicals as acceptor GO:0016721 7.37e-04
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors GO:0016903 0.11362462
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor GO:0016620 0.05077038
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors GO:0016627 0.8554258
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, NAD or NADP as acceptor GO:0016628 0.23955745
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH2 group of donors GO:0016638 0.20143966
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH group of donors GO:0016645 0.20292671
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-OH group of donors, oxygen as acceptor GO:0016899 0.01067355
Drabaniv_totalfam_rapidly_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen GO:0016709 6.844e-12
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! hydroquinone:oxygen oxidoreductase activity GO:0052716 0.0718258
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors GO:0016679 0.15244685
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors, oxygen as acceptor GO:0016682 0.13909459
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on NAD(P)H GO:0016651 0.59970682
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor GO:0016655 0.3091611
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen GO:0016709 0.39721889
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen GO:0016701 0.01260806
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of one atom of oxygen (internal monooxygenases or internal mixed function oxidases) GO:0016703 0.0927099
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen GO:0016702 0.02732131
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors GO:0016903 0.0526045
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor GO:0016620 0.03401809
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH2 group of donors GO:0016638 0.28172295
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH2 group of donors, oxygen as acceptor GO:0016641 0.17684557
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase complex GO:1990204 0.53212756
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, disulfide as acceptor GO:0016671 2.785e-04
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors GO:0016667 1.369e-05
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, NAD(P) as acceptor GO:0016668 0.01279135
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, oxygen as acceptor GO:0016670 0.020732
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen GO:0016709 0.27147736
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water GO:0016717 5.202e-04
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors GO:0016627 0.8101796
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase complex GO:1990204 0.45554868
Oxydig_totalfam_all_contracted_genesets.ermine.results:! hydroquinone:oxygen oxidoreductase activity GO:0052716 0.06084413
Oxydig_totalfam_all_contracted_genesets.ermine.results:! intramolecular oxidoreductase activity GO:0016860 0.57724077
Oxydig_totalfam_all_contracted_genesets.ermine.results:! intramolecular oxidoreductase activity, interconverting aldoses and ketoses GO:0016861 0.36351032
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors GO:0016679 0.11193556
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors, oxygen as acceptor GO:0016682 0.0916894
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen GO:0016701 0.20115831
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen GO:0016702 0.10506282
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, disulfide as acceptor GO:0016624 0.03470258
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors GO:0016903 2.37e-03
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor GO:0016620 0.0215512
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors GO:0016627 0.74251006
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors, oxygen as acceptor GO:0016634 0.15806621
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH2 group of donors GO:0016638 0.26002518
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH2 group of donors, oxygen as acceptor GO:0016641 0.15806621
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase activity, acting on the CH-NH group of donors GO:0016645 0.26002518
Oxydig_totalfam_all_contracted_genesets.ermine.results:! oxidoreductase complex GO:1990204 0.24042212
Oxydig_totalfam_all_expanded_genesets.ermine.results:! disulfide oxidoreductase activity GO:0015036 0.28565391
Oxydig_totalfam_all_expanded_genesets.ermine.results:! hydroquinone:oxygen oxidoreductase activity GO:0052716 0.280092
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors, disulfide as acceptor GO:0016671 0.17719331
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on a sulfur group of donors GO:0016667 0.2017709
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors GO:0016679 0.44043589
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on diphenols and related substances as donors, oxygen as acceptor GO:0016682 0.38215124
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen GO:0016701 0.06048587
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen GO:0016702 0.01234741
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors GO:0016903 0.27365276
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor GO:0016620 0.13775197
Oxydig_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on the CH-CH group of donors GO:0016627 0.86718634

grep water List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! water-soluble vitamin biosynthetic process GO:0042364 0.5858575
Cochgro_totalfam_all_contracted_genesets.ermine.results:! water-soluble vitamin metabolic process GO:0006767 0.63087304
Cochgro_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water GO:0016717 0.21712634
Cochgro_totalfam_all_expanded_genesets.ermine.results:! response to water GO:0009415 0.03522927
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! water-soluble vitamin biosynthetic process GO:0042364 0.18809472
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! water-soluble vitamin metabolic process GO:0006767 1.329e-03
Drabaniv_totalfam_rapidly_expanded_genesets.ermine.results:! water-soluble vitamin biosynthetic process GO:0042364 0.02590635
Drabaniv_totalfam_rapidly_expanded_genesets.ermine.results:! water-soluble vitamin metabolic process GO:0006767 0.03735632
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! water-soluble vitamin biosynthetic process GO:0042364 0.43144272
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! water-soluble vitamin metabolic process GO:0006767 0.45317493
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water GO:0016717 5.202e-04
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! response to water deprivation GO:0009414 7.465e-05
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! response to water GO:0009415 9.02e-04
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! water-soluble vitamin biosynthetic process GO:0042364 1.11e-05
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! water-soluble vitamin metabolic process GO:0006767 1.904e-05
Oxydig_totalfam_all_expanded_genesets.ermine.results:! response to water GO:0009415 4.876e-03
Oxydig_totalfam_all_expanded_genesets.ermine.results:! water-soluble vitamin biosynthetic process GO:0042364 0.38215124
Oxydig_totalfam_all_expanded_genesets.ermine.results:! water-soluble vitamin metabolic process GO:0006767 0.42132973

grep "amino acid" List_all_GO.txt

Cochgro_totalfam_all_contracted_genesets.ermine.results:! alpha-amino acid biosynthetic process GO:1901607 0.07268272
Cochgro_totalfam_all_contracted_genesets.ermine.results:! amino acid activation GO:0043038 0.20382498
Cochgro_totalfam_all_contracted_genesets.ermine.results:! amino acid kinase activity GO:0019202 0.23522525
Cochgro_totalfam_all_contracted_genesets.ermine.results:! amino acid transmembrane transporter activity GO:0015171 7.391e-05
Cochgro_totalfam_all_contracted_genesets.ermine.results:! aromatic amino acid family biosynthetic process GO:0009073 0.25558143
Cochgro_totalfam_all_contracted_genesets.ermine.results:! aromatic amino acid metabolic process GO:0009072 0.218666
Cochgro_totalfam_all_contracted_genesets.ermine.results:! aspartate family amino acid biosynthetic process GO:0009067 0.32583803
Cochgro_totalfam_all_contracted_genesets.ermine.results:! aspartate family amino acid metabolic process GO:0009066 0.12586102
Cochgro_totalfam_all_contracted_genesets.ermine.results:! branched-chain amino acid biosynthetic process GO:0009082 0.01139765
Cochgro_totalfam_all_contracted_genesets.ermine.results:! branched-chain amino acid metabolic process GO:0009081 0.04246775
Cochgro_totalfam_all_contracted_genesets.ermine.results:! cellular modified amino acid biosynthetic process GO:0042398 0.47871863
Cochgro_totalfam_all_contracted_genesets.ermine.results:! cellular modified amino acid metabolic process GO:0006575 0.90735567
Cochgro_totalfam_all_contracted_genesets.ermine.results:! glutamine family amino acid biosynthetic process GO:0009084 0.56966575
Cochgro_totalfam_all_contracted_genesets.ermine.results:! glutamine family amino acid metabolic process GO:0009064 0.71781339
Cochgro_totalfam_all_contracted_genesets.ermine.results:! L-amino acid biosynthetic process GO:0170034 0.24996095
Cochgro_totalfam_all_contracted_genesets.ermine.results:! peptidyl-amino acid modification GO:0018193 0.75887518
Cochgro_totalfam_all_contracted_genesets.ermine.results:! proteinogenic amino acid biosynthetic process GO:0170038 0.24996095
Cochgro_totalfam_all_contracted_genesets.ermine.results:! sulfur amino acid biosynthetic process GO:0000097 0.68338578
Cochgro_totalfam_all_contracted_genesets.ermine.results:! sulfur amino acid metabolic process GO:0000096 0.72843748
Cochgro_totalfam_all_expanded_genesets.ermine.results:! alpha-amino acid biosynthetic process GO:1901607 0.66999821
Cochgro_totalfam_all_expanded_genesets.ermine.results:! alpha-amino acid catabolic process GO:1901606 0.61622694
Cochgro_totalfam_all_expanded_genesets.ermine.results:! amino acid activation GO:0043038 0.74754246
Cochgro_totalfam_all_expanded_genesets.ermine.results:! amino acid catabolic process GO:0009063 0.44068672
Cochgro_totalfam_all_expanded_genesets.ermine.results:! amino acid salvage GO:0043102 0.13742534
Cochgro_totalfam_all_expanded_genesets.ermine.results:! amino acid transmembrane transporter activity GO:0015171 0.25104232
Cochgro_totalfam_all_expanded_genesets.ermine.results:! aromatic amino acid family biosynthetic process GO:0009073 0.7001632
Cochgro_totalfam_all_expanded_genesets.ermine.results:! aromatic amino acid family catabolic process GO:0009074 0.14185436
Cochgro_totalfam_all_expanded_genesets.ermine.results:! aromatic amino acid metabolic process GO:0009072 0.53062142
Cochgro_totalfam_all_expanded_genesets.ermine.results:! aspartate family amino acid biosynthetic process GO:0009067 0.11445638
Cochgro_totalfam_all_expanded_genesets.ermine.results:! aspartate family amino acid metabolic process GO:0009066 0.14736326
Cochgro_totalfam_all_expanded_genesets.ermine.results:! branched-chain amino acid biosynthetic process GO:0009082 0.08639179
Cochgro_totalfam_all_expanded_genesets.ermine.results:! branched-chain amino acid metabolic process GO:0009081 0.3114937
Cochgro_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid biosynthetic process GO:0042398 0.40701337
Cochgro_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid catabolic process GO:0042219 1.606e-04
Cochgro_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid metabolic process GO:0006575 0.12298833
Cochgro_totalfam_all_expanded_genesets.ermine.results:! L-amino acid biosynthetic process GO:0170034 0.87219876
Cochgro_totalfam_all_expanded_genesets.ermine.results:! non-proteinogenic amino acid metabolic process GO:0170041 0.05373233
Cochgro_totalfam_all_expanded_genesets.ermine.results:! peptidyl-amino acid modification GO:0018193 0.09852562
Cochgro_totalfam_all_expanded_genesets.ermine.results:! proteinogenic amino acid biosynthetic process GO:0170038 0.87219876
Cochgro_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid biosynthetic process GO:0009070 0.8689807
Cochgro_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid metabolic process GO:0009069 0.9718662
Cochgro_totalfam_all_expanded_genesets.ermine.results:! sulfur amino acid biosynthetic process GO:0000097 0.58976253
Cochgro_totalfam_all_expanded_genesets.ermine.results:! sulfur amino acid metabolic process GO:0000096 0.69045842
Cochgro_totalfam_rapidly_expanded_genesets.ermine.results:! cellular modified amino acid catabolic process GO:0042219 1.222e-13
Cochgro_totalfam_rapidly_expanded_genesets.ermine.results:! cellular modified amino acid metabolic process GO:0006575 1.39e-07
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! alpha-amino acid catabolic process GO:1901606 0.31028195
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! amino acid catabolic process GO:0009063 0.36168334
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! aromatic amino acid family catabolic process GO:0009074 0.14334833
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! aromatic amino acid metabolic process GO:0009072 0.16783557
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! branched-chain amino acid biosynthetic process GO:0009082 0.21932535
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! branched-chain amino acid metabolic process GO:0009081 3.199e-04
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! cellular modified amino acid catabolic process GO:0042219 3.966e-04
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! cellular modified amino acid metabolic process GO:0006575 0.10200567
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! non-proteinogenic amino acid metabolic process GO:0170041 0.25475639
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! pyruvate family amino acid biosynthetic process GO:0009079 0.08864675
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! pyruvate family amino acid metabolic process GO:0009078 0.08864675
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! amino acid catabolic process GO:0009063 0.61878649
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! aromatic amino acid family biosynthetic process GO:0009073 0.35820449
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! aromatic amino acid metabolic process GO:0009072 0.67042987
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! aspartate family amino acid biosynthetic process GO:0009067 0.52841079
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! aspartate family amino acid metabolic process GO:0009066 0.56282194
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid catabolic process GO:0042219 0.15461266
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid metabolic process GO:0006575 0.25809372
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! glutamine family amino acid biosynthetic process GO:0009084 0.35820449
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! glutamine family amino acid metabolic process GO:0009064 0.65632968
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! non-proteinogenic amino acid metabolic process GO:0170041 0.14753595
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! peptidyl-amino acid modification GO:0018193 0.99761048
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! alpha-amino acid biosynthetic process GO:1901607 0.29799566
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! amino acid activation GO:0043038 0.46446082
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! amino acid kinase activity GO:0019202 9.601e-03
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! branched-chain amino acid metabolic process GO:0009081 0.05599521
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! glutamine family amino acid biosynthetic process GO:0009084 8.327e-03
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! glutamine family amino acid metabolic process GO:0009064 0.01960746
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! L-amino acid biosynthetic process GO:0170034 0.24063223
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! proteinogenic amino acid biosynthetic process GO:0170038 0.24063223
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! alpha-amino acid biosynthetic process GO:1901607 0.7506088
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! alpha-amino acid catabolic process GO:1901606 0.12833238
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! amino acid activation GO:0043038 2.867e-06
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! amino acid catabolic process GO:0009063 0.1556155
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! amino acid kinase activity GO:0019202 0.03673755
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! amino acid transmembrane transporter activity GO:0015171 2.574e-12
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! aspartate family amino acid metabolic process GO:0009066 0.27147736
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! branched-chain amino acid biosynthetic process GO:0009082 2.785e-04
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! branched-chain amino acid metabolic process GO:0009081 8.675e-04
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid catabolic process GO:0042219 0.020732
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid metabolic process GO:0006575 0.81584611
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! L-amino acid catabolic process GO:0170035 0.06679051
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! peptidyl-amino acid modification GO:0018193 0.88501021
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! proteinogenic amino acid catabolic process GO:0170040 0.06679051
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid catabolic process GO:0009071 0.01418755
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid metabolic process GO:0009069 0.44223419
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! sulfur amino acid metabolic process GO:0000096 0.33006423
Dryasoct_totalfam_rapidly_expanded_genesets.ermine.results:! amino acid activation GO:0043038 1.569e-03
Dryasoct_totalfam_rapidly_expanded_genesets.ermine.results:! amino acid transmembrane transporter activity GO:0015171 4.623e-22
Oxydig_totalfam_all_contracted_genesets.ermine.results:! amino acid activation GO:0043038 0.32724742
Oxydig_totalfam_all_contracted_genesets.ermine.results:! amino acid transmembrane transporter activity GO:0015171 0.13975964
Oxydig_totalfam_all_contracted_genesets.ermine.results:! amino acid transport GO:0006865 6.4e-03
Oxydig_totalfam_all_contracted_genesets.ermine.results:! aromatic amino acid family biosynthetic process GO:0009073 0.59507157
Oxydig_totalfam_all_contracted_genesets.ermine.results:! aromatic amino acid metabolic process GO:0009072 0.36598418
Oxydig_totalfam_all_contracted_genesets.ermine.results:! aspartate family amino acid biosynthetic process GO:0009067 0.47562813
Oxydig_totalfam_all_contracted_genesets.ermine.results:! aspartate family amino acid metabolic process GO:0009066 0.49773393
Oxydig_totalfam_all_contracted_genesets.ermine.results:! branched-chain amino acid metabolic process GO:0009081 0.06084413
Oxydig_totalfam_all_contracted_genesets.ermine.results:! L-amino acid biosynthetic process GO:0170034 0.87095032
Oxydig_totalfam_all_contracted_genesets.ermine.results:! peptidyl-amino acid modification GO:0018193 0.69402457
Oxydig_totalfam_all_contracted_genesets.ermine.results:! proteinogenic amino acid biosynthetic process GO:0170038 0.87095032
Oxydig_totalfam_all_contracted_genesets.ermine.results:! sulfur amino acid biosynthetic process GO:0000097 0.45255175
Oxydig_totalfam_all_contracted_genesets.ermine.results:! sulfur amino acid metabolic process GO:0000096 0.48679976
Oxydig_totalfam_all_expanded_genesets.ermine.results:! alpha-amino acid catabolic process GO:1901606 0.07259251
Oxydig_totalfam_all_expanded_genesets.ermine.results:! amino acid activation GO:0043038 0.35656025
Oxydig_totalfam_all_expanded_genesets.ermine.results:! amino acid catabolic process GO:0009063 0.1168531
Oxydig_totalfam_all_expanded_genesets.ermine.results:! aspartate family amino acid biosynthetic process GO:0009067 0.49556479
Oxydig_totalfam_all_expanded_genesets.ermine.results:! aspartate family amino acid metabolic process GO:0009066 0.5303422
Oxydig_totalfam_all_expanded_genesets.ermine.results:! cellular modified amino acid metabolic process GO:0006575 0.31863329
Oxydig_totalfam_all_expanded_genesets.ermine.results:! glutamine family amino acid metabolic process GO:0009064 0.24045264
Oxydig_totalfam_all_expanded_genesets.ermine.results:! L-amino acid biosynthetic process GO:0170034 0.4220437
Oxydig_totalfam_all_expanded_genesets.ermine.results:! L-amino acid catabolic process GO:0170035 0.04567173
Oxydig_totalfam_all_expanded_genesets.ermine.results:! proteinogenic amino acid biosynthetic process GO:0170038 0.4220437
Oxydig_totalfam_all_expanded_genesets.ermine.results:! proteinogenic amino acid catabolic process GO:0170040 0.04567173
Oxydig_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid biosynthetic process GO:0009070 0.02986301
Oxydig_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid metabolic process GO:0009069 0.07951264


grep "cellular modified amino acid catabolic process" List_all_GO.txt


grep "amino sugar catabolic process" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! amino sugar catabolic process GO:0046348 3.97e-03
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! amino sugar catabolic process GO:0046348 4.383e-06
Drabaniv_totalfam_rapidly_contracted_genesets.ermine.results:! amino sugar catabolic process GO:0046348 2.248e-10
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! amino sugar catabolic process GO:0046348 0.03294838
Oxydig_totalfam_all_expanded_genesets.ermine.results:! amino sugar catabolic process GO:0046348 4.791e-08


grep "serine family amino acid metabolic process" List_all_GO.txt
Cochgro_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid metabolic process GO:0009069 0.9718662
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid metabolic process GO:0009069 0.44223419
Oxydig_totalfam_all_expanded_genesets.ermine.results:! serine family amino acid metabolic process GO:0009069 0.07951264

#----------------------------
grep "nutrient" List_all_GO.txt

Cochgro_totalfam_all_contracted_genesets.ermine.results:! nutrient reservoir activity GO:0045735 6.319e-05
Cochgro_totalfam_all_contracted_genesets.ermine.results:! response to nutrient levels GO:0031667 0.4371767
Cochgro_totalfam_all_expanded_genesets.ermine.results:! nutrient reservoir activity GO:0045735 0.05619654
Cochgro_totalfam_all_expanded_genesets.ermine.results:! response to nutrient levels GO:0031667 0.32912348
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! nutrient reservoir activity GO:0045735 6.92e-05
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! nutrient reservoir activity GO:0045735 0.12734707
Oxydig_totalfam_all_contracted_genesets.ermine.results:! cellular response to nutrient levels GO:0031669 0.12105595
Oxydig_totalfam_all_contracted_genesets.ermine.results:! nutrient reservoir activity GO:0045735 8.834e-03
Oxydig_totalfam_all_contracted_genesets.ermine.results:! response to nutrient levels GO:0031667 0.26002518


#---------------------------------
# terms of interesting

# GO:0019953 sexual reproduction
grep "GO:0019953" List_all_GO.txt
# Cochgro_totalfam_all_contracted_genesets.ermine.results:! sexual reproduction GO:0019953 0.23522525
# Cochgro_totalfam_all_expanded_genesets.ermine.results:! sexual reproduction GO:0019953 0.05619654
# Cochgro_totalfam_rapidly_contracted_genesets.ermine.results:! sexual reproduction GO:0019953 0.01825699
# Oxydig_totalfam_all_expanded_genesets.ermine.results:! sexual reproduction GO:0019953 1.235e-04
# Oxydig_totalfam_rapidly_expanded_genesets.ermine.results:! sexual reproduction GO:0019953 3.434e-08

# response to wounding
grep "GO:0009611" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! response to wounding GO:0009611 0.41518221
Cochgro_totalfam_rapidly_contracted_genesets.ermine.results:! response to wounding GO:0009611 0.03618763
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! response to wounding GO:0009611 5.081e-06

# oxidative stress
grep "GO:0034599" List_all_GO.txt
Cochgro_totalfam_all_expanded_genesets.ermine.results:! cellular response to oxidative stress GO:0034599 0.13299729
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! cellular response to oxidative stress GO:0034599 0.21932535
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! cellular response to oxidative stress GO:0034599 1.263e-03

# response to hormone
# GO:0009725
grep "GO:0009725" List_all_GO.txt
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! response to hormone GO:0009725 0.77272625
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! response to hormone GO:0009725 4.477e-09

# reproductive process
grep "GO:0022414" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! reproductive process GO:0022414 0.22770592
Cochgro_totalfam_all_expanded_genesets.ermine.results:! reproductive process GO:0022414 0.02401169
Cochgro_totalfam_rapidly_contracted_genesets.ermine.results:! reproductive process GO:0022414 0.15749389
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! reproductive process GO:0022414 0.58646132
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! reproductive process GO:0022414 0.01802574
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! reproductive process GO:0022414 0.07044312
Oxydig_totalfam_all_expanded_genesets.ermine.results:! reproductive process GO:0022414 4.435e-05
Oxydig_totalfam_rapidly_expanded_genesets.ermine.results:! reproductive process GO:0022414 1.739e-04

# 	methylation
grep "GO:0032259" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! methylation GO:0032259 0.52994434
Cochgro_totalfam_all_expanded_genesets.ermine.results:! methylation GO:0032259 0.19112328
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! methylation GO:0032259 0.01465912
Drabaniv_totalfam_rapidly_expanded_genesets.ermine.results:! methylation GO:0032259 2.493e-05
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! methylation GO:0032259 0.04609878
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! methylation GO:0032259 0.76599638
Oxydig_totalfam_all_contracted_genesets.ermine.results:! methylation GO:0032259 0.79247974
Oxydig_totalfam_all_expanded_genesets.ermine.results:! methylation GO:0032259 0.37191038

# biological process involved in interspecies interaction between organisms
grep "GO:0044419" List_all_GO.txt

Cochgro_totalfam_all_contracted_genesets.ermine.results:! biological process involved in interspecies interaction between organisms GO:0044419 0.4136777
Cochgro_totalfam_all_expanded_genesets.ermine.results:! biological process involved in interspecies interaction between organisms GO:0044419 0.99987557
Cochgro_totalfam_rapidly_contracted_genesets.ermine.results:! biological process involved in interspecies interaction between organisms GO:0044419 0.01156649
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! biological process involved in interspecies interaction between organisms GO:0044419 0.94555683
Oxydig_totalfam_all_contracted_genesets.ermine.results:! biological process involved in interspecies interaction between organisms GO:0044419 0.6873559
Oxydig_totalfam_all_expanded_genesets.ermine.results:! biological process involved in interspecies interaction between organisms GO:0044419 1.396e-04

# anther wall tapetum development
grep "GO:0048658" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! anther wall tapetum development GO:0048658 4.992e-04
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! anther wall tapetum development GO:0048658 9.996e-06
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! anther wall tapetum development GO:0048658 1.136e-05

# recognition of pollen
grep "GO:0048544" List_all_GO.txt
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! recognition of pollen GO:0048544 5.824e-09
Oxydig_totalfam_all_expanded_genesets.ermine.results:! recognition of pollen GO:0048544 0.06296949

# regulation of shoot system development
grep "GO:0048831" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! regulation of shoot system development GO:0048831 0.05353342
Cochgro_totalfam_all_expanded_genesets.ermine.results:! regulation of shoot system development GO:0048831 8.043e-03
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! regulation of shoot system development GO:0048831 0.02743667
Oxydig_totalfam_all_expanded_genesets.ermine.results:! regulation of shoot system development GO:0048831 0.08282567

grep "GO:0009909" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! regulation of flower development GO:0009909 0.05353342
Cochgro_totalfam_all_expanded_genesets.ermine.results:! regulation of flower development GO:0009909 8.043e-03
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! regulation of flower development GO:0009909 0.02743667
Oxydig_totalfam_all_expanded_genesets.ermine.results:! regulation of flower development GO:0009909 0.08282567

# regulation of reproductive process
grep "GO:2000241" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! regulation of reproductive process GO:2000241 0.05953276
Cochgro_totalfam_all_expanded_genesets.ermine.results:! regulation of reproductive process GO:2000241 0.01066476
Drabaniv_totalfam_all_contracted_genesets.ermine.results:! regulation of reproductive process GO:2000241 0.03055843
Oxydig_totalfam_all_expanded_genesets.ermine.results:! regulation of reproductive process GO:2000241 0.08282567

# response to water deprivation
grep "GO:0009414" List_all_GO.txt
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! response to water deprivation GO:0009414 7.465e-05

grep "GO:0009269" List_all_GO.txt
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! response to desiccation GO:0009269 7.465e-05

grep "GO:0009415" List_all_GO.txt
Cochgro_totalfam_all_expanded_genesets.ermine.results:! response to water GO:0009415 0.03522927
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! response to water GO:0009415 9.02e-04
Oxydig_totalfam_all_expanded_genesets.ermine.results:! response to water GO:0009415 4.876e-03

grep "defense" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! defense response to fungus GO:0050832 0.18575091
Cochgro_totalfam_all_contracted_genesets.ermine.results:! defense response to other organism GO:0098542 0.4136777
Cochgro_totalfam_all_expanded_genesets.ermine.results:! defense response to other organism GO:0098542 0.99987557
Cochgro_totalfam_all_expanded_genesets.ermine.results:! regulation of defense response GO:0031347 0.01426031
Cochgro_totalfam_all_expanded_genesets.ermine.results:! regulation of defense response to fungus GO:1900150 7.685e-04
Cochgro_totalfam_rapidly_contracted_genesets.ermine.results:! defense response to fungus GO:0050832 1.375e-03
Cochgro_totalfam_rapidly_contracted_genesets.ermine.results:! defense response to other organism GO:0098542 0.01156649
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! defense response to other organism GO:0098542 0.94555683
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! regulation of defense response GO:0031347 0.03806259
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! regulation of defense response to fungus GO:1900150 0.01515314
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! defense response to symbiont GO:0140546 0.45317493
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! regulation of defense response GO:0031347 0.14416788
Dryasoct_totalfam_all_expanded_genesets.ermine.results:! defense response to symbiont GO:0140546 0.03184547
Oxydig_totalfam_all_contracted_genesets.ermine.results:! defense response to other organism GO:0098542 0.6873559
Oxydig_totalfam_all_contracted_genesets.ermine.results:! defense response to symbiont GO:0140546 0.1935215
Oxydig_totalfam_all_contracted_genesets.ermine.results:! regulation of defense response GO:0031347 0.02187143
Oxydig_totalfam_all_contracted_genesets.ermine.results:! regulation of defense response to fungus GO:1900150 8.834e-03
Oxydig_totalfam_all_expanded_genesets.ermine.results:! defense response to fungus GO:0050832 8.734e-03
Oxydig_totalfam_all_expanded_genesets.ermine.results:! defense response to other organism GO:0098542 1.396e-04
Oxydig_totalfam_all_expanded_genesets.ermine.results:! defense response to symbiont GO:0140546 0.10097045

grep "DNA repair" List_all_GO.txt
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! regulation of DNA repair GO:0006282 6.026e-03

grep "recombination" List_all_GO.txt
Drabaniv_totalfam_all_expanded_genesets.ermine.results:! double-strand break repair via homologous recombination GO:0000724 0.04814268
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! DNA recombination GO:0006310 0.6706545
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! double-strand break repair via homologous recombination GO:0000724 0.14416788
Dryasoct_totalfam_all_contracted_genesets.ermine.results:! recombinational repair GO:0000725 0.16066531
Oxydig_totalfam_all_expanded_genesets.ermine.results:! DNA recombination GO:0006310 0.09331759
Oxydig_totalfam_all_expanded_genesets.ermine.results:! double-strand break repair via homologous recombination GO:0000724 0.10097045
Oxydig_totalfam_all_expanded_genesets.ermine.results:! recombinational repair GO:0000725 0.11907696

grep "telomere" List_all_GO.txt
Cochgro_totalfam_all_contracted_genesets.ermine.results:! telomere maintenance GO:0000723 0.70677431
Cochgro_totalfam_all_contracted_genesets.ermine.results:! telomere organization GO:0032200 0.70677431
Cochgro_totalfam_all_expanded_genesets.ermine.results:! telomere maintenance GO:0000723 0.94191746
Cochgro_totalfam_all_expanded_genesets.ermine.results:! telomere organization GO:0032200 0.94191746
Oxydig_totalfam_all_expanded_genesets.ermine.results:! telomere maintenance GO:0000723 0.47757064
Oxydig_totalfam_all_expanded_genesets.ermine.results:! telomere organization GO:0032200 0.47757064

########################################################

# format files
grep '^!' Oxydig_all_expanded_genesets.ermine.results > Oxydig_all_expanded_genesets_forR.ermine.results
# remove header

# Plotting
# https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

install.packages("GOplot")

library(GOplot)

# Load the dataset
path="/home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/"
ermineJ_results <- base::as.data.frame(utils::read.table(paste0(path,"ermine.results/Oxydig_all_expanded_genesets_forR.ermine.results"),sep="\t", header = FALSE, check.names = FALSE))
GO_mappings <- base::as.data.frame(utils::read.table(paste0(path,"Oxydig_GO_mappings.ermineJ.txt"),  sep="\t", header = TRUE, check.names = FALSE))


# Get a glimpse of the data format of the results of the functional analysis... 
head(EC$david) # ermine J currently is not outputting the list of genes? 
head(EC$genelist) #DE genes and DMRs could be put here

# generate matrix of relations
chord_dat(data, genes, process)

# Generate the plotting object
circ <- circle_dat(EC$david, EC$genelist)

# Generate a simple barplot
GOBar(subset(circ, category == 'BP'))


# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')


# Facet the barplot, add a title and change the colour scale for the z-score
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))


# Generate the bubble plot with a label threshold of 3
GOBubble(circ, labels = 3)


# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)

# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)

# Venn diagram comparing groups
l1 <- subset(circ, term == 'heart development', c(genes,logFC))
l2 <- subset(circ, term == 'plasma membrane', c(genes,logFC))
l3 <- subset(circ, term == 'tissue morphogenesis', c(genes,logFC))
GOVenn(l1,l2,l3, label = c('heart development', 'plasma membrane', 'tissue morphogenesis'))


#################
# Future
# Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)

# Generate a circular visualization of selected terms
IDs <- c('GO:0007507', 'GO:0001568', 'GO:0001944', 'GO:0048729', 'GO:0048514', 'GO:0005886', 'GO:0008092', 'GO:0008047')
GOCircle(circ, nsub = IDs)

# Generate a circular visualization for 10 terms
GOCircle(circ, nsub = 10)

# Define a list of genes which you think are interesting to look at. The item EC$genes of the toy 
# sample contains the data frame of selected genes and their logFC. Have a look...
head(EC$genes)

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-8], nlfc = 0)

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

# Clustering
GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)
GOCluster(circ, EC$process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))








