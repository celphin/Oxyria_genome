#########################
# Running CODEML
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/CODEML_WORKFLOW
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# Steps
# Make CODEML control files
# Run CODEML
# Explore results

######################

# Load modules
module load cufflinks/2.2.1 
module load cdbfasta/2017-03-16
module load orthofinder/2.5.2 
module load diamond/2.0.4 
module load clustalo/1.2.4
module load pal2nal.v14 
module load paml/4.9h

##########################
# Understanding Codeml parameters
    # runmode = 0 means I am specifying user trees
    # seqtype = 1 means I am using codon sequences
    # CodonFreq = 0 means I am assuming equal codon frequencies (codon substitution model)
    # model = 0 whether ω should be allowed to vary among branches in the tree
    # NSsites = 0 1 2 7 8 whether ω should be allowed to vary among sites , several of these can be specified at once.
    # icode = 0 is specifies that I want to use the universal genetic code
    # fix_kappa = 0 – initial transition/transversion ratio value set to 0
    # kappa = 0 – assumption starting value for expectation of transitions/transversions
    # fix_omega = 1 – fixes omega to test for significance between current and foreground branch
    # omega = 0.001 – not 100% sure here, but likely the initial value set for omega in current branch
    # RateAncestor = 0 – ancestral sequence not constructed
    # clock = 0 – assumption of a molecular clock is not allowed, rates are free to vary from branch to branch
    # Small_Diff = 5e-08 – supposed to use multiple values to make sure your results are not sensitive to this number uses 1e-8 to 1e-5.
    # Malpha = 0 – when substitution rates are assumed to vary from site to site, Malpha specifies whether one gamma distribution will be applied across all sites
    # method = 0 – nucleotide substitution model
    # ndata = 1 – number of separate datasets


###############################
# Prepare Codeml control files and run Codeml

cd 20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

# Generate the control files appropriately named, but on one line
for f in OG*pal2nal; do echo "seqfile = ../"$f" \t treefile = ../"${f%_*}"\t outfile = "${f%.*}".out  \t verbose = 1 \t runmode = 0 \t seqtype = 1 \t CodonFreq = 0 \t model = 0 \t NSsites = 0 1 2 7 8 \t icode = 0 \t RateAncestor = 0 \t clock = 0 \t  \t method = 0 \t ndata = 1" > ${f%.*}"codeml.ctl";done

# Move orthgroup named files to new folders, 
# rename the control files to codeml.ctl, as codeml only works on codeml.ctl

for f in OG*ctl; do mkdir ${f}dir; mv $f ${f}dir/codeml.ctl; done# Transform the single line of tab separated values to a newline separated document

# Transform the single line of tab separated values to a newline separated document
sed -i 's/\\t/\n/g' *dir/codeml.ctl

# Run codeml – Creates an sh script that can be separated to submit to multiple nodes.
for f in *ctldir; do echo "ml paml/4.9h-m5kkb2e; cd "$f"; codeml; cd ../ ";done >IterateCodeml.sh


###############################
# Checking results

# Which orthogroups are demonstrating positive selection
    # ntime is the number of branch lengths
    # np is the number of parameters
    # The third value is your log likelihood value

cd 20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

# This prints the file name and grabs the values we want (np and Log Liklihood)
grep --with-filename "lnL"  OG*_treecodeml.ctldir/*out |awk '{print $1,$4,$5}' |sed 's/)://g' |tr " " "\t" |awk '{ if (NR % 5 == 0) {print $0"#"} else {print}}'  |tr "\n" "\t" |sed 's/#/\n/g' |awk '{print $1,$2,$3,$5,$6,$8,$9,$11,$12,$14,$15}' |less

#------------------
# Create separate files for each test, 
# as the number of parameters varys among and between models.

# See the rest of the tutorial ...









