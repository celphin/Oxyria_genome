"""

@author: chocobo

"""

# import modules

import os
import glob

'''
*******************************************************************************
This script starts the process of dividing the orthogroup fastas into 
orthogroup subsets - ultimately resulting in alignments with just one gene copy 
per species. The process is slightly complicated (and just an approximation), 
but it starts with defining a species of interest in our dataset (here 
"Poptre" for Populus tremula) and then dividing each orthogroup into subsets 
based on the gene copies of that species. For instance, the script will make 
three orthogroup subsets called OG0000334_Poptre_001.fasta 
OG0000334_Poptre_002.fasta and OG0000334_Poptre_003.fasta if there are three 
gene copies from "Poptre" in orthogroup OG0000334. To find which gene copies 
from the other species that have the smallest genetic distance to each of these 
gene copies, we place each gene copy together with the rest of the orthogroup 
sequences and calculate protein distance. (Note that the species of interest 
is always put last in the file to make file processing easier downstream).

*******************************************************************************
'''

### DEFINE SOME VARIABLES ###

# A: Enter path that contains the .aligned.fa

path = "/cluster/home/siribi/work/CODEML_WORKFLOW_TEST/orthofasta/"

# B: Define the species abbreviation to search for. Protein distances will be calculated in relation to this species.

pattern = r'Prunavi'


def get_gene_of_interest(path,fname,pattern):
 

    '''
    example:
    path = "C:/Users/yellow_chocobo/Desktop/siri_python/version7/"
    fname = "OG0000113.aligned.fa"
    pattern = r'DnA'
    get_gene_of_interest(path,fname,pattern)     

    '''
    

    #change path
    os.chdir(path)
    

    #create empty lists to store line1 and line2
    line1 = []

    

    # NB: Here, I delete the existing fasta files (if this code was run before)
    # I do that because of the append function in open(paths + fasta_file,'a').
    # Every time you will run this code, it will append the same information, again
    # and again. By making sure that no fasta files are created from before, we avoid
    # repetitions. Be careful to not put other important files which ends with .fasta
    # in the paths directory. Otherwise it will be deleted by the three next lines
    
    #fasta_f = glob.glob("*" + pattern + "*.sorted.aligned.fa")

    #for file in fasta_f:
    #    os.remove(file)


    # number of lines in the file        
    num_lines = sum(1 for line in open(path + fname,'r'))
        

    # we first get the gene of interest (in the example all genes containing CbS)
    flag = 0
    nline = ""
    count = 0

    
    with open(path + fname,'r') as search:

        for line in search:
            count = count + 1

            if ((flag == 1) and (count == num_lines)):

                nline += line
                line1.append(nline)

            # if line starts with >
            elif line.startswith('>'):

                #if the flag was previously to 1
                if flag == 1:

                    #save everything that was into nline in line1
                    line1.append(nline)

                # new entry, flag is set to zero
                flag = 0

                # if it ends with the gene of interest, we set the flag to 1                 
                if line.startswith(">" + pattern):

                    flag = 1
                    nline = line
                    
            else:

                if flag == 1:
                    nline += line

    
    # get the name of all the genes e.g.: ['m.23033_CbS', 'm.30750_CbS', 'm.31135_CbS', 'm.31145_CbS']

    patterns = []
    
    for i in range(len(line1)):
        patterns.append(line1[i].split('\n')[0].split('>')[1])

            
    # loop through the aligned.fa    

    '''

    if the pattern is detected then nothing is written in the new files
    because it will be written in the second part of the code (next green
    explanation)

    '''

        

    for idx, p in enumerate(patterns):

        #npattern = re.compile(p)

        flag = 0
        fname2 = fname.split('.')[0] + '_' + pattern + '_' + str(idx).zfill(3) + '.sorted.aligned.fa'

               

        with open(path + fname,'r') as search, open(path + fname2 ,'a') as myfile:

            for line in search:

                if line.startswith('>' + pattern):
                    flag = 1

                elif flag == 1:
                    if line.startswith('>'):

                        if line.startswith('>' + pattern):
                            flag = 1

                        else:
                            flag = 0
                            myfile.writelines(line)

                else:
                    myfile.writelines(line)
                  

    '''

    In this part we set the gene of interest at the end

    '''                

    for idx, p in enumerate(patterns):

        

        fname2 = fname.split('.')[0] + '_' + pattern + '_' + str(idx).zfill(3) + '.sorted.aligned.fa'

        with open(path + fname2 ,'a') as myfile:

            myfile.writelines(line1[idx])

    
    print ("DONE")

    
'''

*******************************************************************************

'''


#change directory to path
os.chdir(path)

#get all .aligned.fa files
f = glob.glob("*[0-9].aligned.fa")

# loop through all fasta files

for fname in f:

    get_gene_of_interest(path,fname,pattern)

    

'''

*******************************************************************************

'''





