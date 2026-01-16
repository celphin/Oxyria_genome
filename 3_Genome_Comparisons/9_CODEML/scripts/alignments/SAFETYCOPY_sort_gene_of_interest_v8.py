"""

@author: chocobo

"""

# import modules

import os
import glob

'''

*******************************************************************************

'''



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

    

    # Here, I delete the existing fasta files (if this code was run before)
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



#path that contains the .aligned.fasta

path = "/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/orthofasta/"

#change directory to path
os.chdir(path)

#get all .aligned.fasta files
f = glob.glob("*[0-9].aligned.fa")



# give pattern name

pattern = r'Poptre'

# loop through all fasta files

for fname in f:

    get_gene_of_interest(path,fname,pattern)

    

'''

*******************************************************************************

'''





