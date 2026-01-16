"""
@author: chocobo
"""

# import modules
import numpy as np, os
import re
import time
import glob

def chocobo(start,end):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

start_time = time.time()

### A: Enter folder containing the concatenated peptide file
path = '/cluster/home/siribi/work/EVOTREE/CODEML/alignments/All40/super_pep/' #enter directory with backslash at the end

### B: Enter folder where the new files will be saved to (created by this script)
paths ='/cluster/home/siribi/work/EVOTREE/CODEML/alignments/All40/orthofasta/' #enter directory with backslash at the end 

# change to directory containing the .fasta file
os.chdir(path)

### C: Enter filename of concatenated fasta file (peptide)
fdna = "All40_concat.pep"

# create plots directory if not existing
if not os.path.exists(paths):
    os.makedirs(paths)

os.chdir(paths)

# Here, I delete the existing fasta files (if this code was run before)
# I do that because of the append function in open(paths + fasta_file,'a').
# Every time you will run this code, it will append the same information, again
# and again. By making sure that no fasta files are created from before, we avoid
# repetitions. Be careful to not put other important files which ends with .fasta
# in the paths directory. Otherwise it will be deleted by the three next lines
fasta_f = glob.glob("*.fasta")

for file in fasta_f:
    os.remove(file)

### D: Enter directory where you can find all the Orthogroup index files
pathi = "/cluster/home/siribi/work/EVOTREE/CODEML/alignments/All40/index_files/"

# change directory to pathi
os.chdir(pathi)

# get all the files that starts with "OG" and finish with ".txt"
files = sorted(glob.glob('OG*.txt'))

# loop through all OG* files
for file in files:
    print file
    # load the data from OG* .txt
    index = np.genfromtxt(pathi + file ,delimiter=" ",dtype=None)

    # name of the fasta file to be created
    fasta_file = file.split(".")[0] + ".fasta"
        
    # each file contains several entries, so we also have to loop through the indexes
    for ix in index:
        flag = 0
        #register the pattern
        pattern = re.compile(ix)

        with open(path + fdna,'r') as search, open(paths + fasta_file,'a') as myfile:
            for line in search:
                if flag == 2:
                    break
                elif flag != 2:
                    if flag == 1:
                        myfile.writelines(line)
                        flag = 2
                    elif (pattern.search(line)) and (flag == 0):
                        new_line = line
                        myfile.writelines(new_line)
                        flag = 1


chocobo(start_time,time.time())



