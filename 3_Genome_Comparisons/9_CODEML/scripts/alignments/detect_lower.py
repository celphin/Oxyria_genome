# -*- coding: utf-8 -*-
"""

@author: chocobo

"""

# import modules

import os
import glob
import numpy as np

'''
*******************************************************************************
'''

def detect_lower06(path, thresh):
    
    #change to directory
    os.chdir(path)

    # get filenames
    files = glob.glob('*PRANK.Guidance2_res_pair_*')
    
    # create empty list
    lista = []

    # loop through files
    for f in files:
        
        # load data of one file
        data = np.genfromtxt(f,skip_header=1)
        second_column = data[:,1]
        
        # get values where it is under the threshold
        threshold_indices = np.where(second_column < thresh)[0]
        print (f)
        print (threshold_indices)
        print (threshold_indices.shape[0])
        # if not equal to zero then add to the list
        if threshold_indices.shape[0] == 0:
            None
        else:
            lista.append(f.split(".")[0])
            
    
    # save file to badseqscores.txt
    np.savetxt('badseqscores.txt',lista,fmt='%s')
    
    print ("DONE")
    
    
    
'''
*******************************************************************************
'''

# specify path and value of the threshold
path = '/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/GUIDANCE/Seq_Scores_Guidance/'
thresh = 0.6

# run command
detect_lower06(path, thresh)
