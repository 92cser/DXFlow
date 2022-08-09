#!/usr/bin/python

# Removing unnecessary results when the writing interval is too small
# Results of the instants before 'reduceUntil' and after 'reduceAfter' are removed according to 'deltaT' 

import os
import shutil as shl

# time interval for store
deltaT = 0.5
# results before this instant will be deleted
reduceUntil = 67
# results after this instant will be deleted
reduceAfter = 74

# input the number of processors
nOfProc = input("The number of processors: ")

# current path
path = os.path.abspath(".")

# enter into each processor
for i in range(0, int(nOfProc)):
    currPath = path + "/processor" + str(i)
    print("Handling processor"+str(i)+"...")
    # all the documents in the path
    allDocu = os.listdir(currPath)
    # loop 
    for j in range(0, len(allDocu)):
    	# avoid constant file 
        if(allDocu[j] != "constant"):
            currInstant = float(allDocu[j])
            remainder = currInstant % deltaT
            # remove 
            if(remainder != 0 and (currInstant < reduceUntil or currInstant > reduceAfter)):
                shl.rmtree(currPath + "/" + allDocu[j])
