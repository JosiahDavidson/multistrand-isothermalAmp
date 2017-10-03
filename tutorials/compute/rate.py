from anneal import computeDissociationAndWriteToCL as computeDissociation
from anneal import computeAndWriteToCL as computeAnneal 

import sys, time

if(len(sys.argv) < 2):
    print("Please provide a DNA sequence as commandline argument")
    print("Specify the experiment in the first argument: <dissociation> or <association>  ")
    print("Add -bootstrap to do a boostrap ")
    print("Example: rate.py dissociation ATGCAGT -bootstrap")
    exit()

start_time = time.time()

type = sys.argv[1]
mySequence = sys.argv[2]

doBootstrap = False
if(len(sys.argv) > 3):
    if(str(sys.argv[3])=="-bootstrap"):
        doBootstrap = True

if (type ==  "dissociation"):

    result = computeDissociation(mySequence, doBootstrap )

if type == "association":
    
    result = computeAnneal(mySequence, doBootstrap)


print ("Computing took %.4f s" % (time.time() - start_time))