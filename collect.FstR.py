#!/usr/bin/python

import re
import os

#DIR = "/home/users/xiaodong/Documents/Project/Seal/trouble_shooting_2/sfs/F0999P001/data/fst"
DIR = "/home/users/xiaodong/Documents/Project/Seal/Fst/data/fst"

POPS = [re.sub(".fst","",f) for f in os.listdir(DIR) if f.endswith("fst")]
dictA = {}
dictSort = {}
for pop in POPS:
    infile = os.path.join(DIR,pop+".fst")
    with open( infile, 'r' ) as f:
        s = f.read().rstrip('\n')
        weight = s
        dictA[pop] = [weight]
    f.close()


for key in dictA.keys():
    pop1,pop2 = key.split("_")
    key2 = "_".join (sorted([pop1,pop2]))
    dictSort[key2] = dictA[key]
    

    
pairlist= [f.split("_") for f in POPS]
pops = sorted ( {x for l in pairlist for x in l} )
#print(pops)
header = ['rowName'] + pops
print( "\t".join(header) )

for i in pops:
    print (i,end="\t")
    row = []
    for j in pops:
        if i==j:
            row.append("NA")
        elif i<j:
            skey = i+"_"+j
            if skey in dictSort:
                row.append( dictSort[skey][0] )
            else:
                row.append('Not found')
        else:
            skey = j+"_"+i
            if skey in dictSort:
                row.append( dictSort[skey][0] )
            else:
                row.append('Not found')
    print ("\t".join(row))



            
