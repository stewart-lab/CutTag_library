# parse bed file to get 300bp + TSS

import sys

input_bed = open(sys.argv[1], "r")
output= open(str(sys.argv[1])+"_300bpTSS.bed", "w")
for line in input_bed:
    L=line.strip().split("\t")
    if L[5]=="+":
        newTSS = int(L[1])+300
        # check gene
        if newTSS > int(L[2]):
            continue
        else:
            newline = "\t".join([str(L[0]),str(newTSS)]+L[2:])
            #print(newline)
            output.write(newline+"\n")
    elif L[5]=="-":
        newTSS =int(L[2])-300
        if newTSS < int(L[1]):
            continue
        else:
            newline = "\t".join([str(L[0]),str(L[1]),str(newTSS)]+L[3:])
            #print(newline)
            output.write(newline+"\n")
    else:
        print(L[5],"direction can't be determined")
        
output.close()
input_bed.close()
        

