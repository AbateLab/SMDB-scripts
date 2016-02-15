# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 12:18:45 2014

pileup reader - for reading pileups

@author: Freeman
"""

data = {"*":{}, "A":{}, "B":{}, "C":{}, "D":{}, "E":{}, "F":{}, "G":{}, "H":{}, "I":{}, "J":{}, "L":{} }
#template is key, dictionary of {position:[bases]}

def input_pileup (path, limit = 0):
    infile = open(path, "r")
    for line in infile:
        lines = line.split()    #lines[3] is the number for coverage
        data[lines[0]][lines[1]] = list(lines[4])        #list of covering bases append to the correct position
        limit -= 1
        if limit == 1:
            break
    return

def input_ref(infile):
    reference = {}
    for line in infile:
        if line[0] == ">":
            name = line[1:].strip()
            seq = infile.next().strip()
            reference[name] = seq
    return reference
    
def output_coverages_ends (outpath):
    for key in data:
        with open("%s%scov.txt"%(outpath,key), "w") as outfile:                   
            for p in data[key]:
                outfile.write("%s\t%d\n"%(p,len(data[key][p])))
    return

def percent_cov(): #prints out the of bases that the assemblies covered of the reference
    lengths = {"A":0, "B":0, "C":5227, "D":5001, "E":5002, "F":0, "G":3185, "H":4774, "I":4876, "J":3920, "L":4361}
    cov = 0
    total = 0
    for key in data:
        if key == "*":
            continue
        for pos in data[key]:
            cov += 1   
    for key in lengths:
        total += lengths[key]
    print "%f fraction of references covered by these contigs is"%(1.*cov/total)
    return

def GC_pos(ref, binsize, outfile): #prints out the GC percentage as function of positions
        n = 0
        bases = []
        for i in range(len(ref)+1):
            n += 1
            try:
                bases.append(ref[i-1])
            except:
                print "error with reference file!"
            if n >= binsize:  #once bin size has reached, write it all out and restart
                outfile.write("%d\t%f\n"% (i+n/2, GCcontent(bases)))
                bases = []
                n = 0    
        return

def GCcontent (seq): #returns the GC content of the input sequence - internal function
    GC = 0
    for c in seq:
        if c == "G" or c == "C" or c == "g" or c == "c":
            GC += 1
    return 1.*GC/len(seq)
    
def consensus(): #outputs a dictionary of {pos: consensus} and an associated 
#confidence score which is fraction consensus base/ total # bases
    consensus = {} #consensus = {pos:(base, score)}
    for key in data:
        for pos in data[key]:
            SNPs = {"a":0, "c":0,"g":0,"t":0}
            for base in data[key][pos]:
                SNPs[base] += 1
            call = keywithmaxval(SNPs)
            score = 1.*SNPs[call]/len(data[key][pos])
            print (call, score)
    return

def keywithmaxval(d):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value"""  
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]
    
            
if __name__ == "__main__":
    path = "./"
    #input_pileup("%sbowtiecontigs.pileup" %path, limit = 0)    
    with open("%s8templates.fasta"%path, "r") as infile:    
        references = input_ref(infile)    
    print "all data imported"
    
    for key in references:    
        with open("%sGC_pos%s.txt"%(path, key), "w") as outfile: 
            GC_pos(references[key], 100, outfile)

    #consensus()
    #percent_cov()
    #output_coverages_ends(path)    
    print "done"
    
