# -*- coding: utf-8 -*-
"""
R4 Parser - for dealing with Barcoded Reads in R4 format for Lib17
Created on Fri May 23 13:46:02 2014

@author: Freeman
"""
Qscore = dict((chr(i),i-33) for i in range(33,90))


def outpickle(exp, outfile):
    import cPickle as pickle
    pickle.dump(exp, outfile)
    print "experiment pickled"
    outfile.close()
    return

def inpickle(infile):
    import cPickle as pickle
    exp = pickle.load(infile)
    print "experiment imported"
    infile.close()
    return exp
    
    
class Read (object):
    '''this reads class holds a paired end read'''
    def __init__ (self, ID = "default", seqA = "none", qualityA = None, seqB = "none", qualityB = None):
        self.ID = ID    #this is a barcode + a number (unique)
        self.seqA = seqA
        self.qualityA = qualityA
        self.seqB = seqB
        self.qualityB = qualityB
        return
        
    def __str__ (self):
        return "ID: %s \t SeqA: %s \t SeqB: %s" % (self.ID, self.seqA, self.seqB)

def in_put (infileA, infileB, limit = 0):      #this one takes R4A and R4B files and puts reads into experiment
    experiment = {}       #new dictionary for reads with barcodes as the key
    count = 0
    
    for line in infileA:  #each line is a new read
        if line[0] != "@":        #ignore the line if it is not a new read
            print "@ not found, check format"
            continue  
        
        count += 1        
        ID = line.strip()[:18]          #the first 18 chars, the last char indicates read pair        
        seqA = infileA.next().strip()
        infileA.next() #ignoring the comment line
        #infileA.next()        
        #qualityA = None
        qualityA = infileA.next().strip()
        
        if infileB.next().strip()[:18] == ID: #check synchronization w/ 2 files
            seqB = infileB.next().strip()
            infileB.next()
            #infileB.next()    
            #qualityB = None
            qualityB = infileB.next().strip()
        else: 
            print "infileA and B desynchronized!"
            return
        
        #now assign each segment of the split into a read
        read = Read (ID = ID, seqA = seqA, qualityA = qualityA, seqB = seqB, qualityB = qualityB)
        #add to the dictionary

        if read.ID[1:16] in experiment:       #ID[1:16] is the barcode sequence
            experiment[read.ID[1:16]].append(read)
        else:
            experiment[read.ID[1:16]] = [read]
        
        if limit != 0:
            if count >= limit:
                break
    
    infileA.close()
    infileB.close()
    print "%d reads inputted" % count
    return experiment

def find_orphans(experiment, cutoff = 1):     #identifies the barcodes that only have x reads as orphans
    orphans = {}          #clears out the orphans dictionary    
    for key in experiment:
        if len(experiment[key]) <= cutoff:
            orphans[key] = experiment[key]
    for key in orphans:
        experiment.pop(key)
    
    print "number of orphans discovered %d" %len(orphans)
    return orphans
  
def adopt_orphans(parents, orphans, dist):
    adopted = 0
    count = 0
    total = len(orphans)    
    for bar in orphans.keys():
        for parent in parents.keys():
            ham = hamming_distance(bar, parent)
            if  ham <= dist and ham != 0:
                parents[parent] = parents[parent] + orphans[bar]
                adopted += 1
                orphans.pop(bar)
                break
        count = count + 1
        #progress indicator        
        if count % 1000 == 0:
            print 1.00*count/total
 
    print "%d adopted!"%adopted
    return parents

def read_counts_histo (experiment, outfile, count_reads = False):  #prints into a file the # of reads for each barcode to make histogram
    for barcode in experiment:
        if count_reads == False:
            outfile.write("%d\n"%len(experiment[barcode]))
        elif count_reads == True:     #this prints out a data point for each read
            for read in experiment[barcode]:
                outfile.write("%d\n"%len(experiment[barcode]))
    return

def export_fasta(outpath, experiment):  #takes an experiment, and exports the barcode groups as separate fasta files for each group
    for key in experiment:
        count = 0
        outfile = open("%s%s%s.fasta" % (outpath, key, str(len(experiment[key]))), "w")
        for read in experiment[key]:
            outfile.write("%s:1-%d\n%s\n"%(read.ID, count, read.seqA))
            outfile.write("%s:2-%d\n%s\n"%(read.ID, count, read.seqB))
            count += 1
        outfile.close()
    print "fasta exported"
    return

def export_fastq(outpath, experiment):  #takes an experiment, and exports the barcode groups as separate fasta files for each group
    for key in experiment:
        count = 0
        outfile = open("%s%s%s.fastq" % (outpath, key, str(len(experiment[key]))), "w")
        for read in experiment[key]:
            outfile.write("%s:1-%d\n%s\n+\n%s\n"%(read.ID, count, read.seqA, read.qualityA))
            outfile.write("%s:2-%d\n%s\n+\n%s\n"%(read.ID, count, read.seqB, read.qualityB))
            count += 1
        outfile.close()
    print "fastq exported"
    return

def export_fastq_twofiles(outfile1, outfile2, experiment): #exports reads all into two files
    for key in experiment:
        count = 0
        for read in experiment[key]:
            outfile1.write("%s:1-%d\n%s\n+\n%s\n"%(read.ID, count, read.seqA, read.qualityA))
            outfile2.write("%s:2-%d\n%s\n+\n%s\n"%(read.ID, count, read.seqB, read.qualityB))
            count += 1
    print "fastq exported into two files"
    return
    
def export_fasta_onefile(outfile, experiment): #exports reads all into one file
    for key in experiment:
        count = 0
        for read in experiment[key]:
            outfile.write("%s:1-%d\n%s\n"%(read.ID, count, read.seqA))
            outfile.write("%s:2-%d\n%s\n"%(read.ID, count, read.seqB))
            count += 1
    outfile.close()
    print "fasta exported"
    return
    
def barcomposition(experiment): #prints out the ACGT composition of barcodes
    A, C, G, T = 0, 0, 0, 0
    for bar in experiment:
        for letter in list(bar):
            if letter == "A":
                A += 1
            elif letter == "C":
                C += 1
            elif letter == "G":
                G += 1
            elif letter == "T":
                T += 1
    
    print "Barcode composition A: %s" % (1.0*A/(A+C+G+T))
    print "Barcode composition C: %s" % (1.0*C/(A+C+G+T))
    print "Barcode composition G: %s" % (1.0*G/(A+C+G+T))  
    print "Barcode composition T: %s" % (1.0*T/(A+C+G+T))
    return

def output_bar(experiment, outfile): #prints out each unique barcode to a file
    for bar in experiment:
        outfile.write("%s\n" %bar)
    print "barcodes written!"
    outfile.close()
    return

def remove_PCR_duplicates(experiment): #removes all the reads that are exactly the same in each bargroup    
    count = 0
    uniques = 0
    for bar in experiment.keys():
        existing_seqs = {} #stores all the seqs seen so far for this bargroup - dictionary format is faster lookup
        for read in experiment[bar]:
            seq = read.seqA + "-" + read.seqB
            if seq in existing_seqs:
                experiment[bar].remove(read)
                count += 1
            else:
                existing_seqs[seq] = 1
                uniques += 1
    print "%d pcr duplicate reads removed!"%count
    return experiment
    
#internal function to calculate hamming distances                
def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    
def get_hammings(barcode, experiment): #return a list of hamming dist of the barcode
#vs every barcode in the experiment
    hammings = []
    for bar in experiment:
        hammings.append(hamming_distance(barcode, bar))
    hammings.remove(0)          #theres always a 0 value corresponding to hamming vs itself
    return hammings
    
def output_hammings(num, experiment, outfile, lowest = False): #returns the hamming dists in a list 
#for num number of subsampld bargroups
    count = 0
    for bar in experiment:
        hammings = get_hammings(bar, experiment)
        if lowest == False:        
            for dist in hammings:
                outfile.write("%s\n"%dist)
        else:
            outfile.write("%s\n"%min(hammings))
        count += 1
        if count >= num:
            break
    return

    
def mask_lowQ(experiment, cutoff, truncate = False): #this function masks low quality bases
   count = 0
 
   for bar in experiment:
        for read in experiment[bar]:
            for i in range (len(read.seqA)-1):
                if Qscore[read.qualityA[i]] < cutoff:
                    count += 1                    
                    if truncate:   #break the looping if we truncate
                        read.seqA = read.seqA[:i]
                        read.qualityA = read.qualityA[:i]
                        break
                    else:
                        read.seqA = read.seqA[:i] + "n" + read.seqA[i+1:]
            for i in range (len(read.seqB)-1):
                if Qscore[read.qualityB[i]] < cutoff:
                    count += 1                    
                    if truncate:
                        read.seqB = read.seqB[:i]
                        read.qualityB = read.qualityB[:i]
                        #print read.seqB
                        break
                    else:
                        read.seqB = read.seqB[:i] + "n" + read.seqB[i+1:]    
   print "%d low quality bases masked as n or truncated" %count
   return
        


def reads_bar(experiment, outfile, unique = False): #outputs the barcode vs the number of reads for that group
    if unique:
        remove_PCR_duplicates(experiment)
    for bar in experiment:
        outfile.write("%s\t%d\n"%(bar, len(experiment[bar])))
    print "reads_bar file written!"
    return

def rarefaction(inpathA, inpathB, datapoints, threshold, outfile, unique=True): #this does a rarefaction curve of readgroups over x reads
    for p in datapoints:
        with open(inpathA, "r") as infileA, open(inpathB, "r") as infileB:
            experiment = in_put(infileA, infileB, limit = p)
        if unique:
            remove_PCR_duplicates(experiment)
        #mask_lowQ(experiment, 10, truncate = False) 
        find_orphans(experiment, threshold)
        outfile.write("%d\t%d\n"%(p,len(experiment)))
    print "rarefaction data exported!"
    return
    
if __name__ == "__main__":
  
    import timeit
    tic = timeit.default_timer()
    
    path = "./"

    #infileA = open("%sR4A.fastq" %path)
    #infileB = open("%sR4B.fastq" %path)
    #experiment = in_put(infileA, infileB, limit = 0)    
    #print "%d barcodes in experiment"%len(experiment)
    #remove_PCR_duplicates(experiment)
   # outpickle(experiment, open("%sexperiment-nodups.pickle"%path, "w"))  
   # with open("%sreadsbar.txt" %path, "w") as outfile:
   #     reads_bar(experiment, outfile, unique = False)
   # with open("%sreadscounts.txt"%path, "w") as outfile:
   #     read_counts_histo(experiment, outfile, True)
    #find_orphans(experiment, 500)
    #print "%d barcodes in experiment"%len(experiment)
    #adopt_orphans(experiment, experiment, 1)
              

    #with open("%sreadsbar500.txt" %path, "w") as outfile:
    #    reads_bar(experiment, outfile, unique = False)
    #with open("%sreadsbar500unique.txt" %path, "w") as outfile:
    #    reads_bar(experiment, outfile, unique = True)

    #mask_lowQ(experiment, 20, truncate = False)
    #export_fasta("%sfastas/"%path, experiment)
    #with open("%sexperiment.fasta"%path, "w") as outfile:   
    #    export_fasta_onefile(outfile, experiment)

    with open("%srarefaction.txt"%path, "w") as outfile:
        rarefaction("%sR4A.fastq"%path, "%sR4B.fastq"%path, (1e4, 2e4, 4e4, 8e4, 2e5, 4e5, 8e5, 2e6, 4e6, 8e6), 100, outfile, False)   
    

  
    toc = timeit.default_timer()
    print "%d seconds processing time" % (toc-tic)
    