# -*- coding: utf-8 -*-
"""
Created on Fri Dec 05 11:29:49 2014
Lib17 Barcode Sampler/Examiner
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
    
def extract_barcodes(infile, limit, qcutoff):
    count = 0
    passed = 0
    experiment = {}
    #extracts the barcodes into a dictionary from RI read
    print "extracting barcodes: %s" %infile

    for line in infile:  #each line is a new read
        if line[0] != "@":        #ignore the line if it is not a new read
            print "@ not found, check format"
            continue  
                
        barcode = infile.next().strip()
        infile.next() #ignoring the comment line
        quality = infile.next().strip()
        count += 1
        
        #this skips the read if quality score is bad
        lowqual = False        
        for s in quality:
            if Qscore[s] < qcutoff:
                lowqual = True
                break
        if lowqual == True:
            continue
        
        #store barcode and associated quality in the dictionary
        if barcode in experiment:
            experiment[barcode].append(quality)
        else:
            experiment[barcode] = [quality]
        passed += 1

        #this breaks the loop at the set limit of reads
        if limit != 0:        
            if count >= limit:
                break
        #print barcode
    print "%d reads analyzed %d reads passed quality cutoff" % (count, passed)
    return experiment

def input_barlist(infile): #just inputs barcodes from a file listing it
    experiment = {}
    for bar in infile:
        experiment[bar] = 1
    print len(experiment)
    return experiment

def GC_histo(experiment, outfile): #outputs the base compositions for each bargroup for histogram
    for bar in experiment:
        outfile.write("%f\n"%base_composition({bar}))
    print "GC contents printed!"
    return

def base_composition(experiment): #prints out the ACGT composition of barcodes
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
    return 1.0*(G+C)/(A+C+G+T)

def plot_Q_scores(experiment): #plots out the aggregate Qscores for each position
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.pyplot as plot

    data = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]  #contain a list of qscores for each position 
        
    for barcode in experiment:
        for quality in experiment[barcode]:
            for i in range(15):
                data[i].append(Qscore[quality[i]])
    plot.figure(figsize=(20,12))           
    plot.boxplot(data, whis=1, widths = 0.2)
    plot.savefig("%sbarQscores.png"%path, format = 'png')
    return

def bar_overlap(infile1, infile2): #determines how many barcodes overlapped between two barcode lists
    barlist1 = input_barlist(infile1)
    barlist2 = input_barlist(infile2)
    overlaplist = { k: barlist1.get(k, 0) + barlist2.get(k, 0) for k in set(barlist1) | set(barlist2) }
    overlapped = 0
    for bar in overlaplist:
        if overlaplist[bar] > 1:
            overlapped += 1
            print bar
    print "%d number of total barcodes and %d barcodes overlapped"%(len(overlaplist), overlapped)
    return

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

def reads_bar(experiment, outfile): #outputs the barcode vs the number of reads for that group
    for bar in experiment:
        outfile.write("%s\t%d\n"%(bar, len(experiment[bar])))
    print "reads_bar file written!"
    return


if __name__ == "__main__":

    path = "./"    
        
    import timeit
    tic = timeit.default_timer()        

    infile = open("%sRI.fastq" % path)
    experiment = extract_barcodes(infile, 0, 20)

    with open("%sreadsbartotal.txt"%path,"w") as outfile:
        reads_bar(experiment, outfile)    

    find_orphans(experiment, 10)
    print "%d barcodes in experiment"%len(experiment)

    with open("%sreadsbar10.txt"%path,"w") as outfile:
        reads_bar(experiment, outfile)    

    with open("%shamminglist10.txt"%path, "w") as outfile:
        output_hammings(1000, experiment, outfile, True)

    adopt_orphans(experiment, experiment, 1)
    print "%d barcodes in experiment"%len(experiment)

    with open("%shamminglist10adop.txt"%path, "w") as outfile:
        output_hammings(1000, experiment, outfile, True)
    
    with open("%sGChisto.txt"%path, "w") as outfile:
        GC_histo(experiment, outfile)    

    #plot_Q_scores(experiment)
    base_composition(experiment)    

    toc = timeit.default_timer()
    
    print "%d seconds processing time" % (toc-tic)