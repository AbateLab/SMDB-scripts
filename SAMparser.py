# -*- coding: utf-8 -*-
"""
SAM parser for Lib17 Data
Created on Fri May 23 13:46:02 2014
@author: Freeman
"""

#SAM Analyzer
#For storing SAM data and determining mapping correctness

experiment = {}
orphans = {}

import cPickle as pickle


def outpickle(exp, outfile):
    pickle.dump(exp, outfile)
    print "experiment pickled"
    outfile.close()
    return

def inpickle(infile):
    exp = pickle.load(infile)
    print "experiment imported"
    infile.close()
    return exp
    
    
class Read (object):
    '''this reads class is specific to SAM files, contains a bunch of attribtues:
ID(barcode), Sequence, Mapping position, Mapping quality, CIGAR string, and FLAG'''
    def __init__ (self, ID = "default", seq = "none", MapPos = -1, MapQ = -1, CIGAR = "none", FLAG = "none", align = "none"):
        self.ID = ID
        self.seq = seq
        self.MapPos = MapPos
        self.MapQ = MapQ
        self.CIGAR = CIGAR
        self.FLAG = FLAG
        self.align = align
        return
        
    def __str__ (self):
        return "ID: %s \t Seq: %s \t MapPos: %d \t CIGAR: %s \t FLAG: %s \t" % (self.ID, self.seq, self.MapPos, self.CIGAR, self.FLAG)

def input_ref(infile):
    reference = {}
    for line in infile:
        if line[0] == ">":
            name = line[1:].strip()
            seq = infile.next().strip()
            reference[name] = seq
    return reference
    
def in_put1 (infile, limit = 0):      #this one takes only one file as input
    reads = {}       #new dictionary for reads with barcodes as the key
    count = 0
    
    for line in infile:  #each line is a new read
        if line[0] == "@":        #ignore the line if it is a header
            continue       
        
        count += 1        

        split = line.split("\t")    #split line into separate words
        #now assign each segment of the split into a read
        try:
            read = Read (ID = split[0][:15], align = split[2], seq = split[9], MapPos = int(split[3]), MapQ = split[4], CIGAR = split[5], FLAG = split[1])
        except:
            print split[3]
            print "read could not be stored"
        
        if read.ID in reads:
            reads[read.ID].append(read)
        else:
            reads[read.ID] = [read]
        
        if limit != 0:
            if count >= limit:
                break
    print "%d reads inputted" % count
    return reads

def map_positions (experiment, path = None, freq = False): #returns a histogram of positions and # of reads that match, and prints it out if path given
    hist = {}
    
    #generate frequency histogram
    for key in experiment:
            for read in experiment[key]:
                pos = read.MapPos
                if pos in hist:
                    hist[pos] = hist[pos] + 1
                else:
                    hist[pos] = 1
 
   #print out into a file   
    if path != None:
            outfile = open(path, "w")
            if freq == True:           #this is the frequency table format
                    for key in hist:
                        outfile.write("%s\t" % key)
                        print >> outfile, hist[key]
                        outfile.write("\n")
            if freq == False:           #this just prints out every map position once for each read
                for key in experiment:
                    for read in experiment[key]:
                        outfile.write("%d\n" %read.MapPos)
            print "Outfile printed to %s" % path
    return hist


def remove_unmapped (experiment): #returns an experiment with all the unmmapped reads removed
    mapped = {}
    for key in experiment:
        for read in experiment[key]:
            if int(read.MapPos) > 0 and int(read.MapQ) >= 20:
                if read.ID in mapped:
                    mapped[read.ID].append(read)
                else:
                    mapped[read.ID] = [read]
    print "reminder: remove unmapped returns a new experiment, does not edit original!"
    return mapped

def remove_repeats(experiment): #looks at reads in the same barcode, and removes the reads that map to the same start position
    num = 0
    for key in experiment:
        
        positions = []
        remove = []
        for read in experiment[key]:
            if read.MapPos in positions:
                num = num + 1
                remove.append(read)
            else:
                positions.append(read.MapPos)
        for read in remove:
            experiment[key].remove(read)
    print "%d reads removed" %num
    return experiment

def mapping_identities(experiment): #prints out the references mapped in a table
    histo = {"*":0, "A":0, "B":0, "C":0, "D":0, "E":0, "F":0, "G":0, "H":0, "I":0, "J":0, "L":0} #used to keep track of frequencies of mapping for each reference
    for bar in experiment:
        for read in experiment[bar]:
            histo[read.align] += 1
    for key in histo:
        print "%s - %d" %(key, histo[key])
    return

def group_identities(experiment, outfile, cutoff = 0.1): #determines how many different templates were in each bargroup    
    #only if a template has cutoff percent of reads does it count as a template in the group    
    for bar in experiment:
        count = 0
        histo = mapping_histo(experiment[bar])
        basecutoff = cutoff*len(experiment[bar])  #determine the cutoff for having enough reads        
        for key in histo:      #counts the # of keys that have >0 reads, not counting unmapped reads
            if histo[key] > basecutoff and key != "*":
                count += 1
        outfile.write("%s\t%d\t%d\n" % (bar, len(experiment[bar]), count)) #barcode, # of reads, templates count
    return
    
def mapping_histo(bargroup): #internal function for generating mapping histo for a set of reads
        histo = {"*":0, "A":0, "B":0, "C":0, "D":0, "E":0, "F":0, "G":0, "H":0, "I":0, "J":0, "L":0} #used to keep track of frequencies of mapping for each reference      
        for read in bargroup:
            histo[read.align] += 1
        return histo

#this is an internal function, not directly called
def purity(bargroup): #returns purity score defined as % of reads mapping to the major template
    histo = mapping_histo(bargroup)
    maxval = 0
    maxkey = "none"
    for key in histo:   #determine the highest read templates
        if histo[key] > maxval:
            maxkey = key
            if key != "*":    #only count it as a valid maxval if it is not unmapped reads max
                maxval = histo[key]
                
    score = 1.0*maxval/len(bargroup)
    return (maxkey, score)


def purity_histo(experiment, outfile): #outputs the histogram of purity scores in the experiment
    #also prints out the frequencies of barcode groups by major templates
    histo = {"*":0, "A":0, "B":0, "C":0, "D":0, "E":0, "F":0, "G":0, "H":0, "I":0, "J":0, "L":0} #used to keep track of frequencies of mapping for each reference              
    for key in experiment:
        maxkey, score = purity(experiment[key])
        histo[maxkey] += 1
        outfile.write("%s \t %s \n"%(maxkey, score))
    print "Template counts by bar group"
    for key in histo:
        print "%s - %d" %(key, histo[key])    
    return

    ''' here for reference, defunct function
def purity_old (experiment, outfile, XY = False, extract = 0, Pvalue = False):   #prints into a file the ratio % of ss or tm mapping reads for each barcode
    #the 2D option generates X: Size of barcode groups Y: Purity of said group
    extracted = {}      #for extraction    
    if Pvalue == True:
        print "Writing Pvalues!"
        
    for barcode in experiment:
        tm = 0
        ss = 0
        for read in experiment[barcode]:
            if read.align == "tm":
                tm = tm + 1
            elif read.align == "ss":
                ss = ss + 1
        if ss == 0 and tm == 0:
            print "unmapped reads only"
            continue
        ratio = max(ss, tm)*100./(ss+tm)
        
        if extract != 0:      #extract ratios > x into a separate experiment
            if ratio >= extract:
                extracted[barcode] = experiment[barcode]
            continue
        
        if Pvalue == True:
            import scipy.stats as stats
            p = stats.binom.cdf(min(ss,tm), ss + tm , 0.5)
            outfile.write("%s\n" %p)
        elif XY == False:
            outfile.write("%s\n" % ratio)
        elif XY == True:
            outfile.write("%s\t%s\n" % (len(experiment[barcode]), ratio))

    outfile.close()
    return extracted
    '''

def read_counts_histo (experiment, outfile, count_reads = False):  #prints into a file the # of reads for each barcode to make histogram
    for barcode in experiment:
        if count_reads == False:
            outfile.write("%d\n"%len(experiment[barcode]))
        elif count_reads == True:
            for i in range(len(experiment[barcode])-1):
                outfile.write("%d\n"%len(experiment[barcode]))
    return

def export_fasta(outpath, experiment):  #takes an experiment, and exports the barcode groups as separate fasta files for each group
    for key in experiment:
        count = 0
        outfile = open("%s%s%s.fasta" % (outpath, key, str(len(experiment[key]))), "w")
        for read in experiment[key]:
            outfile.write("@%s-%d\n%s\n"%(read.ID, count, read.seq))
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

def find_orphans(experiment, cutoff = 1):     #identifies the barcodes that only have x reads as orphans
    orphans = {}          #clears out the orphans dictionary    
    for key in experiment:
        if len(experiment[key]) <= cutoff:
            orphans[key] = experiment[key]
    for key in orphans:
        experiment.pop(key)
    
    print "number of orphans discovered %d" %len(orphans)
    return orphans

def calculate_entropy(experiment, outfile, subsample = False, addbar = False): #takes an experiment, and calculates its entropies and outputs into a file
#entropy is S = Sum(PiLn(Pi)) where Pi is the probability for each bin
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.pyplot as plot

    if subsample: 
        print "subsample mode!"
    entropy_dict = {}    #to store barcode:entropy and returns it    
    for bar in experiment:
        #histo to keep track of mappositions for each template
        histo = {"*":[], "A":[], "B":[], "C":[], "D":[], "E":[], "F":[], "G":[], "H":[], "I":[], "J":[], "L":[]} #used to keep track of mapping for each reference              
        lengths = {"A":5000, "B":5315, "C":5227, "D":5001, "E":5002, "F":5143, "G":3185, "H":4774, "I":4876, "J":3920, "L":4361}
    
        for read in experiment[bar]:
                histo[read.align].append(read.MapPos)  #add the MapPos to the correct location   
            
        entropies = {}
        numreads = {}
        for key in histo:
            if key == "A" or key == "*":    #ignore these two templates they dont'exist
                continue
            S = matplotlib.mlab.entropy(histo[key], range(0, lengths[key], 50))
            entropies[key] = S
            numreads[key]=(len(histo[key]))
        #set the true entropy to the highest of all the templates
        key = keywithmaxval(entropies)
        S = entropies[key]
        numreads = numreads[key]
        entropy_dict[bar] = S
        
        if S < 0:        
            print S
            for read in experiment[bar]:
                print read.MapPos
                
        #subsampling a bargroup
        if subsample:
            if 8<S<9 and 1e4<numreads<1e6:
                print "subsample found! %s entropy, %s reads"%(S, numreads)                
                export_coverage(experiment[bar], key, outfile)
                break
        elif addbar == True:    
            outfile.write("%s\t %d \t %s \t %f\n"%(bar, numreads, key, S))
        else:
            outfile.write("%f\n"%(S))
    outfile.close()
    return entropy_dict

def keywithmaxval(d):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value"""  
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]
     
def export_coverage(reads, key, outfile, GC = False, binsize = 100, references = None):  #takes a list of reads and exports their coverage as histogram
    lengths = {"A":5000, "B":5315, "C":5227, "D":5001, "E":5002, "F":5143, "G":3185, "H":4774, "I":4876, "J":3920, "L":4361}
    coverage = []
    unusedread = 0
    print key
    for i in range(lengths[key]+1):
        coverage.append(0)   #initialize coverage array

    for read in reads:
        if read.align == key and read.MapPos > 0:
            for i in range(len(read.seq)):
                try:
                    coverage[read.MapPos+i-1] += 1
                except:
                    print "read out of range of template"
        else:
            unusedread += 1
    if GC == False:      
        for i in range(len(coverage)):
            outfile.write("%d\t%d\n"%(i, coverage[i]))
        print "%d unused reads"%unusedread
    else:  #put data into bins of x basepairs and output GC content vs coverage in bin
        #putting data into bins
        n = 0
        bases = []
        bincov = 0
        for i in range(len(coverage)):
            n += 1
            bincov += coverage[i-1]
            try:
                bases.append(references[key][i-1])
            except:
                print "error with reference file!"
            if n >= binsize:  #once bin size has reached, write it all out and restart
                outfile.write("%f\t%d\n"% (GCcontent(bases), bincov))
                bases = []
                bincov = 0
                n = 0
    return

def GCcontent (seq): #returns the GC content of the input sequence - internal function
    GC = 0
    for c in seq:
        if c == "G" or c == "C" or c == "g" or c == "c":
            GC += 1
    return 1.*GC/len(seq)
    
def rarefaction (inpath, datapoints, outfile): #this does a rarefaction curve of readgroups over 50unique mapping reads
   
    for p in datapoints:
        with open(inpath, "r") as infile:
            experiment = in_put1(infile, limit = p)
        remove_repeats(experiment)
        experiment = remove_unmapped(experiment)
        find_orphans(experiment,50)
        outfile.write("%d\t%d\n"%(p,len(experiment)))
    print "rarefaction data exported!"
    return
        
if __name__ == "__main__":
  
    import timeit
    tic = timeit.default_timer()
    path = "./"
    with open("%s8templates.fasta"%path, "r") as infile:
        reference = input_ref(infile)
 
    infile = open("%sR4.SAM" %path)
    experiment = in_put1(infile, limit = 0)
 
    print "%d barcodes in experiment"%len(experiment)
    find_orphans(experiment, 100)
    print "%d barcodes in experiment"%len(experiment)
    #mapping_identities(experiment)
    #with open("%sgroupidentities500.txt"%path, "w") as outfile:
    #    group_identities(experiment, outfile, cutoff = 0.1)
    #with open("%spurity1000.txt"%path, "w") as outfile:        
    #    purity_histo(experiment, outfile)
    with open("%sentropies.txt"%path, "w") as outfile:    
        entropy_dict = calculate_entropy(experiment, outfile, subsample = False, addbar = True)
    #with open("%scovS8-9.txt"%path, "w") as outfile:    
    #     entropy_dict = calculate_entropy(experiment, outfile, subsample = True, addbar = True)
    #with open("%sentropy_dict.pickle"%path, "w") as outfile:
    #    outpickle(entropy_dict, outfile)
   # outfile = open("%srarefaction"%path, "w")
    #rarefaction("%sR4.SAM"%path, (1e4, 2e4, 3e4, 4e4, 5e4, 1e5, 2e5, 3e5, 4e5, 5e5, 1e6, 2e6, 3e6, 4e6, 5e6), outfile)   
        
    toc = timeit.default_timer()
    print "%d seconds processing time" % (toc-tic)
    