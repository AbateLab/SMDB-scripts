# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 10:43:42 2015
VCF Parser for SNP call data
@author: Freeman
"""

import os
cwd = os.getcwd()
experiment = []     #a list of all the bargroups
#reference sequence
ref = {"SS":"AATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATGGTCCCGGCGGCGCAGCAGACGGCTACGGCACCGGATGCAGCACTGACGTTCCCGGAAGGCTTTCTGTGGGGCTCGGCGACGGCAAGTTATCAGATTGAAGGTGCAGCAGCAGAAGATGGTCGTACGCCGTCCATCTGGGACACCTACGCGCGTACGCCGGGTCGTGTGCGTAATGGTGATACCGGCGACGTTGCCACGGATCATTATCACCGTTGGCGCGAAGACGTGGCACTGATGGCAGAACTGGGTCTGGGTGCATACCGTTTTAGTCTGGCTTGGCCGCGTATTCAGCCGACCGGTCGTGGTCCGGCACTGCAAAAAGGTCTGGATTTCTATCGTCGCCTGGCGGACGAACTGCTGGCCAAAGGCATCCAGCCGGTCGCCACCCTGTACCATTGGGATCTGCCGCAAGAACTGGAAAATGCAGGCGGTTGGCCGGAACGTGCTACGGCAGAACGCTTTGCGGAATATGCTGCGATTGCCGCAGATGCCCTGGGTGACCGCGTTAAAACCTGGACCACGCTGAACGAACCGTGGTGCAGTGCGTTCCTGGGCTACGGTTCCGGCGTTCACGCACCGGGTCGTACCGATCCGGTCGCTGCGCTGCGCGCCGCACATCACCTGAACCTGGGTCATGGCCTGGCAGTTCAGGCTCTGCGTGATCGTCTGCCGGCAGACGCACAATGTAGCGTCACCCTGAATATTCATCACGTGCGTCCGCTGACGGATTCTGACGCTGATGCGGACGCCGTGCGTCGCATCGATGCACTGGCTAACCGCGTTTTTACCGGTCCGATGCTGCAGGGCGCATATCCGGAAGATCTGGTCAAAGACACCGCTGGTCTGACGGATTGGTCTTTTGTGCGTGATGGTGACCTGCGTCTGGCACACCAAAAACTGGATTTCCTGGGTGTCAACTATTACTCACCGACCCTGGTGTCGGAAGCGGATGGTAGCGGCACGCATAATTCTGACGGTCACGGTCGTAGTGCACATTCCCCGTGGCCGGGTGCAGATCGTGTTGCATTCCACCAGCCGCCGGGTGAAACCACGGCAATGGGCTGGGCTGTCGATCCGAGCGGCCTGTATGAACTGCTGCGTCGCCTGAGCTCTGACTTTCCGGCGCTGCCGCTGGTGATTACCGAAAATGGTGCTGCGTTCCATGATTATGCCGACCCGGAAGGCAACGTTAATGATCCGGAACGTATTGCATACGTTCGTGACCACCTGGCAGCAGTCCATCGTGCAATCAAAGATGGTTCAGACGTGCGCGGCTATTTTCTGTGGTCGCTGCTGGATAACTTCGAATGGGCACATGGTTACAGCAAACGTTTTGGCGCTGTGTATGTTGATTACCCGACCGGCACGCGCATCCCGAAAGCCAGTGCTCGTTGGTATGCTGAAGTTGCTCGCACGGGCGTTCTGCCGACCGCTGGGGATCCGAATTCGAGCTCCGTCGACAAGCTTGCGGCCGCACTCGAGCACCACCACCACCACCACTGAGATCCGGCTGCTAACAAAGCCCGAAAGGAAGCTGAGTTGGCTGCTGCCACCGCTGAGCAATAAC", "TM":"AATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATGCATCACCATCATCATCACATGAACGTGAAAAAGTTCCCTGAAGGATTCCTCTGGGGTGTTGCAACAGCTTCCTACCAGATCGAGGGTTCTCCCCTCGCAGACGGAGCTGGTATGTCTATCTGGCACACCTTCTCCCATACTCCTGGAAATGTAAAGAACGGTGACACGGGAGATGTGGCCTGCGACCACTACAACAGATGGAAAGAGGACATTGAAATCATAGAGAAACTCGGAGTAAAGGCTTACAGATTTTCAATCAGCTGGCCAAGAATACTTCCGGAAGGAACAGGAAGGGTGAATCAGAAAGGACTGGATTTTTACAACAGGATCATAGACACCCTGCTGGAAAAAGGTATCACACCCTTTGTGACCATCTATCACTGGGATCTTCCCTTCGCTCTTCAGCTGAAAGGAGGATGGGCGAACAGAGAAATAGCGGATTGGTTCGCAGAATACTCAAGGGTTCTCTTTGAAAATTTCGGTGATCGTGTGAAGAACTGGATCACCTTGAACGAACCGTGGGTTGTTGCCATAGTGGGGCATCTGTACGGAGTCCACGCTCCTGGAATGAGAGATATTTACGTGGCTTTCCGAGCTGTTCACAATCTCTTGAGGGCACACGCCAGAGCGGTGAAAGTGTTCAGGGAAACCGTGAAAGATGGAAAGATCGGAATAGTTTTCAACAATGGATATTTCGAACCTGCGAGTGAAAAAGAAGAAGACATCAGAGCGGTGAGATTCATGCATCAGTTCAACAACTATCCTCTCTTTCTCAATCCGATCTACAGAGGAGATTACCCGGAGCTCGTTCTGGAATTTGCCAGAGAGTATCTACCGGAGAATTACAAAGATGACATGTCCGAGATACAGGAAAAGATCGACTTTGTTGGATTGAACTATTACTCCGGTCATTTGGTGAAGTTCGATCCAGATGCACCAGCTAAGGTCTCTTTCGTTGAAAGGGATCTTCCAAAAACAGCCATGGGATGGGAGATCGTTCCAGAAGGAATCTACTGGATCCTGAAGAAGGTGAAAGAAGAATACAACCCACCAGAGGTTTACATCACAGAGAATGGGGCTGCTTTTGACGACGTAGTTAGTGAAGATGGAAGAGTTCACGATCAAAACAGAATCGATTATTTGAAGGCCCACATTGGTCAGGCATGGAAGGCCATACAGGAGGGAGTGCCGCTTAAAGGTTACTTCGTCTGGTCGCTCCTCGACAATTTCGAATGGGCAGAGGGATATTCCAAGAGATTTGGTATTGTGTATGTAGACTACAGCACTCAAAAACGCATCGTAAAAGACAGTGGGTACTGGTACTCGAATGTGGTTAAAAACAACGGTCTGGAAGACTGAGAATTCGAGCTCCGTCGACAAGCTTGCGGCCGCACTCGAGCACCACCACCACCACCACTGAGATCCGGCTGCTAACAAAGCCCGAAAGGAAGCTGAGTTGGCTGCTGCCACCGCTGAGCAATAAC"}

class bargroup(object):
    def __init__ (self, barcode = "", SNPs = [], numreads = 0):
        self.barcode = barcode
        self.SNPs = SNPs
        self.numreads = int(numreads)
        return
    def __str__(self):
        return "bargroup %s with %d reads and %d SNPs"%(self.barcode, self.numreads, len(self.SNPs))

class SNP(object):   #this stores information for one SNP
    def __init__ (self, template = "*", pos = 0, allele = "*", quality = "*", genotype = "*", readdepth = 0):
        self.pos = pos
        self.allele = allele
        self.quality = quality
        self.template = template
        self.genotype = genotype
        self.readdepth = readdepth
        return
    
    def __str__(self):
        return "%s\t%s\t%s\t%s\n"%(self.template, self.pos, self.allele, self.quality)
    
def outpickle(exp, outfile):
    import cPickle as pickle
    pickle.dump(exp, outfile)
    print "experiment pickled"
    return

def inpickle(infile):
    import cPickle as pickle
    exp = pickle.load(infile)
    print "experiment imported"
    return exp

#VCF file format:
#SS+TM   116     .       C       T       221.999 .       DP=50;VDB=3.82253e-06;SGB=-0.693147;MQSB=0.999423;MQ0F=0;AF1=1;AC1=2;DP4=0,0,36,13;MQ=40;FQ=-170.988    GT:PL

def in_putVCF(infile): #this one takes a VCF file and returns a list of SNPs out of it
    count = 0
    snplist = []    
    if infile.readline().strip() == "##fileformat=VCFv4.2":    
        for line in infile:
            if line[0] == "#": #skips the header lines
                continue
            line = line.split()
            INFO = line[7].split(";")
            FORMAT = line[-1].split(":")
            genotype = FORMAT[0]
            readdepth = int(INFO[0][3:])

            snplist.append(SNP(line[0], int(line[1]), line[4][0], line[5], genotype, readdepth)) #adds the snp to list
            count += 1
    else:
        print "fileformat error! is it VCFv4.2?"
        print infile.readline().strip()
        print infile
        return
    print "%d snps processed in file %s"%(count, infile)
    return snplist        

def in_putall(path): #this one goes through the directory and inputs all the VCFs
    count = 0
    for f in os.listdir(path):
        if os.path.isdir("%s%s"%(path,f)): #ignore if is a folder
            continue
        elif f.find(".vcf") != -1:  #if is a vcf file: open this file for processing
            #parse the file name
            rev = f.find(".fasta")
            numreads = f[15:rev]
            bar = f[0:15]
            #process the vcf file
            infile = open("%s%s"%(path, f))
            snplist = in_putVCF(infile)
            #make this bargroup and add to the list of bargroups
            experiment.append(bargroup(bar,snplist,numreads))
            count += 1
    print "%d files processed"%count
    return experiment

def pos_vs_qual(outfile): #exports the quality vs position of each snp
    count = 0
    for bar in experiment:
        for SNP in bar.SNPs:
            outfile.write("%s\t%s\n"%(SNP.pos, SNP.quality))
            count +=1 
    print "%d snps written out"%count
    return

def qual_filter(cutoff): #removes SNPs under cutoff qual value
    count = 0
    for bar in experiment:
        newsnps = []
        for SNP in bar.SNPs:
            if float(SNP.quality) > cutoff:
                newsnps.append(SNP)
                count += 1
        bar.SNPs = newsnps           #newnsps only contains pass cutoff SNPS
    print "%d SNPS remained w/ cutoff %d"%(count, cutoff)
    return

def genotype_filter(): #removes all non 1/1 variants (because it is assumed haploid)
    count = 0
    for bar in experiment:
        newsnps = []
        for SNP in bar.SNPs:
            if SNP.genotype == "1/1":
                newsnps.append(SNP)
                count += 1
        bar.SNPs = newsnps
    print "%d SNPS remained after homozygous genotype filter"%(count)
    return
            
def generateconsensus(bar, template = "SS"): #generates a consensus w/ the SNP calls for the bar input

    consensus = list(ref[template])
    numsnps = 0
    #go through each position, if it exists in SNPS, replace w/ the SNP
    for SNP in bar.SNPs:
        if SNP.template == template:
            consensus[int(SNP.pos)-1]=SNP.allele
            numsnps += 1
    bar.consensus = "".join(consensus)
    return bar.consensus, numsnps
    
def in_putnobar(infile): #inputs the vcf for the nobar sample
    experiment.append(bargroup("all", in_putVCF(infile), 9999))
    return experiment

def variants_per_bar(outfile): #outputs the barcode and the number of variants in it
    for bar in experiment:
        outfile.write("%s\t%d\n"%(bar.barcode, len(bar.SNPs)))
    print "variants per bar outputted!"
    return

def qual_readdepth(outfile): #outputs the quality vs the read depth
    for bar in experiment:
        for SNP in bar.SNPs:
            outfile.write("%s\t%d\n"%(SNP.quality,SNP.readdepth))
    print "qual vs readdepth outputted!"
    return
    
def mutation_types(): #outputs the frequenceis of every mutation type
    mutations = {} #types of mutations
    for bar in experiment:
        for SNP in bar.SNPs:
            muttype = ref[SNP.pos-1] + "->" + SNP.allele
            if muttype in mutations:            
                mutations[muttype] += 1
            else:
                mutations[muttype] = 1
    print mutations
    return mutations

def consensus_output(outfile, template = "SS", duplicate = False, multiple = 1): #outputs the consensus sequences w/o duplicates
    consensus = {}
    count = 0
    for bar in experiment:
        seq, numsnps = generateconsensus(bar, template)
        
        if duplicate and numsnps >= multiple:
            outfile.write(">%s\n%s\n"%(template,seq))
            continue
    
        if not duplicate and numsnps >= multiple:              
            if seq in consensus:
                consensus[seq] += 1
            else:
                consensus[seq] = 1
        
    if not duplicate:
      for seq in consensus.keys():
        count += 1
        outfile.write(">%d\n%s\n"%(count,seq))
    print "%d number of unique consensus seqs"%count
    return consensus
    
def co_coverage(outfile, template, threshold = 5): 
#this one outputs the co-coverage information for a heatmap
#for each position in a 2d array of len(seq) vs len(seq), we calculate the number of barcodes
#that contained reads over a certian threshold, covering both of these positions
    import numpy as np
    count = 0
    length = len(ref[template])
    array = np.zeros((length, length))
    for bar in experiment:
        count += 1
        SNPs = [SNP for SNP in bar.SNPs if SNP.template == template] #filter SNPS on template and readdepth
        SNPs = [SNP for SNP in SNPs if SNP.readdepth >= threshold] #filter SNPS on threshold
        for SNPi in SNPs: #now add all these co-cov positions into the array
            for SNPj in SNPs:
                array[SNPi.pos][SNPj.pos] += 1
        #if count >= 100:
            #break
    print "%d array length, number of bars looked at %d"%(len(array), count)
    outpickle(array, outfile)    
    return array
        
if __name__ == "__main__":
    import timeit

    tic = timeit.default_timer()
    
    path = "./"
    in_putall("%sfastas/"%path)
    #with open("%snobar.vcf"%path) as infile:
    #    experiment = in_putnobar(infile)
    #print len(experiment)
    tic = timeit.default_timer()
    #with open("%sco-cov.pickle"%path, "w") as outfile:
    #    co_coverage(outfile,"SS", 5)
    genotype_filter()
    #with open("%sSNP_consensus-multi.fasta"%path, "w") as outfile:
    #    consensus_output(outfile, "SS", False, 2)
    #mutation_types()
    with open("%sSNP_posqualallhomo.txt"%path, "w") as outfile:
        pos_vs_qual(outfile)
    qual_filter(50)
        
    with open("%sSNP_varperbarhomo.txt"%path, "w") as outfile:
        variants_per_bar(outfile)
        
    #with open("%sSNP_qualdepth.txt"%path, "w") as outfile:
    #    qual_readdepth(outfile)
    
    toc = timeit.default_timer()
    print "%d seconds processing time" % (toc-tic)

    