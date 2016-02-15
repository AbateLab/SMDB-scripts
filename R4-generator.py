# -*- coding: utf-8 -*-
"""
Created on Fri Dec 05 17:38:34 2014
R4generator for Lib17 -> two options - generate R4A R4B files or generate R4 folder

Created on Thu Sep 18 17:35:56 2014
*********oh god, this is a VERY fragile program*********** Note the assumptions!
R1-Barcode mapper V2 - it maps the barcode onto the R1 R2 reads as its ID
This one is much faster, but it assumes the reads in R1 R2 RI are in order, and if there is one mismatch, it simply skips it
@author: Freeman
"""

Qscore = dict((chr(i),i-33) for i in range(33,90))
import os

def mapit (R1, R2, RI, outfileA, outfileB, barcutoff = 30, limit = 0): #accepts a path, opens that file and loads contents 
        #barcutoff is the quality score at which we cut off if any base in the read is below that score
        num = 0
        count = 0
        filtered = 0        #number of reads that failed filter
        passed = 0          #number of reads that passed filter
        desynch = 0
        
        try:
            for line in RI:     #this is for getting the barcode
                num += 1
                passfilter = True #resets pass filter for determining if we should skip this read
                
                if limit != 0:
                    if num > limit:
                        break

                line = line.strip()   #get rid of 
                if line[0] != "@":    #ignores the line if it doesn't start with an @ character
                    print "RI: line ignored because no @"    
                    print line
                    continue
                
                else:    #look in bar-read for barcode
                    ID = line #just the clusterID
                    bar = RI.next().strip()   #next line is the sequence, we extract barcode
                    RI.next()                  #ignore next line of comments
                    barscore = RI.next().strip()  #adds the quality scores lines
                    #check barcode quality score, if fail, then go on to next read                   
                             
                    for s in barscore:
                        if Qscore[s] < barcutoff:
                            passfilter = False
                            filtered = filtered + 1
                            break
                    else:
                        passed = passed + 1
                        
                if passfilter == False:             #forget this read if it does not pass filter AND skip the read in the Sequence file
                    continue                    
                
                else:   #if pass quality filter, then we look to the actual read to see if we get the same cluster
                    for line in R1:
                        if line[0] != "@":    #ignores the line if it doesn't start with an @ character
                            print "R1: line ignored because no @"  
                            print line
                            continue
                        
                        if line[21:44] == ID[21:44]:     #check we are synchronized in both files
                            #go ahead and store in the sequence data                        
                            seqA = R1.next().strip()
                            R1.next()                   #skip the comment line
                            scoreA = R1.next().strip()
                            outfileA.write("@%s-%d:1\n%s\n+\n%s\n" % (bar, count, seqA, scoreA))                            
                            count += 0.5                            
                            break
                        else: 
                            desynch += 0.5
                            R1.next()
                            R1.next()
                            R1.next()
                        
                    for line in R2:
                        if line[0] != "@":    #ignores the line if it doesn't start with an @ character
                            print "R2: line ignored because no @"                    
                            print line                            
                            continue
                        
                        if line[21:44] == ID[21:44]:     #check we are synchronized in both files
                            #go ahead and store in the sequence data                        
                            seqB = R2.next().strip()
                            R2.next()                   #skip the comment line
                            scoreB = R2.next().strip()
                            outfileB.write("@%s-%d:2\n%s\n+\n%s\n" % (bar, count, seqB, scoreB))
                            count += 0.5   #now we record number of reads written
                    
                            break    #if we find the read correctly then we break and go on to the next bar read                            
                        else: 
                            desynch += 0.5
                            R2.next()
                            R2.next()
                            R2.next()                

            
            print "%d number of barcode reads looked at" % num
            print "%d number of reads outputted" % count      
            print "%d number of bars passed filter" % passed
            print "%d number of bars failed filter" % filtered
            print "%d number of times desynchronized" % desynch            
            R1.close()
            R2.close()            
            RI.close()
            outfileA.close()
            outfileB.close()
            return

        except (IOError):
            #print num
            print "bad file path"
            return "bad file path!"

#this one creates a new fasta file for each barcode
def map_separate (R1, R2, RI, barcutoff = 30, limit = 0): #accepts a path, opens that file and loads contents 
        #barcutoff is the quality score at which we cut off if any base in the read is below that score
        num = 0
        count = 0
        filtered = 0        #number of reads that failed filter
        passed = 0          #number of reads that passed filter
        desynch = 0
        if os.listdir("./R4/") != []:  #R4 folder must be empty to start
            print os.listdir("./R4/")            
            print "*****R4 folder is NOT EMPTY!!!**** Exiting..."
            return            
            
        try:
            for line in RI:     #this is for getting the barcode
                num += 1
                passfilter = True #resets pass filter for determining if we should skip this read
                
                if limit != 0:
                    if num > limit:
                        break

                line = line.strip()   #get rid of 
                if line[0] != "@":    #ignores the line if it doesn't start with an @ character
                    print "RI: line ignored because no @"    
                    print line
                    continue
                
                else:    #look in bar-read for barcode
                    ID = line #just the clusterID
                    bar = RI.next().strip()   #next line is the sequence, we extract barcode
                    RI.next()                  #ignore next line of comments
                    barscore = RI.next().strip()  #adds the quality scores lines
                    #check barcode quality score, if fail, then go on to next read                   
                             
                    for s in barscore:
                        if Qscore[s] < barcutoff:
                            passfilter = False
                            filtered = filtered + 1
                            break
                    else:
                        passed = passed + 1
                        
                if passfilter == False:             #forget this read if it does not pass filter AND skip the read in the Sequence file
                    continue                    
                
                else:   #if pass quality filter, then we look to the actual read to see if we get the same cluster
                    for line in R1:
                        if line[0] != "@":    #ignores the line if it doesn't start with an @ character
                            print "R1: line ignored because no @"  
                            print line
                            continue
                        
                        if line[21:44] == ID[21:44]:     #check we are synchronized in both files
                            #go ahead and store in the sequence data                        
                            seqA = R1.next().strip()
                            R1.next()                   #skip the comment line
                            scoreA = R1.next().strip()
                            write_read("./R4/", bar, "@%s-%d:1\n%s\n+\n%s\n" % (bar, count, seqA, scoreA))         
                            count += 0.5                            
                            break
                        else: 
                            desynch += 0.5
                            R1.next()
                            R1.next()
                            R1.next()
                        
                    for line in R2:
                        if line[0] != "@":    #ignores the line if it doesn't start with an @ character
                            print "R2: line ignored because no @"                    
                            print line                            
                            continue
                        
                        if line[21:44] == ID[21:44]:     #check we are synchronized in both files
                            #go ahead and store in the sequence data                        
                            seqB = R2.next().strip()[66:]
                            R2.next()                   #skip the comment line
                            scoreB = R2.next().strip()[66:]
                            write_read("./R4/", bar, "@%s-%d:2\n%s\n+\n%s\n" % (bar, count, seqB, scoreB))
                            count += 0.5   #now we record number of reads written
                    
                            break    #if we find the read correctly then we break and go on to the next bar read                            
                        else: 
                            desynch += 0.5
                            R2.next()
                            R2.next()
                            R2.next()                

            
            print "%d number of barcode reads looked at" % num
            print "%d number of reads outputted" % count      
            print "%d number of reads passed filter" % passed
            print "%d number of reads failed filter" % filtered
            print "%d number of times desynchronized" % desynch            
            R1.close()
            R2.close()            
            RI.close()
            return

        except (IOError):
            #print num
            print "bad file path"
            return "bad file path!"    
            
#this is an internal function that writes the read to harddrive based on the bar
def write_read(path, bar, string): 
    #search for the bar file, if it exists, it appends to the file, if not it creates one
    if "%s.fastq"%bar in os.listdir(path):
        with open("%s/%s.fastq"%(path,bar), "a") as myfile:
            myfile.write("%s"%string)  
    else:
        with open("%s/%s.fastq"%(path,bar), "w") as myfile:
            myfile.write("%s"%string)  
    return
    
if __name__ == "__main__":
    import timeit
    tic = timeit.default_timer()
    
    path = "./"    
    
    R1 = open("%sR1.fastq" % path)
    R2 = open("%sR2.fastq" % path)
    RI = open("%sRI.fastq" % path)
    R4A = open("%sR4A.fastq"%path, "w")
    R4B = open("%sR4B.fastq"%path, "w")
    mapit(R1, R2, RI, R4A, R4B, 20, 0)    
    
    #map_separate(R1, R2, RI, limit = 10, barcutoff = 30)
    
    toc = timeit.default_timer()
    print "%d seconds processing time" % (toc-tic)