# -*- coding: utf-8 -*-
"""
Created on Tue May 05 10:53:04 2015

Non-figures plotter for Lib17

@author: Freeman
"""

#histogram plotter for barcodes and reads, for DOLMOS paper

import matplotlib.pyplot as plot
import numpy as np


def in_put(path):
    data = []
    infile = open(path, "r")
    for line in infile:
             data.append(float(line))
    return data

def in_put2(path):
    X = []
    Y = []
    infile = open(path, "r")
    for line in infile:
        X.append(line.split()[0])
        Y.append(line.split()[1])
    return (X,Y)
    
def in_put3(path):
    X = []
    Y = []
    Z = []
    infile = open(path, "r")
    for line in infile:
        X.append(line.split()[0])
        Y.append(line.split()[1])
        Z.append(line.split()[2])
    return (X,Y,Z)

def in_put4(path):
    W = []    
    X = []
    Y = []
    Z = []
    infile = open(path, "r")
    for line in infile:
        W.append(line.split()[0])
        X.append(line.split()[1])
        Y.append(line.split()[2])
        Z.append(line.split()[3])
    return (W,X,Y,Z)
    
    
def templatecounts(path):
    X,Y,Z = in_put3("%sgroupidentities500.txt"%path)
    Y = [int(i) for i in Y]
    Z = [int(i) for i in Z]
    n, bins, patches = plot.hist(Z, bins=9, normed = True, range = (0,9), align = 'left', color = ("#82a53c"), histtype = "bar")
    plot.xticks(range(0,9)) 
    plot.xlabel("number of molecules in barcode group")
    plot.ylabel("fraction of barcode groups")
    plot.savefig("%stemplatecounts500.svg"%path, format = "svg")
    return

def purity(path):
    X,Y = in_put2("%spurity500.txt"%path)
    Y = [float(i)*100 for i in Y]
    
    count = 0    
    for i in Y:
        if i >= 90:
            count += 1
    print "The fraction of bars over 0.9 purity is: %f"%(1.*count/len(Y))
    fig, ax1 = plot.subplots()
    ax2 = ax1.twinx()
    n, bins, patches = ax1.hist(Y, bins=100, normed = True, range = (0,100), cumulative = False, align = 'left', color = ("#82a53c"), histtype = "stepfilled")
    n, bins, patches = ax2.hist(Y, bins=100, normed = True, range = (0,100), cumulative = True, align = 'left', color = ("#82a53c"), histtype = "step")    
    plot.xlim(0,100)  
    plot.title("Purity Histogram") 
    plot.xlabel("percent purity")
    plot.show()
    plot.savefig("%spurityhisto500.svg"%path, format = 'svg')
    return

def purity_by_bar(path):
    X,Y = in_put2("%spurity.txt"%path)
    Y = [float(i)*100 for i in Y]
    Y.sort()
    plot.bar(range(len(Y)), Y, color = ("#82a53c"), lw = 0, width = 1)
    plot.xlim(0,len(Y))  
    plot.title("Percent Major Template") 
    plot.xlabel("Barcode Cluster")
    plot.ylabel("Percent Purity")
    plot.savefig("%spuritybybar.svg"%path, format = "svg")
    return
    
def coverage(path):
    templates = ["Ccov", "Dcov", "Ecov", "Gcov", "Hcov", "Icov", "Jcov", "Lcov"]   
    covlist = []    
    num = 0
    for name in templates:
        num += 1   
        X,Y = in_put2("%s%s.txt"%(path,name))
        X = [float(i) for i in X]   
        Y = [float(i) for i in Y]
        Ynew = []    
        ymax = max(Y)
        for y in Y:
            Ynew.append(1.*y/ymax)  
        plot.plot(X, Ynew, mew = 0, ms = 0, marker = "d", lw = 2, label = name) 
        covlist += Ynew
        #plot.hist(covlist)
        plot.legend()
        plot.yticks((0,0.2,0.4,0.6,0.8,1.0,1.1))
        plot.xlim(0,5500)
        if num%2 == 0:
            plot.savefig("%scoverages%s.svg"%(path, name), format = 'svg')
            plot.close()
    return

def coverage_subplots(path):
            
    X,Y = in_put2("%scovS4-5.txt"%(path))
    X = [float(i) for i in X]   
    Y = [float(i) for i in Y]
    X, Y = zip(*sorted(zip(X, Y)))
    Ynew = []    
    ymax = max(Y)
    for y in Y:
        Ynew.append(1.*y/ymax)  
    plot.plot(X, Ynew, mew = 0, ms = 0, marker = "d", lw = 1, color = "#82a53c" ) 
    #plot.hist(covlist)
    #plot.xticks((-50,0,200,400,600,800,1000,1200,1400,1600, 1700))
    plot.yticks((0,0.2,0.4,0.6,0.8,1.0,1.1))
    plot.savefig("%scoveragesubplots4-5.svg"%path, format = 'svg')
    return
    
def readcountslog(path):
    X,Y = in_put2("%sreadsbar100.txt"%path)
    Ynew = []
    count = 0
    for point in Y:
        count = count + 1
        Ynew.append(int(point))
    
    print "Number of reads processed: %d" %count
    #plot.figure(figsize=(20,12))
    bins=[1,2,3,4,5,6,7,8,9]+list(10 ** np.linspace(np.log10(10), np.log10(1e5), 50))
    print bins
    n, bins, patches = plot.hist(Ynew, bins, color = ("#82a53c"), histtype = "stepfilled", normed = True)
    plot.xlim(1e2,1e5)
    plot.xscale("log")
    plot.ylim(0,1e-3)    
    #plot.yticks((0.2e-3,0.4e-3,0.6e-3,0.8e-3,1e-3,1.2e-3,1.4e-3,1.6e-3,1.8e-3,2e-3,2.2e-3))
    #plot.xlim(1,1000)
    #plot.ylim(0, 120000)
    #plot.yticks((10e3, 20e3,30e3,40e3,50e3, 60e3))
    #plot.xticks((1, 50, 200, 400, 600, 800, 1000))    
    '''#plot the cumulative function
    y = [0]
    i = 0
    for num in n:
        y.append(y[i]+num)
        i+=1
    plot.close()
    plot.plot(bins, y, lw = 2, color = "#82a53c")
    #plot.xlim(0, 1000)
    #plot.ylim(0, 1.4e6)
    '''
    plot.savefig("%sreadcountsnormedlogbinsinset.svg"%path, format = 'svg')
    return    
    
def readcounts(path):
    X,Y = in_put2("%sreadsbartotal.txt"%path)
    Ynew = []
    count = 0
    for point in Y:
        count = count + 1
        Ynew.append(int(point))
    print "Number of reads processed: %d" %count
    #plot.figure(figsize=(20,12))
    n, bins, patches = plot.hist(Ynew, 1e5, color = ("#82a53c"), histtype = "stepfilled", normed = False)
    plot.xlim(1e0,1e5)
    #plot.ylim(0,2e-3)    
    #plot.yticks((0.2e-3,0.4e-3,0.6e-3,0.8e-3,1e-3,1.2e-3,1.4e-3,1.6e-3,1.8e-3,2e-3,2.2e-3))
    #plot.xlim(1,1000)
    #plot.ylim(0, 120000)
    #plot.yticks((10e3, 20e3,30e3,40e3,50e3, 60e3))
    #plot.xticks((1, 50, 200, 400, 600, 800, 1000))    
    #plot the cumulative function
    '''
    y = [0]
    i = 0
    for num in n:
        y.append(y[i]+num)
        i+=1
    plot.close()
    plot.plot(bins, y, lw = 2, color = "#82a53c")
    plot.xscale("log")
    #plot.xlim(0, 1000)
    #plot.ylim(0, 1.4e6)
    
    plot.savefig("%sreadcounts.svg"%path, format = 'svg')
    '''
    return    

def GC_histo(path):
    data = in_put("%sGChisto.txt"%path)
    data = [p * 100 for p in data]
    n, bins, patches = plot.hist(data, 16, range = (0, 100), color = ("#82a53c"), histtype = "stepfilled", normed = True)
    data = in_put("%sGC_histotest.txt"%path)
    data = [p * 100 for p in data]    
    n, bins, patches = plot.hist(data, 16, range = (0, 100), color = ("blue"), alpha = 0.5, histtype = "stepfilled", normed = True)    
    plot.savefig("%sGChisto.svg"%path, format = 'svg')    
    return

def qscores(path):
    X, Y = in_put2("%scontigqscores.txt"%path)
    plot.plot(X,Y, mew = 0, ms = 0, color = '#82a53c', lw = 1 )
    #plot.ylim(20, 50)
    plot.savefig("%sqscores.svg"%path, format = 'svg')
    return

def coverage_identity(path): #percent match of long contigs
    X, Y = in_put2("%scontigsidentity.txt"%path)
    plot.scatter(X, Y, color = "#82a53c")
    #plot.xlim(0.60, 1.02)    
    #plot.ylim(-150, 4000)
    #plot.xlim(0.996, 1.0005)
    #plot.xticks((0.996, 0.997, 0.998, 0.999, 1))
    plot.yscale("log")
    plot.savefig("%spercentmatchinset.svg"%path, format = 'svg')        
    return   
    
def read_lengths(path): #histogram of read lengths of contigs
    data = in_put("%sallcontig-length.txt"%path)
    count = 0
    for point in data:
        count = count + 1
    print "Number of contigs processed: %d" %count
    plot.hist(data, bins=50, normed = False, color = ("#82a53c"), histtype = "stepfilled")  
    plot.savefig("%sreadlengths.svg"%path, format = 'svg')
    return
 
def entropies(path):
    W,X,Z,Y = in_put4("%sentropies.txt"%path)
    Ynew = []
    count = 0    
    for p in Y:
        count += 1
        Ynew.append(float(p))
    print "data points processed: %d" %count
    #plot.hist(Ynew, bins=50, color = ("#82a53c"), histtype = "stepfilled")   
    plot.plot(Z, Y, lw = 0, ms = 2, marker = "D", color = "#82a53c")    
    plot.show()
    plot.xscale("log")
    plot.xlabel("# of reads")
    plot.ylabel("Entropy")
    #plot.savefig("%sentropy.svg"%path, format = 'svg')
    return

def reads_success(path):
    X,Y = in_put2("%sreads_success.txt"%path)
    count = 0
    for p in X:
        count += 1
    print "data points processed: %d" %count
    plot.scatter(X, Y, color = "#82a53c")
    plot.xlabel("# of reads")
    plot.ylabel("Assembly Success Rate")
    return

def entropy_success(path):
    X,Y = in_put2("%sentropy_success.txt"%path)
    X, Y = zip(*sorted(zip(X, Y)))
    count = 0
    for p in X:
        count += 1
    print "data points processed: %d" %count
    plot.plot(X, Y, color = "#82a53c", marker = "D", lw = 1, ms = 3)
    plot.ylim(-0.2,1.2)
    plot.xlabel("entropy of bargroup")
    plot.ylabel("Assembly Success Rate")
    plot.savefig("%sentropysuccess17.svg"%path, format = 'svg')
    return

def hammings(path):
    data = in_put("%shamminglist100adopted.txt"%path) 
    count = 0
    num = 0
    for p in data:
        count += 1
        if p == 1:
           num += 1
    print "Number of data points processed: %d" %count
    print "Hammings 1 is %d out of %d total" %(num, count)
    n, bins, patches = plot.hist(data, bins=16, normed = True, range = (0,16), align = 'left', color = ("#82a53c"), histtype = "stepfilled")
    
    data = in_put("%shamminglisttest1624.txt"%path) 
    count = 0
    num = 0
    for p in data:
        count += 1
        if p == 1:
           num += 1
    print "Number of data points processed: %d" %count
    print "Hammings 1 is %d out of %d total" %(num, count)
    n, bins, patches = plot.hist(data, bins=16, normed = True, range = (0,16), align = 'left', color = ("blue"), histtype = "step", lw = 5, alpha = 0.5)
    
    plot.xlim(0, 8)
    #plot.xticks((1,2,3,4))
    #plot.ylim(0,1e-4)    
    '''
    y = [0]
    i = 0
    for num in n:
        y.append(y[i]+num)
        i+=1
    plot.close()
    plot.plot(bins, y, lw = 1, color = "#82a53c")
    plot.xlim(0, 4)
    plot.xticks(range(4))
    plot.ylim(0, 0.0003)

    '''
    plot.savefig("%shammingsadop.svg"%path, format = 'svg')
    return

def readcountsbybar(path):
    X, Y = in_put2("%sreadsbar500unique.txt"%path)
    Ynew = []
    count = 0
    for point in Y:
        count = count + 1
        Ynew.append(int(point))
    print "Number of reads processed: %d" %count 
    Ynew.sort()
    plot.plot(Ynew, color = "#82a53c")

    X, Y = in_put2("%sreadsbar500.txt"%path)
    Ynew = []
    count = 0
    for point in Y:
        count = count + 1
        Ynew.append(int(point))
    print "Number of reads processed: %d" %count 
    Ynew.sort()
    plot.plot(Ynew)
    

    plot.yscale("log")
 
    plot.savefig("%sreadcountsbybar.svg"%path, format = 'svg')
    return    

def readcountscumplot(path):
    X, Y = in_put2("%sreadsbartotal.txt"%path)
    Ynew = []
    count = 0
    for point in Y:
        count = count + 1
        Ynew.append(int(point))
    print "Number of reads processed: %d" %count
    #plot.figure(figsize=(20,12))
    #make it based on # of reads instead of # of bars
    yreads = []
    for p in Ynew:
        for i in range(p):
            yreads.append(p)
            
    n, bins, patches = plot.hist(yreads, 1e4, color = ("#82a53c"), histtype = "step", normed = True, cumulative = True)
    '''#plot the cumulative function for % of data
    y = [0]      #cumulative data for each bin
    i = 0
    for num in n:
        y.append(y[i]+num)
        i+=1
    plot.close()
    plot.plot(bins, y, lw = 2, color = "#82a53c")
    '''
    plot.xscale("log")
    #plot.ylim(0, 1.4e6)
    plot.axvline(x=500)
    plot.xlabel("Size of Bargroup")
    plot.ylabel("# of reads cumulative")

    plot.savefig("%sreadcountscum.svg"%path, format = 'svg')
    return    
    
    
def read_pie(path):
    X, Y = in_put2("%sreadsbar5k.txt"%path)
    Ynew = []
    count = 0
    for p in Y:
        Ynew.append(int(p))
        count += 1
    print "%d data points"%count
    Ynew.sort()
    plot.pie(Ynew, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plot.axis('equal')
    #plot.savefig("%spie5.svg"%path, format = "svg")
    return

def rarefaction(path):
    X, Y = in_put2("%srarefaction-100.txt"%path)
    plot.plot(X,Y, color = "#82a53c")
    X, Y = in_put2("%srarefaction-300.txt"%path)
    plot.plot(X,Y, color = "b")
    X, Y = in_put2("%srarefaction-500.txt"%path)
    plot.plot(X,Y, color = "g")
    X, Y = in_put2("%srarefaction-1000.txt"%path)
    plot.plot(X,Y, color = "cyan")
    plot.ylabel("number of barcode clusters over threshold")
    plot.xlabel("total reads")
    plot.savefig("%srarefaction.svg"%path, format = "svg")
    return

def GC_cov(path):
    for key in ["C"]:
        X, Y = in_put2("%sGC-cov%s.txt"%(path, key))
        plot.scatter(X, Y)
    return

def GC_pos(path):
    n = 0
    for key in ["C", "D", "E", "G", "H", "I", "J", "L"]:
        X, Y = in_put2("%sGC_pos%s.txt"%(path, key))
        plot.plot(X, Y, label = key)
        plot.legend()
        plot.xlim(0, 5000)
        plot.ylim(0.25, 0.75)
        n += 1
        if n >= 2:
            plot.savefig("%s%sGC_pos.svg"%(path,key), format = "svg")
            plot.close()
            n = 0
    return

def coverages(path):
    templates = ["C", "D", "E", "G", "H", "I", "J", "L"]   
    num = 0
    for name in templates:
        num += 1   
        X,Y = in_put2("%s%scov.txt"%(path,name))
        X = [float(i) for i in X]   
        Y = [float(i) for i in Y]
        Ynew = []    
        ymax = max(Y)
        for y in Y:
            Ynew.append(1.*y/ymax)
        fig, ax1 = plot.subplots()
        ax2 = ax1.twinx()
        ax1.plot(X, Ynew, mew = 0, ms = 0, lw = 1, label = name) 
        plot.yticks((0,0.2,0.4,0.6,0.8,1.0,1.1))
        ax2.plot()
        X, Y = in_put2("%sGC_pos%s.txt"%(path, name))
        ax2.plot(X, Y, label = name)
        plot.ylim(0.25, 0.75)
        plot.legend()
        plot.xlim(0, 5500)

        #if num%2 == 0:
        plot.savefig("%scoverages%s.svg"%(path, name), format = 'svg')
        plot.close()
    return

def contam(path):
    X,Y = in_put2("%scontam300.txt"%path)
    Ynew = []
    count = 0    
    for p in Y:
        count += 1
        Ynew.append(float(p))
    print "data points processed: %d" %count
    Ynew.sort()
    X = [1.*x/len(Y)*100 for x in range(1, len(Y)+1)]
    plot.bar(X, Ynew, color = "#82a53c")    
    plot.ylim(0,10)
    plot.xlim(0,100)
    plot.show()
    plot.xlabel("Fraction of Bargroups (%)")
    plot.ylabel("Contaminating Reads (%)")
    plot.savefig("%scontam.svg"%(path), format = 'svg')
    return
    
if __name__ == "__main__":
    path = "D:\\NGS\\lib17\\"
    #path = "C:\\NGS working directory\\lib16C\\"
    #templatecounts(path)
    #purity(path)
    #purity_by_bar(path)
    #coverages(path)
    #coverage_subplots(path)
    #readcountslog(path)
    #readcounts(path)
    #readcountsbybar(path)
    #read_pie(path)
    #GC_histo(path)
    #qscores(path)
    #coverage_identity(path)
    #read_lengths(path)
    #entropies(path)
    #reads_success(path) 
    #entropy_success(path)
    #hammings(path)
    rarefaction(path)
    #readcountscumplot(path)
    #GC_cov(path)
    #GC_pos(path)
    #contam(path)
    print "done"

    
