#!/usr/bin/env python
__version__='0.6.0'
last_update='2021-10-19'
author='Damien Marsic, damien.marsic@aliyun.com'

import argparse,sys,glob,gzip,os,time
import numpy as np
from Bio.Seq import Seq
from collections import defaultdict
from matplotlib import pyplot as plt

def version():
    print('\n  Project: '+sys.argv[0][max(sys.argv[0].rfind('\\'),sys.argv[0].rfind('/'))+1:-3]+'\n  version: '+__version__+'\n  Latest update: '+last_update+'\n  Author: '+author+'\n  License: GNU General Public v3 (GPLv3)\n')

def override(func):
    class OverrideAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            func()
            parser.exit()
    return OverrideAction

def check_file(filename,strict):
    try:
        f=open(filename,'r')
    except IOError:
        if strict:
            print("\n  File "+filename+" could not be found!\n")
            sys.exit()
        else:
            return False
    else:
        return True

def readfasta(filename):
    f=open(filename,'r')
    seq=1
    for line in f:
        if line and line[0]=='>':
            if seq!=1:
                print('\n  more than one sequence is present in the genome file!\n')
                sys.exit()
            seq=''
            continue
        if seq!=1:
            l=line.lower().strip()
            l=''.join(x for x in l if not x.isdigit() and x!=' ')
            seq+=l
    f.close()
    if not seq or seq==1:
        print('\n  Genome sequence could not be found in '+filename+' !\n')
        sys.exit()
    return seq

def rename(name):
    if glob.glob(name):
        t=str(time.time())
        n=name[:name.rfind('.')]+'-'+t[:t.find('.')]+name[name.rfind('.'):]
        os.rename(name,n)
        print('\n  Existing '+name+' file was renamed as '+n+'\n  Creating new '+name+' file...\n')

def main():
    parser=argparse.ArgumentParser(description="Detection of large deletions in populations of circular genomes. For full documentation, visit: https://delfind.readthedocs.io")
    parser.add_argument('-v','--version',nargs=0,action=override(version),help="Display version")
    subparser=parser.add_subparsers(dest='command',required=True)
    parser_a=subparser.add_parser('map',help="Map read pairs to genome, arguments are optional if no ambiguity (if 2 arguments: read files, if 1 argument: genome sequence file)")
    parser_a.add_argument('R1',nargs='?',default='',type=str,help="File containing the R1 reads in fastq or fastq.gz format")
    parser_a.add_argument('R2',nargs='?',default='',type=str,help="File containing the R2 reads in fastq or fastq.gz format")
    parser_a.add_argument('genome',nargs='?',default='',type=str,help="File containing the genome in fasta format")
    parser_a.add_argument('-p','--probe',type=int,default=None,help="Size in nt of read fragment to be mapped (default: 75%% of read length)")
    parser_a.add_argument('-s','--slide',type=int,default=30,help="Range in nt of sliding the read fragment to generate a probe to be mapped (default: 30 nt)")
    parser_a.add_argument('-l','--limit',type=int,default=None,help="Stop after a number of read pairs have beem processed, instead of processing everything (default: no limit)")
    parser_b=subparser.add_parser('analyze',help="Analyze mapped read pairs to detect and quantify large deletions")
    parser_b.add_argument('map',nargs='?',default='',type=str,help="File containing the read pair maps in csv format (optional if no ambiguity)")
    parser_b.add_argument('genome',nargs='?',default='',type=str,help="File containing the genome in fasta format (optional if no ambiguity)")
    parser_b.add_argument('-c','--clean',type=int,default=100,help="Filter out deletions larger than value expressed in %% of genome size (default: 100%% = no filtering)")
    parser_b.add_argument('-t','--threshold',type=int,default=-1,help="Minimum distance between paired reads for the pair to be considered a deletion (default: -1 = autodetect)")
    parser_b.add_argument('-m','--merge',type=int,default=-1,help="Maximum distance for merging deletions (default: -1 = autodetect)")
    parser_b.add_argument('-e','--exclude',type=int,default=1,help="Deletions defined by a number of reads equal or lower than that number will be excluded from deletion frequency plot (default: 1)")
    parser_b.add_argument('-f','--format',type=str,default='png',help="File format for figures. Choices: svg, png, jpg, pdf, ps, eps, pgf, raw, rgba, tif. Default: png")

    args=parser.parse_args()
    if args.command=='map':
        rmap(args)
    if args.command=='analyze':
        analyze(args)

def rmap(args):
    global overlap,probe,slide,ref,cnt,Genome,rlength
    overlap=0 #   autodetect from a few hundred reads ***********************************************************
    rlength=150 # autodetect from first few reads !!  ***********************************************************
    R1=args.R1
    R2=args.R2
    Genome=args.genome
    limit=args.limit
    probe=args.probe
    slide=args.slide
    if args.R1 and not args.R2 and not args.genome:
        R1=''
        Genome=args.R1
    if not R1:
        x=glob.glob('*1.f*q*')
        if len(x)==1:
            R1=x[0]
        x=glob.glob('*2.f*q*')
        if len(x)==1:
            R2=x[0]
    if not R1 or not R2:
        print('\n  Reads file could not be detected unambiguously!\n')
        sys.exit()
    if not Genome:
        x=glob.glob('*.f*a')
        if len(x)==1:
            Genome=x[0]
        else:
            print('\n  Genome file could not be detected unambiguously!\n')
            sys.exit()
    if not probe:
        probe=round(rlength*0.75)
    if not slide:
        slide=1
        print('\n  Slide value automatically changed to 1 (can not be 0)')
    if slide<0:
        print('\n  Slide value can not be negative!\n')
        sys.exit()
    print('\n  Starting delfind map with the following arguments:')
    print('  R1: '+R1+'\n  R2: '+R2+'\n  Genome: '+Genome+'\n  Probe size: '+str(probe)+'\n  Slide range: '+str(slide)+'\n  Limit: '+str(limit))
    check_file(R1,True)
    check_file(R2,True)
    check_file(Genome,True)
    gfile=Genome[:Genome.rfind('.')]
    Genome=readfasta(Genome)
    ref=Genome+Genome[:rlength]
    if R1[-3:]=='.gz':
        f1=gzip.open(R1,'rt')
    else:
        f1=open(R1,'r')
    if R2[-3:]=='.gz':
        f2=gzip.open(R2,'rt')
    else:
        f2=open(R2,'r')
    z=2
    cnt=[0,0,0,0]
    mapped=defaultdict(int)
    print('\n  Read pairs processed   Read pairs mapped')
    while True:
        if limit and cnt[0]>=limit:
            break
        while z<4:
            l1=f1.readline().strip()
            l2=f2.readline().strip()
            if not l1 or not l2:
                break
            z+=1
            if z==3:
                l1=l1[:l1.find(' ')]
                l2=l2[:l2.find(' ')]
                if l1!=l2 or l1[0]!='@':
                    print('\n  Missing reads or wrong reads file format! Please use raw read files!\n')
                    sys.exit()
        if not l1 or not l2:
            break
        z=0
        l1=l1.lower()
        l2=l2.lower()
        b1=-1
        b2=-1
        cnt[0]+=1
        a1=locate(l1,0)
        a2=locate(l2,0)
        if a2!=-1:
            b1=locate(str(Seq(l1).reverse_complement()),1)
        if a1!=-1:
            b2=locate(str(Seq(l2).reverse_complement()),1)
        if a1!=-1 and b1!=-1:
            cnt[2]+=1
        if a2!=-1 and b2!=-1:
            cnt[2]+=1
        temp=[]
        if a1!=-1 and b2!=-1:
            temp.append([a1,b2])
        if a2!=-1 and b1!=-1:
            temp.append([a2,b1])
        if not temp:
            continue
        for i in range(len(temp)):
            if temp[i][0]>=len(Genome):
                temp[i][0]-=len(Genome)
            if temp[i][0]>temp[i][1]+rlength+1:
                temp[i][1]+=len(Genome)
            temp[i][1]-=temp[i][0]
        if len(temp)>1:
            cnt[3]+=1
            if temp[0][1]<temp[1][1]:
                del temp[1]
            elif temp[0][1]>temp[1][1]:
                del temp[0]
        if len(temp)>1:
            continue
        cnt[1]+=1
        mapped[(temp[0][0],temp[0][1])]+=1
        x=str(cnt[0])
        y=str(cnt[1])
        print('\r'+' '*(22-len(x))+x+' '*(20-len(y))+y,end='')
    f1.close()
    f2.close()
    print('\n\n  Read pairs processed: '+str(cnt[0]))
    print('  Read pairs mapped: '+str(cnt[1]))
    print('  Reads with more than one match: '+str(cnt[2]))
    print('  Read pairs with more than one match: '+str(cnt[3])+'\n')
    if limit:
        gfile+='_'+str(limit)
    gfile+='_map.csv'
    rename(gfile)
    g=open(gfile,'w')
    for n in sorted(mapped):
        g.write(str(n[0])+','+str(n[1])+','+str(mapped[n])+'\n')
    g.write('\n')
    g.close()
    print('  Map was saved into file: '+gfile+'\n')

def locate(l,s):
    x=0
    while x<overlap+slide:
        if s:
            y=l[x:x+probe]
        else:
            y=l[-x-probe:len(l)-x]
        n=ref.find(y)
        if n!=-1 and (ref.count(y)==1 or (y in Genome[:rlength] and Genome.count(y)==1)):
            if not s:
                n+=probe
            break
        if n==-1:
            x+=1
            continue
        cnt[2]+=1
        w=0
        while x+probe+w+1<=len(l):
            w+=1
            if s:
                y=l[x:x+probe+w]
            else:
                y=l[-x-probe-w:len(l)-x]
            n=ref.count(y)
            if n==1:
                n=ref.find(y)
                if not s:
                    n+=probe+w
                break
            if not n:
                n=-1
                break
        break
    return n

def analyze(args):
    f=args.map
    Genome=args.genome
    clean=args.clean
    format=args.format
    threshold=args.threshold
    merge=args.merge
    exclude=args.exclude
    if f[-6:]=='.fasta' or f[-3:]=='.fa':
        Genome=args.map
        f=''
    if args.genome[-4:]=='.csv':
        f=args.genome
    if not f:
        x=glob.glob('*_map.csv')
        if len(x)==1:
            f=x[0]
        else:
            print('\n  Map file could not be detected unambiguously!\n')
            sys.exit()
    if not Genome:
        x=glob.glob('*.f*a')
        if len(x)==1:
            Genome=x[0]
        else:
            print('\n  Genome file could not be detected unambiguously!\n')
            sys.exit()
    if clean>100 or clean<=0:
        print('\n  Filter value should be any integer higher than 0 and not higher than 100 (recommended 50)!\n')
        sys.exit()
    if format not in ('svg','png','jpg','jpeg','pdf','ps','eps','pgf','raw','rgba','tif','tiff'):
        print("\n  File format not recognized! Options are svg, png, jpg, pdf, ps, eps, pgf, raw, rgba, tif, tiff.\n")
        sys.exit()
    check_file(f,True)
    check_file(Genome,True)
    gsize=len(readfasta(Genome))
    filter=clean*gsize/100
    rmap=np.genfromtxt(f,delimiter=',').astype(int)
    if rmap[-1][0]>=gsize:
        print("\n  Wrong genome file!\n")
        sys.exit()
    rmap=rmap[~(rmap[:,1]>filter)]
    print('\n  Starting delfind analyze with the following arguments:')
    if threshold<0:
        threshold=round(np.mean(rmap[rmap[:,1]>0][:,1])*2+500)
    if merge<0:
        merge=round(np.std(rmap[rmap[:,1]<threshold][:,1]))
    print('  Map file: '+f+'\n  Genome file: '+Genome+'\n  Clean: '+str(clean)+'%\n  Threshold: '+str(threshold)+' bp'+'\n  Merge: '+str(merge)+' bp'+'\n  Format: '+format+'\n')
    nr=np.sum(rmap[:,2])
    nd=rmap[rmap[:,1]>=threshold][:,0].size
    print('  Total number of read pairs: '+str(nr)+'\n  Deletion read pairs: '+str(nd)+' ('+str(round(nd/nr*100,2))+'%)\n')
    fig=plt.figure(figsize=(12,6.75))
    plt.title('Distribution of read pair distances along genome',size=15,weight='roman')
    plt.xlabel("Genome position")
    plt.ylabel("Distance (bp)")
    plt.xlim(0,gsize)
    plt.plot(rmap[rmap[:,1]<threshold][:,0],rmap[rmap[:,1]<threshold][:,1],label='<'+str(threshold),marker='.',linestyle='',markersize=1)
    plt.plot(rmap[rmap[:,1]>=threshold][:,0],rmap[rmap[:,1]>=threshold][:,1],label='>='+str(threshold),marker='.',linestyle='',markersize=1)
    plt.legend(loc='upper left')
    z='_c'+str(clean)
    if clean==100:
        z=''
    g=f[:f.rfind('_')]+z+'_distances.'+format
    plt.savefig(g,dpi=600)
    print('  Figure was saved into file: '+g+'\n')
    plt.close()
    x=range(gsize)
    y1=[0]*(gsize)
    y2=[0]*(gsize)
    for n in rmap:
        for i in range(n[1]):
            if n[0]+i<gsize and n[1]<threshold:
                y1[n[0]+i]+=n[2]
            elif n[0]+i<gsize and n[1]>=threshold:
                y2[n[0]+i]+=n[2]
            elif n[0]+i>=gsize and n[1]<threshold:
                y1[n[0]+i-gsize]+=n[2]
            elif n[0]+i>=gsize and n[1]>=threshold:
                y2[n[0]+i-gsize]+=n[2]
    fig=plt.figure(figsize=(12,6.75))
    plt.title('Sequencing coverage of read pair distances along genome',size=15,weight='roman')
    plt.xlabel("Genome position")
    plt.ylabel("Coverage")
    plt.xlim(0,gsize)
    plt.stackplot(x,y1,y2,labels=['<'+str(threshold),'>='+str(threshold)])
    plt.legend(loc='upper left')
    w='_c'+str(clean)
    if clean==100:
        w=''
    w+='_t'+str(threshold)
    g=f[:f.rfind('_')]+w+'_coverage.'+format
    plt.savefig(g,dpi=600)
    print('  Figure was saved into file: '+g+'\n')
    plt.close()
    dels=[]
    i=-1
    b=0
    x=[]
    rem=[]
    while i+1<len(rmap):
        i+=1
        if rmap[i][1]<threshold or i in rem:
            continue
        x=(rmap[i][0],)*2+(rmap[i][0]+rmap[i][1],)*2+(rmap[i][2],)
        rem.append(i)
        b=i
        while i+1<len(rmap):
            i+=1
            if rmap[i][1]<threshold or i in rem:
                continue
            if rmap[i][0]-x[1]>merge:
                break
            y=rmap[i][0]+rmap[i][1]
            if abs(y-x[2])>merge and abs(y-x[3])>merge:
                continue
            x=(min(x[0],rmap[i][0]),max(x[1],rmap[i][0]),min(x[2],y),max(x[3],y),x[4])
            rem.append(i)
        dels.append(x)
        i=b

# Merge deletions here !



    y=[k[4] for k in dels]
    i=y.index(max(y))
    x1=dels[i][0]
    x2=dels[i][1]
    y=np.sum(rmap[(rmap[:,0]>=x1)&(rmap[:,0]<=x2)][:,2])
    for i in range(len(dels)):
        dels[i]=dels[i]+(100*dels[i][4]/y,)
    dels.sort(key=lambda x:x[5], reverse=True)
    w+='_m'+str(merge)
    g=f[:f.rfind('_')]+w+'_frequencies.csv'
    np.savetxt(g,dels,delimiter=',',fmt=['%d','%d','%d','%d','%d','%f'])
    print('  Deletion frequencies were saved into file: '+g+'\n')



# Plot distribution of deletion _frequencies




    fig=plt.figure(figsize=(12,6.75))
    plt.title('Deletion frequency as a function of deletion size',size=15,weight='roman')
    plt.xlabel("Deletion size")
    plt.ylabel("Deletion frequency")
    plt.plot([n[2]-n[1] for n in dels],[n[5] for n in dels],marker='.',linestyle='',markersize=1)
    g=g.replace('frequencies.csv','frequency_size.'+format)
    plt.savefig(g,dpi=600)
    print('  Figure was saved into file: '+g+'\n')
    plt.close()

    fig=plt.figure(figsize=(12,6.75))
    plt.title('Deletion frequencies along genome',size=15,weight='roman')
    plt.xlabel("Genome position")
    plt.ylabel("Frequency (%)")
    plt.xlim(0,gsize)
    for n in dels:
        plt.plot((n[1],n[2]),(n[5],n[5]),alpha=0.4)
    g=g.replace('_frequency_size','_locations_frequencies')
    plt.savefig(g,dpi=600)
    print('  Figure was saved into file: '+g+'\n')
    plt.close()




# Show graph all / restrict to more than 1 read per deletions
# X axis: genome coordinate
# Y axis: frequency and deletion size
   # ex: bar chart with Y=frequency and color=deletion Size
   # ex: bubble chart with y=frequency and bubble size= deletion size
   # ex: y=deletion size, bubble size=frequency


if __name__ == '__main__':
    main()
