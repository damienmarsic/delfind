#!/usr/bin/env python
__version__='0.9.6'
last_update='2021-11-12'
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

def pr2(f,t):
    print(t)
    f.write(t+'\n')

def plt0(title,xlabel,ylabel):
    fig=plt.figure(figsize=(12,6.75))
    plt.title(title,size=15,weight='roman')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def plt1(g,h):
    plt.savefig(g,dpi=600)
    pr2(h,'  Figure was saved into file: '+g+'\n')
    plt.close()

def rename(name):
    if glob.glob(name):
        t=str(time.time())
        n=name[:name.rfind('.')]+'-'+t[:t.find('.')]+name[name.rfind('.'):]
        os.rename(name,n)
        print('\n  Existing '+name+' file was renamed as '+n+'\n  Creating new '+name+' file...\n')

def lncount(f):
    def _make_gen(reader):
        b=reader(1024*1024)
        while b:
            yield b
            b=reader(1024*1024)
    f_gen=_make_gen(f.read)
    return sum(buf.count(b'\n') for buf in f_gen)

def main():
    parser=argparse.ArgumentParser(description="Detection of large deletions in populations of circular genomes. For full documentation, visit: https://delfind.readthedocs.io")
    parser.add_argument('-v','--version',nargs=0,action=override(version),help="Display version")
    subparser=parser.add_subparsers(dest='command',required=True)
    parser_a=subparser.add_parser('map',help="Map reads to genome, arguments are optional if no ambiguity (if 2 arguments: read files, if 1 argument: genome sequence file)")
    parser_a.add_argument('R1',nargs='?',default='',type=str,help="File containing the R1 reads in fastq or fastq.gz format")
    parser_a.add_argument('R2',nargs='?',default='',type=str,help="File containing the R2 reads in fastq or fastq.gz format")
    parser_a.add_argument('genome',nargs='?',default='',type=str,help="File containing the genome in fasta format")
    parser_a.add_argument('-m','--minimum',type=int,default=None,help="Minimum size in nt of contiguous read sequence identity for mapping to be accepted (default: 75%% of read length)")
    parser_a.add_argument('-l','--limit',type=int,default=None,help="Stop after a number of read pairs have beem processed, instead of processing everything (default: no limit)")
    parser_b=subparser.add_parser('analyze',help="Analyze mapped read pairs to detect and quantify large deletions")
    parser_b.add_argument('map',nargs='?',default='',type=str,help="Prefix of files containing the read and read pair maps (optional if no ambiguity)")
    parser_b.add_argument('genome',nargs='?',default='',type=str,help="File containing the genome in fasta format (optional if no ambiguity)")
    parser_b.add_argument('-t','--threshold',type=int,default=-1,help="Minimum distance between paired reads for the pair to be considered a deletion (default: -1 = autodetect)")
    parser_b.add_argument('-m','--merge',type=int,default=-1,help="Maximum distance for merging deletions (default: -1 = autodetect)")
    parser_b.add_argument('-f','--format',type=str,default='png',help="File format for figures. Choices: svg, png, jpg, pdf, ps, eps, pgf, raw, rgba, tif. Default: png")
    parser_b.add_argument('-i','--include',type=float,default=0.2,help="Minimum deletion frequency in %% to be included in the report (default: 0.2)")
    parser_c=subparser.add_parser('snapgene',help="Add deletion features to genome file in Genbank format to be used in Snapgene viewer")
    parser_c.add_argument('deletions',nargs='?',default='',type=str,help="File containing the deletion frequencies in csv format (optional if no ambiguity)")
    parser_c.add_argument('genome',nargs='?',default='',type=str,help="File containing the annotated genome in Genbank format (optional if no ambiguity)")
    parser_c.add_argument('-i','--include',type=float,default=0.2,help="Minimum deletion frequency in %% to be included as a feature in the annotated genome file (default: 0.2)")
    parser_c.add_argument('-o','--outfile',type=str,default=None,help="Name of output file. If none provided, the input file name will be used for the output file, and the input file will be renamed.")
    args=parser.parse_args()
    if args.command=='map':
        rmap(args)
    if args.command=='analyze':
        analyze(args)
    if args.command=='snapgene':
        snapgene(args)

def rmap(args):
    global minimum,ref,Genome,rlength
    print('\n  Checking arguments...',end='')
    R1=args.R1
    R2=args.R2
    Genome=args.genome
    limit=args.limit
    minimum=args.minimum
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
    check_file(R1,True)
    check_file(R2,True)
    check_file(Genome,True)
    print('         OK\n\n  Counting reads...',end='')
    if R1[-3:]=='.gz':
        f1=gzip.open(R1,'r')
    else:
        f1=open(R1,'rb')
    if R2[-3:]=='.gz':
        f2=gzip.open(R2,'r')
    else:
        f2=open(R2,'rb')
    nr1=lncount(f1)//4
    nr2=lncount(f2)//4
    f1.close()
    f2.close()
    if nr1!=nr2:
        print('\n  Read files must have the same number of reads! Please use raw read files!\n')
        sys.exit()
    if R1[-3:]=='.gz':
        f1=gzip.open(R1,'rt')
    else:
        f1=open(R1,'r')
    if R2[-3:]=='.gz':
        f2=gzip.open(R2,'rt')
    else:
        f2=open(R2,'r')
    Z=2
    print('             OK\n\n  Determining Read length...',end='')
    cnt=0
    x=defaultdict(int)
    while cnt<100:
        while Z<4:
            l1=f1.readline().strip()
            l2=f2.readline().strip()
            if not l1 or not l2:
                break
            Z+=1
        if not l1 or not l2:
            break
        Z=0
        x[len(l1)]+=1
        x[len(l2)]+=1
        cnt+=1
    rlength=max(x.keys())
    if not minimum:
        minimum=round(rlength*0.75)
    elif minimum>rlength:
        minimum=rlength
    f1.seek(0)
    f2.seek(0)
    print('    OK\n\n  Starting delfind map with the following settings:')
    print('  R1: '+R1+' ('+str(nr1)+' reads)\n  R2: '+R2+' ('+str(nr2)+' reads)\n  Genome: '+Genome+'\n  Read length: '+str(rlength)+'\n  Minimum map size: '+str(minimum)+'\n  Limit: '+str(limit))
    gfile=Genome[:Genome.rfind('.')]
    Genome=readfasta(Genome)
    ref=Genome+Genome[:rlength]
    Z=2
    cnt=[0,0,0,0,0]
    rmap=defaultdict(int)
    pmap=defaultdict(int)
    y=np.arange(0,nr1,nr1/1000)
    x=[round(n) for n in y]
    z=np.arange(0,100,0.1)
    y=[str(round(n,1)) for n in z]
    show=dict(zip(x,y))
    print('\n  Processing reads...       0.0%',end='')
    while True:
        if limit and cnt[0]>=limit:
            break
        while Z<4:
            l1=f1.readline().strip()
            l2=f2.readline().strip()
            if not l1 or not l2:
                break
            Z+=1
            if Z==3:
                if l1[:l1.find(' ')]!=l2[:l2.find(' ')] or l1[0]!='@':
                    print('\n  Missing reads or wrong reads file format! Please use raw read files!\n')
                    sys.exit()
        if not l1 or not l2:
            break
        Z=0
        l1=l1.lower()
        l2=l2.lower()
        cnt[0]+=1
        a1=locate(l1)
        a2=locate(l2)
        b1=locate(str(Seq(l1).reverse_complement()))
        b2=locate(str(Seq(l2).reverse_complement()))
        A=[a1,b1,a2,b2]
        x=sum([1 for n in A if n])
        if cnt[0] in show:
            x=show[cnt[0]]
            print('\r  Processing reads...      '+' '*(4-len(x))+x+'%',end='')
        if not x:
            continue
        if x==4:
            a=a1[0][1]-a1[0][0]+b2[0][1]-b2[0][0]
            b=a2[0][1]-a2[0][0]+b1[0][1]-b1[0][0]
        if (x==4 and a>b) or (x==3 and a1 and b2):
            a2=[]
            b1=[]
        elif (x==4 and b>a) or (x==3 and a2 and b1):
            a1=[]
            b2=[]
        elif x==2 and not (a1 and b2) and not (a2 and b1):
            z=max([n[0][1]-n[0][0] for n in A if n])
            for i in range(len(A)):
                if A[i] and A[i][0][1]-A[i][0][0]<z:
                    A[i]=[]
            a1,b1,a2,b2=A
        if (a1 and b2 and b1 and a2) or (a1 and b2 and len(a1+b2)>2) or (a2 and b1 and len(a2+b1)>2):
            cnt[4]+=1
        elif (a1 and b2) or (a2 and b1):
            cnt[2]+=1
        for n in ((a1,b1),(a2,b2)):
            if len(n[0]+n[1])>1:
                cnt[3]+=1
            elif len(n[0]+n[1])==1:
                cnt[1]+=1
        for n in ((a1,b2),(a2,b1)):
            if not (n[0] and n[1]):
                continue
            for x in n[0]:
                for y in n[1]:
                    z=[x[1],y[0]]
                    if z[0]>=len(Genome):
                        z[0]-=len(Genome)
                    if z[0]>z[1]+rlength+1:
                        z[1]+=len(Genome)
                    z[1]-=z[0]
                    pmap[tuple(z)]+=1
        for n in (a1,b2,a2,b1):
            if not n:
                continue
            for x in n:
                rmap[(x[0],x[1]-x[0])]+=1
    f1.close()
    f2.close()
    print('\b\b\b\b\b\b100.0%\n\n  Read pairs processed: '+str(cnt[0]))
    print('  Reads mapped exactly once: '+str(cnt[1]))
    print('  Read pairs mapped exactly once: '+str(cnt[2]))
    print('  Reads mapped more than once: '+str(cnt[3]))
    print('  Read pairs mapped more than once: '+str(cnt[4])+'\n')
    print('  Total mapped reads: '+str(cnt[1]+cnt[3])+' ('+str(round((cnt[1]+cnt[3])/cnt[0]*50,2))+'%)')
    print('  Total mapped read pairs: '+str(cnt[2]+cnt[4])+' ('+str(round((cnt[2]+cnt[4])/cnt[0]*100,2))+'%)')
    if limit:
        gfile+='_'+str(limit)
    g=open(gfile+'_rmap.csv','w')
    for n in sorted(rmap):
        g.write(str(n[0])+','+str(n[1])+','+str(rmap[n])+'\n')
    g.write('\n')
    g.close()
    print('\n  Read map was saved into file: '+gfile+'_rmap.csv')
    g=open(gfile+'_pmap.csv','w')
    for n in sorted(pmap):
        g.write(str(n[0])+','+str(n[1])+','+str(pmap[n])+'\n')
    g.write('\n')
    g.close()
    print('  Pair map was saved into file: '+gfile+'_pmap.csv\n')

def locate(l):
    L=[]
    a=rlength-minimum
    b=len(l)-a
    if not l[a:b] in ref:
        return L
    x=0
    y=len(l)
    while x<a:
        if l[x:y] in ref:
            break
        x+=10
        y-=10
    else:
        x=a
        y=b
    while x>0:
        x-=1
        if l[x:y] not in ref:
            x+=1
            break
    while y<len(l):
        y+=1
        if l[x:y] not in ref:
            y-=1
            break
    if y-x<minimum:
        return L
    a=0
    z=l[x:y]
    while True:
        n=ref.find(z,a)
        if n==-1:
            break
        a=n+1
        if z not in Genome[:rlength]:
            L.append((n,n+len(z)))
    return L

def analyze(args):
    fname=args.map
    Genome=args.genome
    format=args.format
    threshold=args.threshold
    merge=args.merge
    include=args.include
    if fname[-6:]=='.fasta' or fname[-3:]=='.fa':
        Genome=args.map
        fname=args.genome
    x=glob.glob(fname+'*_rmap.csv')
    y=glob.glob(fname+'*_pmap.csv')
    if len(x)==1 and len(y)==1 and x[0][:x[0].rfind('_')]==y[0][:y[0].rfind('_')]:
        fname=x[0][:x[0].rfind('_')]
    else:
        print('\n  Map files could not be detected unambiguously!\n')
        sys.exit()
    if not Genome:
        x=glob.glob('*.f*a')
        if len(x)==1:
            Genome=x[0]
        else:
            print('\n  Genome file could not be detected unambiguously!\n')
            sys.exit()
    if format not in ('svg','png','jpg','jpeg','pdf','ps','eps','pgf','raw','rgba','tif','tiff'):
        print("\n  File format not recognized! Options are svg, png, jpg, pdf, ps, eps, pgf, raw, rgba, tif, tiff.\n")
        sys.exit()
    if include<=0 or include>100:
        print('\n  Include value should be between 0 and 100!\n')
        sys.exit()
    check_file(fname+'_rmap.csv',True)
    check_file(fname+'_pmap.csv',True)
    check_file(Genome,True)
    gsize=len(readfasta(Genome))
    pmap0=np.genfromtxt(fname+'_pmap.csv',delimiter=',').astype(int)
    rmap0=np.genfromtxt(fname+'_rmap.csv',delimiter=',').astype(int)
    if pmap0[-1][0]>=gsize or rmap0[-1][0]>=gsize:
        print("\n  Wrong genome file!\n")
        sys.exit()
    x=pmap0[pmap0[:,2]>1].tolist()
    pmap1=np.delete(pmap0,2,1)
    for i in range(len(x)):
        for j in range(x[i][2]-2):
            x.append(x[i][:-1])
        del x[i][-1]
    pmap1=np.append(pmap1,x,0)
    del x
    pmap1=pmap1[pmap1[:,0].argsort()]
    ins=round(np.mean(pmap1[:,1]))




###  plan to calculate threshold from distance distribution:

    d=defaultdict(int)
    for n in pmap0:
        d[n[1]]+=n[2]


#### to be continued....

    if threshold<0:
        threshold=round(np.mean(pmap1[pmap1[:,1]>0][:,1])*2+500)  #  Temporary threshold definition !


    ins0=round(np.mean(pmap1[pmap1[:,1]<threshold][:,1]))
    ins1=round(np.mean(pmap1[pmap1[:,1]>=threshold][:,1]))
    sd=round(np.std(pmap1[pmap1[:,1]<threshold][:,1]))
    f=fname
    if threshold>=0:
        f=fname+'_t'+str(threshold)
    if merge>=0:
        f=fname+'_m'+str(merge)
    if merge<0:
        merge=sd
    h=open(f+'_report.txt','w')
    pr2(h,'\n  Starting delfind analyze with the following settings:')
    pr2(h,'  Read map file: '+fname+'_rmap.csv\n  Read pair map file: '+fname+'_pmap.csv\n  Genome file: '+Genome+'\n  Threshold: '+str(threshold)+' bp'+'\n  Merge: '+str(merge)+' bp'+'\n  Format: '+format+'\n')
    pr2(h,'  Number of reads in read file: '+str(np.sum(rmap0[:,2]))+'\n')
    nrp=np.size(pmap1,0)
    nd=np.size(pmap1[pmap1[:,1]>=threshold],0)
    pr2(h,'  Total number of read pairs in read pair file: '+str(nrp)+'\n  Deletion read pairs: '+str(nd)+' ('+str(round(nd/nrp*100,2))+'%)')
    pr2(h,'  Average insert size: '+str(ins)+' bp\n  Average under threshold insert size: '+str(ins0)+' bp\n  Average above threshold insert (deletion) size: '+str(ins1)+' bp\n')
    plt0('Distribution of read mapped sizes','Size (nt)','Number of reads')
    x=defaultdict(int)
    for n in rmap0:
        x[n[1]]+=n[2]
    n=np.array([(i,x[i]) for i in x])
    plt.plot(n[:,0],n[:,1],marker='.',linestyle='',markersize=1)
    g=f+'_mapped-size-distr.'+format
    plt1(g,h)
    plt0('Distribution of read pair distances','Distance (bp)','Number of read pairs')
    x=np.array([(n,d[n]) for n in d if n<threshold])
    plt.plot(x[:,0],x[:,1],label='<'+str(threshold),marker='.',linestyle='',markersize=1)
    x=np.array([(n,d[n]) for n in d if n>=threshold])
    plt.plot(x[:,0],x[:,1],label='>='+str(threshold),marker='.',linestyle='',markersize=1)
    plt.legend(loc='upper right')
    g=f+'_distances-distr.'+format
    plt1(g,h)
    plt0('Distribution of mapped reads along genome','Genome position','Read mapped size (nt)')
    plt.xlim(0,gsize)
    plt.plot(rmap0[:,0],rmap0[:,1],marker='.',linestyle='',markersize=1)
    g=f+'_mapped-reads-genome.'+format
    plt1(g,h)
    plt0('Distribution of read pair distances along genome','Genome position','Distance (bp)')
    plt.xlim(0,gsize)
    plt.plot(pmap1[pmap1[:,1]<threshold][:,0],pmap1[pmap1[:,1]<threshold][:,1],label='<'+str(threshold),marker='.',linestyle='',markersize=1)
    plt.plot(pmap1[pmap1[:,1]>=threshold][:,0],pmap1[pmap1[:,1]>=threshold][:,1],label='>='+str(threshold),marker='.',linestyle='',markersize=1)
    plt.legend(loc='upper left')
    g=f+'_distances-genome.'+format
    plt1(g,h)
    plt0('Sequencing depth of mapped reads along genome','Genome position','Coverage')
    x=range(gsize)
    y=[0]*(gsize)
    for n in rmap0:
        for i in range(n[1]):
            if n[0]+i<gsize:
                y[n[0]+i]+=n[2]
            elif n[0]+i>=gsize:
                y[n[0]+i-gsize]+=n[2]
    plt.xlim(0,gsize)
    plt.stackplot(x,y)
    g=f+'_read-coverage.'+format
    plt1(g,h)
    n=[(x[i],y[i]) for i in range(len(x))]
    g=f+'_read-coverage.csv'
    np.savetxt(g,n,delimiter=',',fmt=['%d','%d'])
    pr2(h,'  Read coverage was saved into file: '+g+'\n')
    plt0('Sequencing depth of large (>='+str(threshold)+') read pair distances (deletions) along genome','Genome position','Coverage')
    x=range(gsize)
    y=[0]*(gsize)
    for n in pmap0[pmap0[:,1]>=threshold]:
        for i in range(n[1]):
            if n[0]+i<gsize:
                y[n[0]+i]+=n[2]
            elif n[0]+i>=gsize:
                y[n[0]+i-gsize]+=n[2]
    plt.xlim(0,gsize)
    plt.stackplot(x,y,color='red')
    g=f+'_distance-coverage.'+format
    plt1(g,h)


    dels=[]
    i=-1
    b=0
    rem=[]
    while i+1<len(pmap0):
        i+=1
        if i in rem or pmap0[i][1]<threshold:
            continue
        x=(pmap0[i][0],)*2+(pmap0[i][0]+pmap0[i][1],)*2+(pmap0[i][2],)
        rem.append(i)
        b=i
        while i+1<len(pmap0):
            i+=1
            if i in rem or pmap0[i][1]<threshold:
                continue
            if pmap0[i][0]-x[1]>merge:
                break
            y=pmap0[i][0]+pmap0[i][1]
            if abs(y-x[2])>merge and abs(y-x[3])>merge:
                continue
            x=(min(x[0],pmap0[i][0]),max(x[1],pmap0[i][0]),min(x[2],y),max(x[3],y),x[4]+pmap0[i][2])
            rem.append(i)
        dels.append(x)
        i=b
    rem=[]
    for i in range(len(dels)):
        if dels[i][4]==1:
            continue
        j=i+1
        while j<len(dels) and dels[j][0]<=dels[i][1]:
            if j not in rem and dels[j][3]>=dels[i][2] and dels[j][2]<=dels[i][3]:
                rem.append(j)
                dels[i]=(min(dels[i][0],dels[j][0]),max(dels[i][1],dels[j][1]),min(dels[i][2],dels[j][2]),max(dels[i][3],dels[j][3]),dels[i][4]+dels[j][4])
            j+=1
    rem.sort(reverse=True)
    for i in rem:
        del dels[i]
    y=[k[4] for k in dels]
    i=y.index(max(y))
    x1=dels[i][0]
    x2=dels[i][1]
    y=np.size(pmap1[(pmap1[:,0]>=x1)&(pmap1[:,0]<=x2)],0)
    for i in range(len(dels)):
        dels[i]=dels[i]+(100*dels[i][4]/y,)
    dels.sort(key=lambda x:x[5], reverse=True)
    g=f+'_del-frequencies.csv'
    np.savetxt(g,dels,delimiter=',',fmt=['%d','%d','%d','%d','%d','%f'])
    pr2(h,'  Deletion frequencies were saved into file: '+g+'\n')



# Add plot distributions of deletion frequencies / deletion sizes ?




    plt0('Deletion frequency as a function of deletion size','Deletion size','Deletion frequency')
    plt.plot([n[2]-n[1] for n in dels],[n[5] for n in dels],marker='.',linestyle='',markersize=1)
    g=f+'_frequency_size.'+format
    plt1(g,h)
    plt0('Deletion frequencies along genome','Genome position','Frequency (%)')
    plt.xlim(0,gsize)
    for n in dels:
        plt.plot((n[1],n[2]),(n[5],n[5]),alpha=0.4)
    g=f+'_locations_frequencies.'+format
    plt1(g,h)
    pr2(h,'  Top deletions (>='+str(include)+'%):\n        Genome location       Size (bp)   Frequency (%)')
    for n in dels:
        if n[5]<include:
            break
        q=n[2]
        if q>gsize:
            q-=gsize
        x='['+str(n[1]+1)+'..'+str(q)+']'
        y=str(n[2]-n[1])
        z=str(round(n[5],2))
        pr2(h,' '*(23-len(x))+x+' '*(16-len(y))+y+' '*(16-len(z))+z)
    pr2(h,'')
    h.close()
    print('  Report was saved into file: '+f+'_report.txt\n')

def snapgene(args):
    dels=args.deletions
    genome=args.genome
    include=args.include
    outfile=args.outfile
    if (genome[-4:]=='.csv' and (not dels or dels[-4:]!='.csv')) or ((dels[-3:]=='.gb' or dels[-4:]=='.gbk') and (not genome or (genome[-3:]!='.gb' and genome[-4:]!='.gbk'))):
        genome,dels=dels,genome
    if dels and dels[-4:]!='.csv':
        print('\n  Deletion frequencies file must be in csv format!\n')
        sys.exit()
    if genome and genome[-3:]!='.gb' and genome[-4:]!='.gbk':
        print('\n  Annotated genome file must be in Genbank format!\n')
        sys.exit()
    if not dels:
        x=glob.glob('*_frequencies.csv')
        if len(x)==1:
            dels=x[0]
        else:
            print('\n  Deletion frequencies file could not be detected unambiguously!\n')
            sys.exit()
    if not genome:
        x=glob.glob('*.gb*')
        if len(x)==1:
            genome=x[0]
        else:
            print('\n  Genome file could not be detected unambiguously!\n')
            sys.exit()
    check_file(dels,True)
    check_file(genome,True)
    if include<=0 or include>100:
        print('\n  Include value should be between 0 and 100!\n')
        sys.exit()
    if not outfile:
        outfile=genome
    if outfile[-3:]!='.gb' and outfile[-4:]!='.gbk':
        outfile+='.gb'
    D=[]
    f=open(dels,'r')
    for line in f:
        ln=line.strip().split(',')
        if float(ln[5])<include:
            break
        D.append((int(ln[1])+1,ln[2],str(round(float(ln[5]),2))))
    D.sort()
    f.close()
    f=open(genome,'r')
    x=f.read()
    a=x.find('FEATURES')
    b=x[a:].find('\n')+1
    y=x[:a+b]
    w=x[a+b:].split('\n')
    end=False
    for ln in w:
        if ln[:6]=='ORIGIN':
            end=True
        if not end and ln[:6]!='      ':
            a=int(''.join(filter(lambda i: i.isdigit(),ln[:ln.find('..')])))
            if D and D[0][0]<a:
                b=D.pop(0)
                a='     misc_feature    '+str(b[0])+'..'+b[1]+'\n                     /label=Del '+b[2]+'%\n                     /note="color: #ffff00; direction: BOTH"\n'
                if not a in x:
                    y+=a
        y+=ln+'\n'
    f.close()
    rename(outfile)
    with open(outfile,'w') as f:
        f.write(y)
    print('  Modified annotated genome was saved into file: '+outfile+'\n')


#  ADD COVERAGE !!!


if __name__ == '__main__':
    main()
