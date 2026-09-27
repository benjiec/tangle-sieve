#!/usr/bin/env python3

import sys
import numpy as np
import argparse
import re
from Bio.PDB import PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser



def pDockQ(structure,cutoff,disocut1=50,disocut2=70,disocut3=90):
    i=0
    tiny=1.e-20
    chains=[]
    NumRes=[]
    IF_plDDT=[]
    plDDT=[]
    NumDiso=[]
    NumOverlap=[]
    pdockq=[]
    length=[]
    maxlen=10000
    for chain in structure.get_chains():
        chains+=[chain]
        i+=1
        #print (chains)
    
    for c in range(len(chains)):
        interface_residues1=[]
        interface_residues2=[]
        k=tiny
        b=0
        for d in range(len(chains)):
            if c==d: continue
            for res1 in chains[c]:
                for res2 in chains[d]:
                    test=False
                    for i in res1:
                        if test:break
                        for j in res2:
                            dist = np.linalg.norm(i.coord - j.coord)
                            if dist < cutoff:
                                interface_residues1.append(res1.id[1])
                                interface_residues2.append(res2.id[1]+maxlen*d)
                                test=True
                                break
                            elif dist > 2*cutoff: # To speed up things
                                test=True
                                break
            if2=np.unique(interface_residues2)
            for res in chains[d]:
                if res.id[1] in if2:
                    b+=res['CA'].get_bfactor()
                    k+=1
        if1=np.unique(interface_residues1)
        b1=0
        i1=0
        NumDiso1=[0,0,0,0]
        for res in chains[c]:
            b1+=res['CA'].get_bfactor()
            i1+=1
            if res['CA'].get_bfactor()>disocut3: # >90
                NumDiso1[0]+=1
            elif res['CA'].get_bfactor()>disocut2: # 70-90
                NumDiso1[1]+=1
            elif res['CA'].get_bfactor()>disocut1: # 50-70
                NumDiso1[2]+=1
            else: # <50
                NumDiso1[3]+=1
            if res.id[1] in if1:
                b+=res['CA'].get_bfactor()
                k+=1
        length+=[i1]
        NumRes+=[if1.shape[0]+if2.shape[0]]
        IF_plDDT+=[b/k]
        plDDT+=[b1/i1]
        NumDiso+=[NumDiso1]
        pdockq+=[sigmoid(np.log(NumRes[c]+tiny)*IF_plDDT[c],*popt)]
        
    return (NumRes,IF_plDDT,plDDT,NumDiso,pdockq,length)


def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return (y)

#popt=[7.07140240e-01, 3.88062162e+02, 3.14767156e-02, 3.13182907e-02]

arg_parser = argparse.ArgumentParser(description="Calculates pDockQ from NumRes and IF_plDDT")
group = arg_parser.add_mutually_exclusive_group(required=True)
group.add_argument("-p","--pdb",type=argparse.FileType('r'),help="Input pdb file pLddt values in bfactor columns")
group.add_argument("-c","--cif",type=argparse.FileType('r'),help="Input cif file plddt values in bfactor columns")
arg_parser.add_argument("-v","--verbose", action='store_true', required=False,help="Verbose output")
arg_parser.add_argument("-C","--cutoff", type=float,required=False,default=10.0,help="Cutoff for defining distances")
args = arg_parser.parse_args()
cutoff=args.cutoff

args = arg_parser.parse_args()
#cutoff=10
cutoff2=3
Numoverlap=0
overlap=0
Name=""
if (args.cif):
    name=args.cif.name
    bio_parser = MMCIFParser()
    structure_file = args.cif
    structure_id = args.cif.name[:-4]
    structure = bio_parser.get_structure(structure_id, structure_file)
        
elif (args.pdb):
    Name=args.pdb.name
    bio_parser = PDBParser()
    structure_file = args.pdb
    structure_id = args.pdb.name[:-4]
    structure = bio_parser.get_structure(structure_id, structure_file)


if cutoff <=5:
    #5 SpearmanrResult(correlation=0.7647585237390458, pvalue=4.749030057232305e-280)
    popt=[6.96234405e-01, 2.35483775e+02, 2.25322970e-02, 2.88445245e-02]
    #0.7805034405869632
elif cutoff <=6:
    #6 SpearmanrResult(correlation=0.7708834427476546, pvalue=2.707297682746201e-287)
    popt=[7.02605033e-01, 2.91749822e+02, 2.70621128e-02, 2.25416051e-02]
    #0.7871982094514278
elif cutoff <=7:
    #7 SpearmanrResult(correlation=0.7709518988131879, pvalue=2.2402500804327052e-287)
    popt=[7.06385097e-01, 3.32456259e+02, 2.97005237e-02, 2.24488132e-02]
    #0.7859609807320201
elif cutoff <=8:
    #8 SpearmanrResult(correlation=0.7632969367380509, pvalue=2.3583905451705336e-278)
    popt=[7.18442739e-01,3.60791204e+02,3.01635944e-02, 2.04076969e-02]
    #0.7764648775754815
elif cutoff <=9:
    #9 SpearmanrResult(correlation=0.7496303495195178, pvalue=4.539049646719674e-263)
    popt=[7.23328534e-01, 3.80036094e+02, 3.06316084e-02, 1.98471192e-02]
    #0.7608417399783565
elif cutoff <=10:
    #10 SpearmanrResult(correlation=0.7330653937901442, pvalue=7.988440779428826e-246)
    popt=[7.20293782e-01, 3.95627723e+02, 3.15235037e-02, 2.37304238e-02]
    #0.7431426093979494
elif cutoff <=11:
    #11 SpearmanrResult(correlation=0.71288058226417, pvalue=1.7542846392453894e-226)
    popt=[7.22015998e-01, 4.09095024e+02, 3.11905555e-02, 2.59467513e-02]
    #0.7219615906164123
else:
    #12 SpearmanrResult(correlation=0.6938911161134763, pvalue=9.284495013784153e-210)
    popt=[7.20555781e-01, 4.21033584e+02, 3.09024241e-02, 2.88659629e-02]
    #0.7023000652310362
    
    
NumRes,IF_plDDT,plDDT,NumDiso,pdockq,length=pDockQ(structure,cutoff)
NumResOverlap,IF_plDDToverlap,plDDToverlap,NumDisooverlap,pdockqoverlap,lengthoverlap=pDockQ(structure,cutoff2)



Name=re.sub(r'.*/','',Name)
Name=re.sub(r'.pdb$','',Name)
Name=re.sub(r'.cif$','',Name)

#print (NumRes,tiny,IF_plDDT, popt)

#print (NumRes,IF_plDDT,plDDT,NumDiso,pdockq,NumResOverlap)

if (args.verbose):
    i=0
    print ("Name,id,Chain,pDockQ,NumRes,IF_plDDT,plDDTNumDiso1+90,NumDiso1-70-90,NumDiso1-50-70,NumDiso1-50,Length" )
    for chain in structure.get_chains():
        id=chain.get_id()
        print ("%s,%d,%s,%f,%d,%f,%f,%d,%d,%d,%d,%d,%d" %( Name,i,id,pdockq[i],NumRes[i],IF_plDDT[i],plDDT[i],NumDiso[i][0],NumDiso[i][1],NumDiso[i][2],NumDiso[i][3],NumResOverlap[i],length[i] ))
        i+=1
else:
    i=0
    for chain in structure.get_chains():
        id=chain.get_id()
        print (i,id,pdockq[i])
        i+=1

#    #print ("PdockQ: %.3f\nNumRes: %d\nIF_plDDT: %.2f\n" % (pDockQ,NumRes,IF_plDDT))
#    print ("Name,pDockQ,NumRes,IF_plDDT,plDDT1,plDDT2,NumDiso1+90,NumDiso1-70-90,NumDiso1-50-70,NumDiso1-50,NumDiso2+90,NumDiso2-70-90,NumDiso2-50-70,NumDiso2-50,NumOverlap,len1,len2" )
#    print ("%s,%f,%d,%f,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d" %( Name,pDockQ,NumRes,IF_plDDT,plDDT1,plDDT2,Diso1[0],Diso1[1],Diso1[2],Diso1[3],Diso2[0],Diso2[1],Diso2[2],Diso2[3],NumResOverlap,len1,len2 ))
#else:
#    print (np.round(pDockQ,3))
