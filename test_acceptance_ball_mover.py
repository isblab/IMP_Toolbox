
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.crosslinking
import os,sys,string

statFile = sys.argv[1]

# read the stat file
sf =open(statFile,'r')


for i,ln in enumerate(sf.readlines()):
    lndict = eval(ln.strip())
    
    if i==0: 
        keysAsked = []		

        for k in lndict:
            if lndict[k].startswith('MonteCarlo_Acceptance_BallMover'):
                
             print(lndict[k])
             keysAsked.append(k)
        
        print("Keys",keysAsked)

    else:
        avgAcceptance =0.0
        
        for k in keysAsked:
            avgAcceptance = avgAcceptance + float(lndict[k])
    
        avgAcceptance = avgAcceptance/float(len(keysAsked))
        
        print(avgAcceptance)

sf.close()

