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

    if i==0: # first line
        # get the keys corresponding to these fields
        bead_keys = []
        rb_keys = []
        srb_keys = []

        for k in lndict:
            if lndict[k].startswith('MonteCarlo_Acceptance_BallMover'):
                bead_keys.append(k)

            if lndict[k].startswith('MonteCarlo_Acceptance_P'):
                rb_keys.append(k)

            if lndict[k].startswith('MonteCarlo_Acceptance_Transform'):
                srb_keys.append(k)


    else:
        avg_bead_acceptance = 0.0
        avg_rb_acceptance = 0.0
        avg_srb_acceptance = 0.0

        for k in bead_keys:
            avg_bead_acceptance = avg_bead_acceptance + float(lndict[k])
        avg_bead_acceptance = avg_bead_acceptance/float(len(bead_keys))

        for k in rb_keys:
            avg_rb_acceptance = avg_rb_acceptance + float(lndict[k])
        avg_rb_acceptance = avg_rb_acceptance/float(len(rb_keys))

        for k in srb_keys:
            avg_srb_acceptance = avg_srb_acceptance + float(lndict[k])

        avg_srb_acceptance = avg_srb_acceptance/float(len(srb_keys))

        print ("{:.2f}".format(avg_bead_acceptance), "{:.2f}".format(avg_rb_acceptance), "{:.2f}".format(avg_srb_acceptance), sep='\t')

sf.close()
