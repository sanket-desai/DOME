import sys
import pysam
from aacids import *
#mutation table file format:
#key	gene	pos	statscore	sigdomscores	sigdomrawscores	sigcaddscores	sigcaddrawscores	sigfuncsitescores	residue	statsigmut	statsigres	domainid	domainname
#WNK2_HUMAN_252	WNK2_HUMAN	252	0.6	0.228172977970538	138648	1.11692307274047	3.88770483333333	0	P	MP2K1_HUMAN_124	P	PF00069	Pkinase

cfile="../data/domerefdata/cadd_missense_aref_aalt_ens_gene_apos_craw_cphred2_sorted_entry.tsv.gz"

def total(arr):
    t=0
    for a in arr:
        t+=a
    return t

class CADDMissenseScoreFileRecord(object):
    def __init__(self, rec ):
        srec=rec.strip().split("\t")
        if len(srec) != 9:
            print("Format error at record : %s" %(rec.strip()) )
            sys.exit(0)
        self.entryname_=srec[0]
        self.aaref_=srec[1]
        self.aachange_=srec[2]
        self.ensgid_=srec[3]
        self.enstid_=srec[4]
        self.genename_=srec[5]
        self.aapos_=int(srec[6])
        self.caddrawscore_=float(srec[7])
        self.caddphredscore_=float(srec[8])

def main():
    caddfile=cfile
    muttable=sys.argv[1]
    mutfi=open(muttable)
    mutheader=mutfi.readline().strip().split("\t")
    ctab=pysam.TabixFile(cfile)
    aa=AminoAcids()
    print( "%s\t%s\t%s" %("\t".join(mutheader), "sim_aa_caddphred", "diff_aa_caddphred" ) )
    for m in mutfi:
        #print(m)
        sm=m.strip().split("\t")
        entry=sm[1] # uniprot entry
        pos=sm[2]   # position
        result=ctab.fetch("%s:%s-%s" %(entry, pos, pos))
        simaascore=[]
        diffaascore=[]
        while True:
            try:
                r=next(result)
                csrec=CADDMissenseScoreFileRecord(r)
                if aa.is_similar(csrec.aaref_, csrec.aachange_):
                    simaascore.append(csrec.caddphredscore_)
                else:
                    diffaascore.append(csrec.caddphredscore_)
            except StopIteration:
                break
        meansimscore=0
        meandiffscore=0
        if len(simaascore)>0:
            meansimscore=total(simaascore)/len(simaascore)
        if len(diffaascore) > 0:
            meandiffscore=total(diffaascore)/len(diffaascore)
        #similar mutation score, different mutation score
        print("%s\t%f\t%f" %(m.strip(), meansimscore, meandiffscore ))
    mutfi.close()
    ctab.close()

if __name__=='__main__':
    main()
