import sys
from domainpositionsentropy import *

class AminoAcids(object):
    def __init__(self):
        self.positivelycharged_=["K", "R", "H"]
        self.polar_=["S", "T", "C" ,"P" ,"N" ,"Q"]
        self.nonpolar_=["G", "A", "V" ,"L" ,"I" ,"M"]
        self.negativelycharged_=["D", "E"]
        self.aromatic_=["F" ,"Y" ,"W"]
        self.aminoacids_=self.positivelycharged_+self.polar_+self.nonpolar_+self.negativelycharged_+self.aromatic_
        self.aminoacidtypes_=['positivelycharged','negativelycharged','polar','nonpolar','aromatic']
        self.type_aminoacids_={}
        self.type_aminoacids_['positivelycharged']=["K", "R", "H"]
        self.type_aminoacids_["polar"]=["S", "T", "C" ,"P" ,"N" ,"Q"]
        self.type_aminoacids_["nonpolar"]=["G", "A", "V" ,"L" ,"I" ,"M"]
        self.type_aminoacids_["negativelycharged"]=["D", "E"]
        self.type_aminoacids_["aromatic"]=["F" ,"Y" ,"W"]
    def get_aminoacids_from_type(self, atype):
        return self.type_aminoacids_[atype]
    def get_aminoacid_types(self):
        return list(self.type_aminoacids_.keys())
    def get_type_of_aminoacid(self, aa):
        atype=""
        for t in self.type_aminoacids_:
            ta=self.type_aminoacids_[t]
            if aa in ta:
                atype=t
                break
        return atype
    def is_similar(self, aa1, aa2):
        return self.get_type_of_aminoacid(aa1)==self.get_type_of_aminoacid(aa2)
    def get_amino_frequency(self, aalist):
        #print(aalist)
        aacounts={}
        aafreqs={}
        for aa in self.aminoacids_:
            aacounts[aa]=0
        for a in aalist:
            aacounts[a]=aacounts[a]+1
            aafreqs[a]=0
        for aac in aacounts:
            aafreqs[aac] = aacounts[aac]/len(aalist)
        return aafreqs

'''
Aminoacid type -> counts (per position)
Aminoacid type -> AminoAcid type -> count
'''

class AminoAcidCosmicChangeFrequency(object): #For each column recort, it will create the following
    def __init__(self, li, cosmicdbtbx, aminoacidobj ): #li contains list of residues in following format [ "EGFR_HUMAN;790;T", ... ]
        A=aminoacidobj
        self.aminoacid_type_counts_={} #Amino acid type and counts of each type / especially within column
        self.aminoacid_changetype_counts_={} # Aminoacid type-> { changedtype -> count }
        self.aminoacid_changetype_frequency_={}
        for a in A.get_aminoacid_types():
            self.aminoacid_type_counts_[a]=0
            self.aminoacid_changetype_counts_[a]={} #store change aatype and their counts
        #write code to extract data from cosmic
        poscosmcountrecs=[]
        for l in li:
            sl=l.split(';')
            if len(sl) != 3:
                print("Bad input in the protein, position, residue list: %s" %(l))
                #sys.exit(0)
            if len(sl[2])==1:
                name=sl[0]
                cpos=int(sl[1])
                aatype=A.get_type_of_aminoacid(sl[2])
                try:
                    aacount=self.aminoacid_type_counts_[aatype]
                    self.aminoacid_type_counts_[aatype]=aacount+1
                    poscosmcountrecs=poscosmcountrecs+get_position_cosmic_count_records( cosmicdbtbx, name,cpos)
                except KeyError:
                    pass
        #print(len(poscosmcountrecs))
        #print("poscosmcountrecs")
        for p in poscosmcountrecs:
            reftype=A.get_type_of_aminoacid(p.refaa_)
            changeaatype=A.get_type_of_aminoacid(p.altaa_)
            try:
                chcounts=self.aminoacid_changetype_counts_[reftype]
                if changeaatype in chcounts:
                    chcounts[changeaatype]=chcounts[changeaatype]+p.mcount_
                else:
                    chcounts[changeaatype]=p.mcount_
                self.aminoacid_changetype_counts_[reftype]=chcounts
            except KeyError:
                pass
        for rt in self.aminoacid_changetype_counts_:
            ct=self.aminoacid_changetype_counts_[rt]
            changect={}
            totalcmuts=0
            for c in ct:
                totalcmuts+=ct[c]
            for c in ct:
                changect[c]=ct[c]/totalcmuts
            self.aminoacid_changetype_frequency_[rt]=changect
    def get_typechange_frequency(self, refaa, altaa, aminoacidobj): #using amino acids revert with the frequency
        x=self.aminoacid_changetype_frequency_[ aminoacidobj.get_type_of_aminoacid(refaa)]
        try:
            return x[aminoacidobj.get_type_of_aminoacid(altaa)]
        except KeyError:
            return 0

cosmicdbtbx="../data/domerefdata/final/CosmicMutantExport_missense_genomewide_cor_entry_pchangecount_ord.tsv.gz"
protdomainindextbx="../data/domerefdata/final/proteome_domain_entropy_cosmicgw.tsv.gz"

def main():
    cosmictbx=pysam.TabixFile(cosmicdbtbx)
    protdomtbx=pysam.TabixFile(protdomainindextbx)
    a=AminoAcids()
    #x="ABL1_HUMAN;315;T|CSF1R_HUMAN;663;T|TGFR2_HUMAN;325;T|TYK2_HUMAN;978;T|TEX14_HUMAN;326;Y|TESK1_HUMAN;129;T|TNI3K_HUMAN;539;T|RIPK2_HUMAN;95;T|TEC_HUMAN;442;T|BRAF_HUMAN;529;T|TXK_HUMAN;343;T|ACVR1_HUMAN;283;T|EPHAA_HUMAN;723;T|EPHA8_HUMAN;713;T|EPHB1_HUMAN;697;T|BLK_HUMAN;312;T|EPHB4_HUMAN;693;T|CSK_HUMAN;266;T|EPHB3_HUMAN;711;T|EPHA4_HUMAN;699;T|ACK1_HUMAN;205;T|BTK_HUMAN;474;T|TESK2_HUMAN;130;T|YES_HUMAN;348;T|SRC_HUMAN;341;T|SG196_HUMAN;149;T|SRMS_HUMAN;302;T|BMR1A_HUMAN;309;T|BMR1B_HUMAN;279;T|ACV1C_HUMAN;270;S|DDR2_HUMAN;654;T|LIMK1_HUMAN;413;T|NRBP_HUMAN;150;T|PGFRB_HUMAN;681;T|RAF1_HUMAN;421;T|ARAF_HUMAN;382;T|NTRK1_HUMAN;589;F|NTRK3_HUMAN;617;F|NTRK2_HUMAN;617;F|ILK_HUMAN;269;T|KSR2_HUMAN;739;T|PGFRA_HUMAN;674;T|FGR_HUMAN;334;T|PTK6_HUMAN;264;T|HCK_HUMAN;333;T|LIMK2_HUMAN;405;T|ABL2_HUMAN;361;T|GUC2D_HUMAN;619;S|GUC2F_HUMAN;623;T|LCK_HUMAN;316;T|KIT_HUMAN;670;T|MUSK_HUMAN;655;F|M3K20_HUMAN;82;T|EPHB6_HUMAN;748;T|ANPRA_HUMAN;612;T|LYN_HUMAN;319;T|JAK3_HUMAN;902;Q|JAK2_HUMAN;929;Q|ROR1_HUMAN;552;F|ROR2_HUMAN;553;F|LMTK2_HUMAN;214;F|EPHA1_HUMAN;702;T|EPHB2_HUMAN;699;T|ITK_HUMAN;435;F|ANPRB_HUMAN;596;T|IRAK4_HUMAN;262;Y|JAK1_HUMAN;956;E|KSR1_HUMAN;686;T|EPHA3_HUMAN;699;T|FYN_HUMAN;342;T|DDR1_HUMAN;701;T|ERBB2_HUMAN;798;T|EGFR_HUMAN;790;T|FRK_HUMAN;306;T|EPHA2_HUMAN;692;T|EPHA5_HUMAN;753;T|FLT3_HUMAN;691;F|BMX_HUMAN;489;T|ERBB4_HUMAN;796;T|ACVL1_HUMAN;277;T|ERBB3_HUMAN;787;T"
    #x=x.split("|")
    pfam="PF00010"
    index=59
    for p in protdomtbx.fetch(pfam, index-1, index):
        sp=p.strip().split("\t")
        if sp[0]==pfam+";"+str(index):
            x=sp[6]+"|"+sp[10]
            x=x.split("|")
            aa=AminoAcidCosmicChangeFrequency(x, cosmictbx, a)
            print("%.4f" %(aa.get_typechange_frequency("T","T", a)))
if __name__=="__main__":
    main()
