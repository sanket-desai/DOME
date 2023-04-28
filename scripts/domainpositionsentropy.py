from aacids import *
import pysam
from domainmultialignedseq import *
from uniprot import *
#Load domain and calculate entropy scores for individual positions
'''
if column contains > 1 residue
Select residue with highest cosmic count in the column (1)
extract all residues conserved with (1) with the column
output needs to be:
    domain count mutsum mutcounts members
'''
domainalndir="../data/domerefdata/final/mafft_outaln/"

class COSMICGenomeWideCountRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[1]
        self.posaa_=int(srec[7])
        self.refaa_=srec[6]
        self.altaa_=srec[8]
        self.cosmicid_=srec[4]
        self.mcount_=int(srec[0])
#This is where we are

def get_position_cosmic_count_records(cosmicdbtbx , name, cpos):
    ccr=[]
    try:
        tbi=cosmicdbtbx.fetch(name, cpos-1, cpos)
        for i in tbi:
            cre=COSMICGenomeWideCountRecord(i)
            if cre.posaa_ == cpos:
                ccr.append( cre )
    except ValueError:
        pass
    return ccr

def get_position_count(ccrvec):
    tmuts=0
    for c in ccrvec:
        tmuts+=c.mcount_
    return tmuts

class DomainMutationCompile(object):
    def __init__(self, dname, domainalndir , cosmicgwtbx, uniobj):
        domalignment=domainalndir+dname+".aln.fa"
        da=DomainMultiAlignedSequence(dname, domalignment)
        aa=AminoAcids()
        self.outvec=[]
        for i in range(0, len(da)): #i is index
            cmap=da.column_map(i) #column map protein->residue
            #print(cmap)
            pposmap={} #protein position map protein->position
            mutmap={} #protein->cosmic mutation count
            filtcmap={}
            filtpposmap={}
            filtmutmap={}
            vfiltcmap={}
            vfiltpposmap={}
            vfiltmutmap={}
            maxmut=0
            highmutprot="" #protein having highest mutation count in a column
            rootaa=""
            #entrynames=[]
            #unid=[]
            totalmuts=0
            vtotalmuts=0
            alignedresuni=[]
            for c in cmap:
                tabrec=[]
                try:
                    entryname=uniobj.accession_entryname_[c]
                    #entrynames.append(entryname)
                    #unid.append(c)
                    ads=da.get_aligned_domain_sequence(c)
                    ppos=ads.aligned_index_to_protein_position(i)
                    #print(entryname)
                    #print(ppos)
                    pposmap[c]=ppos
                    if ppos>0:
                        mutpos=0
                        alignedresuni.append(c)
                        try:
                            tbi=cosmicgwtbx.fetch(entryname, ppos-1, ppos)
                            for ti in tbi:
                                tabrec.append(COSMICGenomeWideCountRecord(ti))
                            mutpos=get_position_count(tabrec)
                        except ValueError:
                            pass
                        mutmap[c]=mutpos
                        if mutpos > maxmut:
                            maxmut=mutpos
                            highmutprot=c
                            rootaa=ads[i]
                except KeyError:
                    pass            
            index=0
            #check following
            for c in alignedresuni:
                if aa.is_similar(rootaa, cmap[c]):
                    #entryname=entrynames[index]
                    filtcmap[c]=cmap[c] #residue
                    filtmutmap[c]=mutmap[c] #mutation
                    totalmuts+=mutmap[c]
                    filtpposmap[c]=pposmap[c] #position
                    index+=1
                else:
                    vfiltcmap[c]=cmap[c] #residue
                    vfiltmutmap[c]=mutmap[c] #mutation
                    vtotalmuts+=mutmap[c]
                    vfiltpposmap[c]=pposmap[c] #position
            if (len(filtcmap)>1 or len(vfiltcmap)) and vtotalmuts+totalmuts > 0:#if multiple values are seen and total mutations > 0
                mutlist=[]
                domain=dname
                protposresent=[]
                for f in filtcmap:
                    res=filtcmap[f]
                    pos=filtpposmap[f]
                    ent=entryname=uniobj.accession_entryname_[f]
                    protposresent.append(ent+";"+str(pos)+";"+res)
                    mutlist.append(str(filtmutmap[f]))
                #v's are all complementary to the highest cosmic count position conserved AA
                vmutlist=[]
                vprotposresent=[]
                for vf in vfiltcmap:
                    vres=vfiltcmap[vf]
                    vpos=vfiltpposmap[vf]
                    vent=uniobj.accession_entryname_[vf]
                    vprotposresent.append(vent+";"+str(vpos)+";"+vres)
                    vmutlist.append(str(vfiltmutmap[vf]))
                #domain  count   mutsum  entropy members
                #outvec.append("%s\t%d\t%d\t%s\t%s" %(domain, len(filtcmap), totalmuts, ",".join(mutlist), "|".join(protposresent)  ) )
                #print("%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s" %(i, domain, len(filtcmap), totalmuts, ",".join(mutlist), "|".join(protposresent), len(vfiltcmap), vtotalmuts, ",".join(vmutlist), "|".join(vprotposresent)  ) )
                self.outvec.append("%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s" %(i, domain, len(filtcmap), totalmuts, ",".join(mutlist), "|".join(protposresent), len(vfiltcmap), vtotalmuts, ",".join(vmutlist), "|".join(vprotposresent)  ) )
            '''
            if len(filtcmap)>1 and totalmuts > 0:#if multiple values are seen and total mutations > 0
                protposresent={}
                mutlist=[]
                domain=dname
                protposresent=[]
                for f in filtcmap:
                    res=filtcmap[f]
                    pos=filtpposmap[f]
                    ent=entryname=uniobj.accession_entryname_[f]
                    protposresent.append(ent+";"+str(pos)+";"+res)
                    mutlist.append(str(filtmutmap[f]))
                #domain  count   mutsum  entropy members
                #outvec.append("%s\t%d\t%d\t%s\t%s" %(domain, len(filtcmap), totalmuts, ",".join(mutlist), "|".join(protposresent)  ) )
                print("%d\t%s\t%d\t%d\t%s\t%s" %(i, domain, len(filtcmap), totalmuts, ",".join(mutlist), "|".join(protposresent)  ) )
            '''
        #print("\t".join(outvec))
        #print(outvec)

def main():
    outfi=open(sys.argv[2],'w')
    #infi=open("../data/domerefdata/final/pfam.ids")
    infi=open(sys.argv[1])
    pfamids=[]
    for i in infi:
        pfamids.append(i.strip())
    infi.close()
    cosmicdbtbx=pysam.TabixFile("../data/domerefdata/final/CosmicMutantExport_missense_genomewide_cor_entry_pchangecount_ord.tsv.gz")
    domainalndir="../data/domerefdata/final/mafft_outaln/"
    #pfids="PF07714"
    uniprotobj=UniprotProteomeParser("../data/uniprot/uniprot_hs_proteome_12_4_2021.txt")
    print("Loaded references!")
    for p in pfamids:
        d=DomainMutationCompile(p, domainalndir, cosmicdbtbx, uniprotobj)
        if len(d.outvec)>0:
            for o in d.outvec:
                outfi.write("%s\n" %(o))
            print("Mutations compiled for %s" %(p))
        else:
            print("No mutations within %s" %(p))
    outfi.close()
if __name__=="__main__":
    main()
