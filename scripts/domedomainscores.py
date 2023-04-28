from domainmultialignedseq import *
from uniprot import *

#inputs domain alignment file, dome output file (format below)
#returns dome scores for each protein position in the domain - as map of protein->score list
#dome output format:
#tier    protein mutpos  refaa   altaa   domain  statscore       analogtomut     caddscorcosmiccount      colcoscount     colnumofresidues        clinvarsig      uniprotsite     phosphosite      uniprotmutagen  hotspot3dann    X3d_dom_close   X3d_statsig_close       X3d_analog_close X3d_uniprotsite_close   X3d_pspsite_close       X3d_hotspot_close       interfaceposition        contextscore    neutralityscore entropyscore    alterfrequencyscore      biochemcontext  meanscore       mean4   logcosm
#example file location: /mnt/d/code/python/lab/domain_mutations/dome2021/data/domerefdata/io/final/tcga_analysis/tcga_ptyrkin_cadd20.tsv

#algorithm
#for each scored position, check if the protein position exists within the domain, convert position to index and store the score

proteomeswissprotfile="../data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
pfamfile="../data/domerefdata/io/final/tcga_analysis/PF07714.aln.fa"
domescorefile="/mnt/d/code/python/lab/domain_mutations/dome2021/data/domerefdata/io/final/tcga_analysis/tcga_ptyrkin_cadd20.tsv"
class DOMEDomainScores(object):
    def __init__(self):
        multidomain=DomainMultiAlignedSequence("PF07714",pfamfile)
        up=UniprotProteomeParser(proteomeswissprotfile)
        domprotentrynames=[]
        domprotscores={}
        index_cosmcount={}
        for l in range(0,len(multidomain)):
            index_cosmcount[l]=0
        for i in range(0,multidomain.number_of_sequences()):
            try:
                asq=multidomain[i]
                p=up.get_protein_object_by_accession(asq.name_)
                domprotentrynames.append(p.entry_name_)
                domprotscores[p.entry_name_]=[0]*len(multidomain)
            except IndexError:
                pass
                #print(i)
        #print(domprotentrynames)
        fi=open(domescorefile)
        h=fi.readline().strip().split("\t")
        for i in fi:
            si=i.strip().split("\t")
            sc=float(si[30]) #score
            pos=int(si[2])
            prot=si[1]
            cosmiccolcount=int(si[10])
            changed=False
            if prot in domprotentrynames: #protein contains the domain in question
                asqq= multidomain[domprotentrynames.index(prot) ]
                if asqq.is_position_within_domain(pos):
                    aind=asqq.sequence_position_to_aligned_index(pos)
                    pscores=domprotscores[prot]
                    index_cosmcount[aind]=cosmiccolcount
                    if not asqq.is_gap(aind):
                        pscores[aind]=sc
                    else:
                        pscores[aind]=-1
                    domprotscores[prot]=pscores
        length=[]
        ccount=[]
        for l in range(0,len(multidomain)):
            length.append(str(l))
            ccount.append(str(index_cosmcount[l]))
        print("%s" %("\t".join(ccount)))
        print("%s" %("\t".join(length)))
        for d in domprotscores:
            scores=domprotscores[d]
            sscores=[]
            for s in scores:
                sscores.append(str(s))
            print("%s\t%s" %(d, "\t".join(sscores)))

def main():
    dobj=DOMEDomainScores()

if __name__=="__main__":
    main()
