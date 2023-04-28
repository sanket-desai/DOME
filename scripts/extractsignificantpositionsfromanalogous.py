import sys
import os
from functionalsitescore import *
from domainconservationscore import *
from statisticallysignificantmutations import *
from uniprot import *
from caddrecordparser import *
from phosphositemodrescoring import *
import pysam

proteomeswissprotfile="../data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
#genedomainmap="../data/mappings/gene_domain_map.txt"
genedomainmap="../data/uniprot/uniacc_pfam.tsv"
domalnscoredir="../data/domerefdata/final/mafft_cons_score/"
domalndir="../data/domerefdata/final/mafft_outaln/"
functionalsitefile="../data/uniprot/uniprot_hs_proteome_12_4_2021_positions_impsites_withensg_uniq.txt"
phosphositeplusfile="../data/phosphositeplus/phosphositeplus_all_modres_human.tsv"
hgncuniprotmapfile='../data/mappings/hgnc_uniprot_mapping2.tsv'
tcgastatsfile="../data/tcga_stats/Stats_Sig_Protein_db_all_cancer.tsv"
caddmissense="../data/cadd_download/cadd_missense_sorted.tsv"
finalscoreoutfile='/media/param/DATA/code/python/lab/domain_mutations/dome2021/data/final_prot_sign_scores.tsv'
tcga_genename_uniprot_accession_map="../data/uniprot/gene_uniaccs.map"

def main():
    pdr=ProteinDomainRecords(genedomainmap)
    uniobj=UniprotProteomeParser(proteomeswissprotfile)
    finalfile=open(finalscoreoutfile)
    header=finalfile.readline()
    acctcgamap=AccessionTCGAGeneNameMap(tcga_genename_uniprot_accession_map)
    hotspots=ProteinHotspotMutations(tcgastatsfile)
    protein_hotpositions_map={}
    #run throught the hotspot file and store statistically significant positions / uniprot entry name -> positions
    for h in hotspots.protein_hotspot_mutations_:
        try:
            accs=acctcgamap.get_uniprot_accessions(h)
            for a in accs:
                hpos=[]
                hmuts=hotspots.get_hotspot_mutations_by_accession(a)
                pro=uniobj.get_protein_object_by_accession(a)
                if not pro.entry_name_ in protein_hotpositions_map:
                    for h in hmuts:
                        if not h.position_ in hpos:
                            hpos.append(h.position_)
                    protein_hotpositions_map[pro.entry_name_]=hpos
                else:
                    hpos=protein_hotpositions_map[pro.entry_name_]
                    for h in hmuts:
                        if not h.position_ in hpos:
                            hpos.append(h.position_)
                    protein_hotpositions_map[pro.entry_name_]=hpos
        except KeyError as ke:
            print(ke)
            pass
    print(protein_hotpositions_map)
    outfo=open('../data/final_prot_sign_scores_sigpos.tsv','w')
    for fl in finalfile:
        outputlist=[]
        sfl=fl.strip().split('\t')
        oline="\t".join(sfl)
        entryname=sfl[0]
        prot=uniobj.get_protein_object_by_entryname(entryname)
        siggene=''
        sigpos=0
        for a in prot.accessions_:
            if pdr.has_domains(a):
                ps=int(sfl[1])
                if pdr.is_position_within_domain(a, ps):
                    d=pdr.get_domain_for_protein_position(a,ps)
                    #load domain d and check which sig mutation is in the column of the position
                    domainfilename=domalndir+d.domainid_+".aln.fa"
                    dmsa=DomainMultiAlignedSequence(d.domainid_ , domainfilename)
                    try:
                        amap=dmsa.get_analogous_positions(a, ps)
                        for aa in amap:
                            aprot=uniobj.get_protein_object_by_accession(aa)
                            if aprot.entry_name_ in protein_hotpositions_map:
                                hposs=protein_hotpositions_map[aprot.entry_name_]
                                for hh in hposs:
                                    if not hh == amap[aa]:
                                        sigpos=hh
                                        siggene=aprot.entry_name_
                        if sigpos != 0:
                            sfl.append(siggene)
                            sfl.append(str(sigpos))
                            outputlist.append(",".join(sfl))
                            siggene=''
                            sigpos=0
                    except IndexError as ie:
                        print(ie)
                        pass
                    except Exception as e:
                        print(e)
                        pass
        for o in outputlist:
            outfo.write("%s\n" %(o))

if __name__=="__main__":
    main()
