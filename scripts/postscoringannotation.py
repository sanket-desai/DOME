#!/usr/bin/python3
#Assimilate final files from scoring.py and perform following annotations:
#1 extract the significant positions for which analogous positions are listed
#2 extract AA residues (both analogous mutation and statistically significant positions)
#3 if analogous position within domain, extract domain ID and name - add column

import sys
import os
from domainconservationscore import *
from uniprot import *

proteomeswissprotfile="../data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
#genedomainmap="../data/mappings/gene_domain_map.txt"
genedomainmap="../data/uniprot/uniacc_pfam.tsv"
domalnscoredir="../data/domerefdata/final/mafft_cons_score/"
domalndir="../data/domerefdata/final/mafft_outaln/"
functionalsitefile="../data/uniprot/uniprot_hs_proteome_12_4_2021_positions_impsites_withensg_uniq.txt"
hgncuniprotmapfile='../data/mappings/hgnc_uniprot_mapping2.tsv'
tcgastatsfile="../data/tcga_stats/Stats_Sig_Protein_db_all_cancer.tsv"
caddmissense="../data/cadd_download/cadd_missense_sorted.tsv"
phosphositeplusfile="../data/phosphositeplus/phosphositeplus_all_modres_human.tsv"
statanalogpos="../data/scoring_track_files_round3/statsig_analogpos.tsv"
# all the scoring tracks would be of following file format:
# proteinid scores,comma,separated  totalscore length

class FinalScoreFileRecordParser(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.key_=srec[0]
        self.entry_name_=srec[1]
        self.position_=int(srec[2])
        self.statscore_=float(srec[3])
        self.domscore_=float(srec[4])
        self.domrawscore_=float(srec[5])
        self.caddscore_=float(srec[6])
        self.caddrawscore_=float(srec[7])
        self.funcsitescore_=int(srec[8])
    def __str__(self):
        return "\t".join([self.key_, self.entry_name_, str(self.position_), str(self.statscore_), str(self.domscore_), str(self.domrawscore_), str(self.caddscore_), str(self.caddrawscore_), str(self.funcsitescore_)])

class PostScoringAnnotation(object):
    def __init__(self, finalfile, statanaposfile):
        pdr=ProteinDomainRecords(genedomainmap)
        #print("Protein domain records loaded")
        uniprotobj=UniprotProteomeParser(proteomeswissprotfile)
        #print("Uniprot loaded")
        analogkey_statkey_map={}
        askfi=open(statanalogpos)
        for a in askfi:
            sa=a.strip().split("\t")
            analogkey_statkey_map[sa[1]]=sa[0]
        #print("stat analog position file loaded")
        with open(finalfile) as fifi:
            for line in fifi:
                if not line.startswith('gene'):
                    fsr=FinalScoreFileRecordParser(line)
                    anapobj=uniprotobj.get_protein_object_by_entryname(fsr.entry_name_)
                    accs=anapobj.accessions_
                    dacc=''
                    analogres='-'
                    try:
                        analogres=anapobj.sequence_[fsr.position_-1]
                    except IndexError:
                        pass
                    sigres='-'
                    for a in accs:
                        if pdr.has_domains(a):
                            dacc=a
                            break
                    if not dacc == '':
                        if pdr.is_position_within_domain(dacc, fsr.position_):
                            dom=pdr.get_domain_for_protein_position(dacc, fsr.position_)
                            domainname=dom.domainname_
                            domainid=dom.domainid_
                            sigkey='-'
                            sigres='-'
                            try:
                                sigkey=analogkey_statkey_map[fsr.key_]
                                ssigkey=sigkey.split("_")
                                sigentryname=ssigkey[0]+"_"+ssigkey[1]
                                sigpos=int(ssigkey[2])
                                sigpobj=uniprotobj.get_protein_object_by_entryname(sigentryname)
                                sigres=sigpobj.sequence_[sigpos-1]
                            except IndexError:
                                pass
                            except KeyError:
                                pass
                            print(str(fsr)+"\t"+analogres+"\t"+sigkey+"\t"+sigres+"\t"+domainid+"\t"+domainname)
                        else:
                            domainname='-'
                            domainid='-'
                            sigkey='-'
                            sigres='-'
                            print(str(fsr)+"\t"+analogres+"\t"+sigkey+"\t"+sigres+"\t-\t-")
                    else: #if mutation not within domain
                        domainname='-'
                        domainid='-'
                        sigkey='-'
                        sigres='-'
                        print(str(fsr)+"\t"+analogres+"\t"+sigkey+"\t"+sigres+"\t"+domainid+"\t"+domainname)

def main():
    if len(sys.argv)!=3:
        print("Arguments: final_statsig_output.tsv, statsig_analogpos.tsv file")
        sys.exit(0)
    psa=PostScoringAnnotation(sys.argv[1], sys.argv[2])

if __name__=='__main__':
    main()
