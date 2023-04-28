#!/usr/bin/python3
#parse the statistically significant mutations in TCGA and map analogous positions

import sys
import os

#gene.name/protein.name uniprot accession map
tcga_genename_uniprot_accession_map="../data/uniprot/gene_uniaccs.map"

class AccessionTCGAGeneNameMap(object):
    def __init__(self, tcga_gene_unid_map):
        fi=open(tcga_gene_unid_map)
        self.accession_genename_map_={} #accession - > tcga gene name
        self.genename_accession_map_={} #genename -> [accessions]
        accs=[]
        gn=''
        for i in fi:
            si=i.strip().split("\t")
            try:
                gn=si[0]
                ssi=si[1].split(',')
                self.genename_accession_map_[gn]=ssi
                for s in ssi:
                    self.accession_genename_map_[s]=gn
            except IndexError as ie:
                print(i.strip())
    def get_genename(self, unid): #returns individual values
        return self.accession_genename_map_[unid]
    def get_uniprot_accessions(self, gname): #returns an array
        return self.genename_accession_map_[gname]

class HotspotMutation(object):
    def __init__(self, rec):
        srec=rec.strip().split('\t')
        #self.gene_=srec[0]
        self.position_=int(srec[1])
        self.domain_id_=srec[2]
        #self.protein_fdr_=float(srec[3])
        self.protein_fdr_=srec[3]
        #in the newer version qvalue has been replaced by H - hotspot, R - resistance mutation label
        self.tumour_type_=srec[4]
    def is_within_domain(self):
        return self.domain_id_.startswith("PF")

class ProteinHotspotMutations(object):
    def __init__(self,hotfile):
        self.accession_genename_obj_=AccessionTCGAGeneNameMap(tcga_genename_uniprot_accession_map)
        self.protein_hotspot_mutations_={} # protein/gene - > [ hotspot mutation ]
        fi=open(hotfile)
        self.header_=fi.readline().strip().split('\t')
        for l in fi:
            gene=l.strip().split('\t')[0]
            hmut=HotspotMutation(l)
            if not gene in self.protein_hotspot_mutations_:
                self.protein_hotspot_mutations_[gene]=[ hmut ]
            else:
                temp=self.protein_hotspot_mutations_[gene]
                temp.append(hmut)
                self.protein_hotspot_mutations_[gene]=temp
    def get_hotspot_mutations(self, genename):
        return self.protein_hotspot_mutations_[genename]
    def get_hotspot_mutations_by_accession(self, acc):
        return self.get_hotspot_mutations(self.accession_genename_obj_.get_genename(acc))
    def has_hotspot_mutations(self, genename):
        return genename in self.protein_hotspot_mutations_
    def has_hotspot_mutations_by_accession(self, acc):
        hashotspots=False
        try:
            genename=self.accession_genename_obj_.get_genename(acc)
            if genename in self.protein_hotspot_mutations_:
                hashotspots=True
        except KeyError:
            hashotspots=False
        return hashotspots

def main():
    print("hello world")

if __name__=='__main__':
    main()
