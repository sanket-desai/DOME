#!/usr/bin/python3
#extracts functional sites from proteome and mapped to individual proteins in form of a map

import sys
import os

class FunctionalSiteRecord(object):
    def __init__(self,rec):
        srec=rec.strip().split('\t')
        if len(srec) < 9:
            raise FunctionalSiteRecordFormatExcpetion("Record is not in right format: %s" %(rec))
        self.entry_name_=srec[0]
        self.gene_name_=''
        self.synonyms_=[]
        x=srec[1].split(';')
        for i in x:
            i=i.strip()
            if i.startswith('Name='):
                self.gene_name_=i.replace('Name=','').replace(';','').strip()
            elif i.startswith('Synonyms='):
                self.synonyms_=i.replace('Synonyms=','').replace(';','').strip()
        self.uniprot_ids_=srec[2].split(",")
        self.ensg_ids_=srec[3].split(',')
        self.feature_type_=srec[4]
        self.feature_id_=srec[5]
        self.feature_start_=int(srec[6])
        self.feature_end_=int(srec[7])
        self.qualifiers_=srec[8].split(';')
    def has_uniprot_id(self, unid):
        return unid in self.uniprot_ids_
    def has_ensg_id(self, ensg):
        return ensg in self.ensg_ids_

class FunctionalSite(object):
    def __init__(self, funrec):
        self.feature_id_=funrec.feature_id_
        self.feature_start_=funrec.feature_start_
        self.feature_end_=funrec.feature_end_
        self.feature_type_=funrec.feature_type_
        self.feature_qualifiers_=funrec.qualifiers_
    def is_position_within_functional_site(self,pos):
        iswithin=False
        for i in range(self.feature_start_, self.feature_end_+1):
            if pos == i:
                iswithin=True
                break
        return iswithin
    def __len__(self):
        return (self.feature_end_-self.feature_start_)+1

class ProteinFunctionalSites(object):
    def __init__(self, fsitesfile):
        functionalsitefi=open(fsitesfile)
        self.header_=functionalsitefi.readline().strip().split('\t')
        self.entry_name_uniprotids_map_={}
        self.entry_name_ensgids_map_={}
        self.entry_name_gene_name_map_={}
        self.protein_functionalsite_map_={} #protein is entry name key based
        lindex=0
        for l in functionalsitefi:
            lindex+=1
            fsr=FunctionalSiteRecord(l)
            if not fsr.entry_name_ in self.entry_name_ensgids_map_:
                self.entry_name_ensgids_map_[fsr.entry_name_]=fsr.ensg_ids_
                self.entry_name_uniprotids_map_[fsr.entry_name_]=fsr.uniprot_ids_
                self.entry_name_gene_name_map_[fsr.entry_name_]=fsr.gene_name_
            if not fsr.entry_name_ in self.protein_functionalsite_map_:
                fs=FunctionalSite(fsr)
                self.protein_functionalsite_map_[fsr.entry_name_]= [ fs ]
            else:
                fss=self.protein_functionalsite_map_[fsr.entry_name_]
                fs=FunctionalSite(fsr)
                fss.append(fs)
                self.protein_functionalsite_map_[fsr.entry_name_] = fss
        functionalsitefi.close()
    def get_number_of_proteins(self):
        return len(self.protein_functionalsite_map_)
    def get_gene_name(self, unid):
        return self.entry_name_gene_name_map_[unid]
    def get_number_of_functional_sites(self, unid):
        nfunc=0
        if unid in self.protein_functionalsite_map_:
            nfunct=len(self.protein_functionalsite_map_[unid])
        return nfunct
    def get_functional_sites(self, unid):
        fsites=[]
        if unid in self.protein_functionalsite_map_:
            fsites=self.protein_functionalsite_map_[unid]
        return fsites
    def is_position_within_functional_site(self, unid, pos):
        isfunctionalsite=False
        fsi=self.get_functional_sites(unid)
        for f in fsi:
            if f.is_position_within_functional_site(pos):
                isfunctionalsite=True
                break
        return isfunctionalsite

def main():
    p=ProteinFunctionalSites(functionalsitefile)
    print(p.get_number_of_proteins())

if __name__=='__main__':
    main()
