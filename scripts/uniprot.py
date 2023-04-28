#!/usr/bin/python3
#Script to create track scores for indivual mutation residue scoring
#1 domain conservation score for individual position
#2 cadd score for individual mutation
#3 mutated or non-mutated statistically significant (analogous to one of such mutations)
#4 Functionally important positions to be scored as 1 or 0 based on presence or abscence


#Entry name is the key for the map

import sys
import os
import Bio.SwissProt as sp

proteomeswissprotfile="../data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
class Protein(object):
    def __init__(self, prec):
        #gn=prec.gene_name[0:prec.gene_name.find(';')]
        gn=prec.gene_name
        sgn=gn.strip().split(';')
        self.genenames_=[]
        for s in sgn:
            s=s.strip()
            ss=s
            if s.find('Name=') > -1:
                s=s.replace('Name=','')
                if s.find("{") > -1:
                    s=s[:s.find(" {")]
                if not s == '':
                    self.genenames_.append(s)
            if ss.find("Synonyms=") > -1:
                ss=ss[ss.find("Synonyms="):]
                ss=ss.replace('Synonyms=','')
                sss=ss.strip().split(',')
                for ass in sss:
                    ass=ass.strip()
                    if ass.find(";"):
                        ass=ass.replace(";", '')
                    if ass.find("{") > -1:
                        ass=ass[:ass.find(" {")]
                    if not ass == '':
                        self.genenames_.append(ass)
        self.accessions_=prec.accessions
        self.sequence_=prec.sequence
        self.entry_name_=prec.entry_name
        self.phosphositeplusids_=[]
        self.ensgids_=[]
        for cr in prec.cross_references:
            if cr[0].find('Ensembl') != -1:
                for c in cr:
                    if c.find('ENSG') != -1:
                        en=c
                        if c.find('.') != -1:
                            en=c[0:c.find('.')]
                        self.ensgids_.append(en.strip())
            elif cr[0].find('PhosphoSitePlus') != -1:
                self.phosphositeplusids_.append(cr[1])
        tempen=set(self.ensgids_)
        self.ensgids_=list(tempen)
    def __len__(self):
        return len(self.sequence_)

class UniprotProteomeParser(object):
    def __init__(self,fn):
        self.protein_uniprot_info_={}
        self.genename_entryname_={} #created to bind tcga and uniprot records
        self.accession_entryname_={} #created for uniprot and domain record mapping
        self.ensgid_entryname_={}
        self.phosphositeplusid_entryname_={}
        with open(fn) as handle:
            records=sp.parse(handle)
            for r in records:
                self.protein_uniprot_info_[r.entry_name]=Protein(r)
        handle.close()
        for p in self.protein_uniprot_info_:
            info=self.protein_uniprot_info_[p]
            for g in info.genenames_:
                self.genename_entryname_[g]=p
            #self.genename_entryname_[info.gene_name_]=p
            accessions=info.accessions_
            pspids=info.phosphositeplusids_
            for pi in pspids:
                self.phosphositeplusid_entryname_[pi]=p
            for a in accessions:
                self.accession_entryname_[a]=p
            for e in info.ensgids_:
                self.ensgid_entryname_[e]=info.entry_name_
    def get_protein_length(self, entryname):
        return len(self.protein_uniprot_info_[entryname])
    def get_accessions_of_protein(self, entryname):
        pro=self.protein_uniprot_info_[entryname]
        return pro.accessions_
    def get_number_of_proteins(self):
        return len(self.protein_uniprot_info_)
    def get_protein_object_by_genename(self, gene):
        return self.protein_uniprot_info_[self.genename_entryname_[gene]]
    def get_protein_object_by_phosphositeplusid(self, pspid):
        pro=None
        try:
            pro=self.protein_uniprot_info_[self.phosphositeplusid_entryname_[pspid]]
        except KeyError:
            try:
                pro=self.protein_uniprot_info_[self.accession_entryname_[pspid]]
            except KeyError:
                raise KeyError
        return pro

    def get_protein_object_by_entryname(self, uacc):
        return self.protein_uniprot_info_[uacc]
    def get_protein_entrynames(self):
        return self.protein_uniprot_info_.keys()
    def get_protein_object_by_accession(self, uacc):
        return self.protein_uniprot_info_[self.accession_entryname_[uacc]]
    def get_protein_object_by_ensgid(self, en):
        return self.protein_uniprot_info_[self.ensgid_entryname_[en]]
    def get_residue_by_entry(self, entry, position):
        p=self.get_protein_object_by_entryname(entry)
        return p.sequence_[position-1] #index position match
    def sequence_residue_check_by_entry(self, entry, position, residue):#residue given by user should match with the uniprot proteome residue stored in the datastructure
        p=self.get_protein_object_by_entryname(entry)
        return p.sequence_[position-1]==residue #index position match
    def is_entryname_within_uniprot_proteome(self, entryn):
        return entryn in self.protein_uniprot_info_.keys()
'''
def main():
    p=UniprotProteomeParser(proteomeswissprotfile)
    enames=p.get_protein_entrynames()
    for e in enames:
        pobj=p.get_protein_object_by_entryname(e)
        gn=pobj.genenames_[0]
        accs=pobj.accessions_
        print("%s\t%s" %(gn, ','.join(accs)))
if __name__=='__main__':
    main()
'''
