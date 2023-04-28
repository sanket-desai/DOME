import sys
import os.path
from domainmultialignedseq import *
from scipy import stats
import pandas as pd

class DomainRecord(object):
    def __init__(self, domainid, domename, start, end):
        self.domainid_=domainid
        self.domainname_=domename
        self.domainstart_=int(start)
        self.domainend_=int(end)
    def is_position_within_domain(self, pos):
        return pos in range(self.domainstart_, self.domainend_+1)

class DomainRecords(object):
    def __init__(self, genedomainmap):
        self.domainrecords_=pd.read_csv( genedomainmap , header=0, sep="\t")
    def get_domain_name_from_id(self, did):
        dname="-"
        df=self.domainrecords_[self.domainrecords_.pfamid == did ]
        if len(df) > 0:
            dname=df.pfamname
        return dname

class ProteinDomainRecords(object):
    def __init__(self, genedomainmap):
        f=open(genedomainmap)
        self.header_=f.readline()
        self.protein_domain_map_={} #protein -> [ domain1, domain2 ]
        #self.uniprotid_genename_={}
        lindex=0
        for l in f:
            #print(lindex)
            lindex+=1
            try:
                sl=l.strip().split("\t")
                #genename=sl[0]
                uniprot=sl[0]
                #print(len(sl))
                #print(sl) #File is not uniform
                #if not uniprot in self.uniprotid_genename_:
                #    self.uniprotid_genename_[uniprot]=genename
                #if sl[2] != 'NA':
                #dr=DomainRecord(sl[2], sl[3], sl[4], sl[5])
                dr=DomainRecord(sl[3], sl[4], sl[1], sl[2])
                if not uniprot in self.protein_domain_map_:
                    self.protein_domain_map_[uniprot]=[dr]
                else:
                    adr=self.protein_domain_map_[uniprot]
                    adr.append(dr)
                    self.protein_domain_map_[uniprot]=adr
            except IndexError as ie:
                print(ie)
                sys.exit(0)
        f.close()
    def get_number_of_proteins(self):
        return len(self.protein_domain_map_)
    #def get_gene_name(self, unid):
    #    return self.uniprotid_genename_[unid]
    def get_number_of_domains(self, unid):
        return len(self.protein_domain_map_[unid])
    def get_domains_for_protein(self, unid):
        return self.protein_domain_map_[unid]
    def is_position_within_domain(self, unid, pos):
        drecs=self.get_domains_for_protein(unid)
        withindomain=False
        for d in drecs:
            withindomain=d.is_position_within_domain(pos)
            if withindomain == True:
                break
        return withindomain
    def get_domain_for_protein_position(self, unid, pos):
        drecs=self.get_domains_for_protein(unid)
        dom=None
        for d in drecs:
            withindomain=d.is_position_within_domain(pos)
            if withindomain == True:
                dom=d
        return dom
    def is_overlapping(self, unid, posstart, posend):
        overlapping=False
        for i in range(posstart, posend+1):
            if is_position_within_domain(unid, i):
                overlapping=True
                break
        return overlapping
    def has_domains(self, unid):
        hasdomains=False
        if unid in self.protein_domain_map_ :
            if len(self.protein_domain_map_[unid]) > 0:
                hasdomains=True
        return hasdomains

class DomainConservationValueError(Exception):
        def __init__(self, value):
                self.value = value
        def __str__(self):
                return repr(self.value)
#name, conservation score file, alignment file
class DomainConservation(object):
    def __init__(self, dname, fn=None, fn2=None):
        self.header_=['domainalnpos','consresidue','consscore']
        fi=open(fn)
        x=fi.readline()
        self.domain_multi_alignment_=DomainMultiAlignedSequence(dname, fn2)
        self.position_residue_score_=[]
        for l in fi:
            sl=l.strip().split('\t')
            if len(sl) != 3:
                print("Error in Domain Conservation Score file format!")
                sys.exit(0)
            self.position_residue_score_.append(sl)
    def protein_position_residue_conservation_score_list(self, unid):
        pos_res_score_list=[]
        if not self.domain_multi_alignment_.has_sequence(unid):
            pass
            #raise DomainConservationValueError("Protein ID %s not found in the domain alignment file %s" %(unid, self.domain_multi_alignment_.domain_name_))
        else:
            domain_seq=self.domain_multi_alignment_.get_aligned_domain_sequence(unid)
            try:
                prange=range( domain_seq.get_domain_start(), domain_seq.get_domain_end()+1)
                for p in prange:
                    aind=domain_seq.sequence_position_to_aligned_index(p)
                    cold=self.position_residue_score_[aind]
                    colscore=float(cold[2])
                    colres=cold[1]
                    pos_res_score_list.append( (p, colres, colscore) )
            except IndexError as ie:
                pass #when the domain indices are beyond the protein positions / may be isoform
        return pos_res_score_list
    def domain_conservationscore_list(self, unid):
        pcs=self.protein_position_residue_conservation_score_list(unid)
        pps=[]
        for p in pcs:
            pps.append(p[2])
        return pps

    def domain_conservationscore_zscore_list(self, unid):
        pcs=self.domain_conservationscore_list(unid)
        return stats.zscore(pcs)
