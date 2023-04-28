# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys
import re
#from seq import Sequence
from alignedseq import *
from seq import ResidueError
from uniprot import *
#Takes only domain multialigned files
class DomainAlignedSequence(Sequence):
    def __init__(self, fn1=None, ssq1=None):
        sfn1=fn1.strip().split("_")
        if len(sfn1)!=3:
            print("Domain multi alignment file format error!!")
            sys.exit()
        self.dom_start_=int(sfn1[1])
        self.dom_end_=int(sfn1[2])
        super().__init__(sfn1[0], ssq1)
        if not self.is_dna():
            if not self.is_protein():
                print("Sequence not protein or DNA / please check the sequence!!")
    def __getitem__(self,i):
        return self.seq_[i]
    def is_gap(self,pos): #given index returns if there is a gap
        res=True
        #print(pos)
        if self.seq_[pos]!='-':
            res=False
        return res
    '''
    #return -1 if gap is encountered
    def aligned_to_sequence_position(self,apos):
        spos=0
        if self.is_gap(apos):
            spos=-1
        else:
            i=0
            while i < apos:
                if self.seq_[i] != "-":
                    spos=spos+1
                i=i+1
        return spos
    '''
    def aligned_index_to_protein_position(self, apos):
        ppos=self.dom_start_
        i=0
        if self.is_gap(apos):
            ppos = -1
        else:
            while i < apos:
                if not self.is_gap(i):
                    ppos=ppos+1
                    i=i+1
                else:
	                i=i+1
        return ppos
    def get_protein_name(self):
        return self.name_
    def get_entryname_from_uniprot(self, uproteome):
        pn=self.get_protein_name()
        p=uproteome.get_protein_object_by_entryname(pn)
        return p.entry_name_
    def get_domain_start(self):
        return self.dom_start_
    def get_domain_end(self):
        return self.dom_end_
    def is_position_within_domain(self, pos):
        iswithindomain=False
        if pos >= self.dom_start_ and pos <= self.dom_end_:
            iswithindomain=True
        return True
    def get_protein_domain_length(self):
        return (self.dom_end_ - self.dom_start_ + 1)
    def get_domain_length(self):
        return (self.dom_end_ - self.dom_start_ + 1)
    def is_dna(self):
        not re.search(r"[^ATGC-]",self.seq_)
    def is_protein(self):
        if self.is_dna():
            return False
        else:
            for i in self.seq_:
                if i not in ['G','A','V','L','I','P','F','Y','W','S','T','C','M','N','Q','K','R','H','D','E','-','X','Z','B']: #'X' a feature where the identity of the amino acid is unknown (an X is shown at this position in the sequence) and the only information concerning the modification is that the N-terminus is blocked: P80979 (Blocked amino end (Xaa))
#Note: Pyro-Glu is often indicated in papers as ‘pGlu’ and sometimes, in one-letter code as “U”, although this is now used for selenocysteine. In figures of publications, it may be cited as Z, pQ or E
                    raise ResidueError("Residue '%s' cannot be identified as either a nucleotide or amino acid for sequence %s."%(i, self.name_))
            return True
	#Given aligned position returns the actual sequece position
#    def aligned_to_sequence_position(self, apos):
#        spos=0
#        i=0
#        while i < apos:
#            if self.seq_[i] != "-":
#                spos=spos+1
#                i=i+1
#        return spos
	#Given actual sequence position returns the aligned position / column index
    def sequence_position_to_aligned_index(self,spos):
        apos=0
        i=0
        spos=spos-self.dom_start_+1
        #print(spos,self.dom_start_)
        #print(spos)
        while spos > 0:
            if not self.is_gap(i):
                apos = apos+1
                spos=spos-1
                i=i+1
            else:
                apos= apos+1
                i=i+1
        #    print("positions are ",spos,apos)
        #print(apos)
        return apos-1
