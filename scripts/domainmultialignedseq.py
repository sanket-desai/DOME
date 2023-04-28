import sys
import os.path
from domainalignedseq import *
from proteinmutationmapper import *

class ProteinNotFoundError(Exception):
    def __init__(self, value):
        super().__init__(value)

class DomainMultiAlignedSequence(object):
    def __init__(self, dname, fn=None):
        self.domain_name_=dname
        self.aseqs=[]
        if fn is None:
            self.aseqs=[]
        elif fn != None and os.path.exists(fn):
            fi = open(fn,'r')
            rec=[]
            b=False
            for ln in fi:
                fl = ln.rstrip()
                if fl.startswith(">") and b:
                    sn,ss="",""
                    for i in rec:
                        if i.startswith(">"):
                            sn=i[1:]
                        else:
                            ss=ss+i.replace(" ","")
                    asq = DomainAlignedSequence(sn,ss)
                    del rec[:]
                    self.aseqs.append(asq)
                    rec.append(fl)
                elif fl.startswith(">") and not b and not rec:
                    b=True
                    rec.append(fl)
                elif not fl.startswith(">"):
                    rec.append(fl)
            sn,ss="",""
            for i in rec:
                if i.startswith(">"):
                    sn=i[1:]
                else:
                    ss=ss+i.replace(" ","")
            asq = DomainAlignedSequence(sn,ss)
            del rec[:]
            self.aseqs.append(asq)
        else:
            raise FileNotFoundError("File",fn,"does not exist!!")
            #print( "File",fn,"does not exist!!")
            #sys.exit()
    def __len__(self):
        return len(self.aseqs[0])
    def __getitem__(self,i):
        return self.aseqs[i]
    def __setitem__(self,i,j):# i is the index, j is the AlignedSeq object
        self.aseqs[i]=j
    def column(self, i): #Given index return the column in the msa
        col = []
        for as_ in self.aseqs:
            col.append(as_[i])
        return col
    def column_map(self,i):
        colmap={}
        for a in self.aseqs:
            colmap[a.name_]=a.seq_[i]
        return colmap
    def number_of_sequences(self):
        return len(self.aseqs)
    def has_sequence(self, prname):
        hseq=False
        for a in self.aseqs:
            if a.get_protein_name() == prname:
                hseq=True
                break
        return hseq
    def get_aligned_domain_sequence(self, prname):
        ads=None
        if self.has_sequence(prname):
            for a in self.aseqs:
                if a.get_protein_name() == prname:
                    ads=a
                    break
        else:
            raise ProteinNotFoundError("%s not found in the DomainMultiAlignedSequence object of %s !!" %(prname, self.domain_name_))
        return ads
    #returns {protein->[ analogous_position, analogous_residue ] }
    '''
    used in the original code, changing MutationAA to only position
    def get_analogous_positions(self, prname, pquery ):
        apmap={}
        qaseq=self.get_aligned_domain_sequence(prname)
        aamu=MutationAA(pquery)
        aind=qaseq.sequence_position_to_aligned_index(aamu.position_)
        for a in self.aseqs:
            apmap[a.name_]=a.aligned_index_to_protein_position(aind)
        return apmap
    '''
    def get_analogous_positions(self, prname, position ):
        apmap={}
        qaseq=self.get_aligned_domain_sequence(prname)
        aind=qaseq.sequence_position_to_aligned_index(position)
        for a in self.aseqs:
            apmap[a.name_]=a.aligned_index_to_protein_position(aind)
        return apmap
    def get_analogous_residues(self, prname, pquery):
        apmap={}
        qaseq=self.get_aligned_domain_sequence(prname)
        aamu=MutationAA(pquery)
        aind=qaseq.sequence_position_to_aligned_index(aamu.position_)
        return self.column_map(aind)
    def get_index_position(self, prname, pquery):
        qaseq=self.get_aligned_domain_sequence(prname)
        aamu=MutationAA(pquery)
        aind=qaseq.sequence_position_to_aligned_index(aamu.position_)
        return aind
