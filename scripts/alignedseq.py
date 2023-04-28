#!/usr/bin/python
import sys
import re
from seq import Sequence
from seq import ResidueError

class AlignedSequence(Sequence):
	def __init__(self, fn1=None, ssq1=None):
		super().__init__(fn1, ssq1)
		if not self.is_dna():
			self.is_protein()
	def __getitem__(self,i):
		return self.seq_[i]
	def is_gap(self,pos):
		return self.seq_[pos]=="-"
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
	def aligned_to_sequence_position(self, apos):
		spos=0
		i=0
		while i < apos:
			if self.seq_[i] != "-":
				spos=spos+1
			i=i+1
		return spos
	#Given actual sequence position returns the aligned position
	def sequence_to_aligned_position(self,spos):
		apos=0
		i=0
		while spos > 0:
			if not self.is_gap(i):
				apos = apos+1
				spos=spos-1
				i=i+1
			else:
				apos= apos+1
				i=i+1
		return apos
