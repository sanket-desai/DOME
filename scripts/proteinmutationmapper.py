#!/usr/bin/python3
#Construct a MutationMapper, input MAF files
#Structure {uniprot->[proteinmutation1, proteinmutation2, .. ]}
#Alternatively
#uniprotid-> Mutations [[ gene_name, original_sequence, mutation_array[ [],[patientid1_change, patient_id2_change, ...],[] ] ]]
import sys
import re

#MISCOUNTER=0

class FormatError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)
#Takes protein mutated positions and generates a 21mer peptides and creates a multiseq object
class Mutation(object):
	def __init__(self,res=None,pos=None,ch=None):
		self.residue_=res
		self.position_=int(pos)
		self.change_=ch
	def __str__(self):
		return str("Pos:%d,Ref:%s,Alter:%s" %(self.position_,self.residue_,self.change_))
class MutationAA(Mutation):
	def __init__(self,word):
		sf=re.findall('\d+',word)
		if len(sf)>1:
			raise FormatError("Mutation %s contains more than one position, please check." %(word))
		elif len(sf)<1:
			raise FormatError("Mutation %s does not contain position information. Please check." %(word))
		else:
			super(MutationAA, self).__init__(word[word.find(".")+1:word.find(sf[0])], int(sf[0]), word[word.find(sf[0])+len(sf[0]):] )

class Mutations(object):
	def __init__(self, gname, sequence, mutarray):
		self.gene_name_=gname
		self.sequence_=sequence
		self.mutation_array_=mutarray #[sample1_changeaa, sample1_changeaa2, ]
		#print("%s %s" %(self.gene_name_, self.sequence_))
		#print(len(mutarray))
		#print(mutarray)
	def __len__(self):
		return len(self.mutation_array_)
	def __getitem__(self,index):
		return self.mutation_array_[index]
	def __setitem__(self, index, samchg):
		ma=self.mutation_array_[index]
		if type(ma)==str:
			ma=[samchg]
		elif type(ma)==list:
			ma.append(samchg)
		self.mutation_array_[index]=ma
	def __str__(self):
		mtc=self.mutations_to_counts()
		muc=str(mtc).replace(" ","").replace("[","").replace("]","")
		return "%s %s %d" %(self.gene_name_,muc, sum(mtc))
	def mutations_to_counts(self): #returns a position wise count of mutations for a particular protein
		counts=[0]*len(self.sequence_)
		i=0
		for m in self.mutation_array_:
			if type(m)==list:
				counts[i]=len(m)
			i=i+1	
		return counts
	def total_mutations(self):
		return sum(self.mutations_to_counts_())
	def samples_mutated(self):
		samples=[]
		for m in self.mutation_array_:
			if m != '':
				for sa in m:
					samples.append(sa.split("_")[0])
		return samples
	def samples_mutated_by_position(self, pos):
		mp=self.mutation_array_[pos]
		samples=[]
		if m != '':
			for sa in m:
				samples.append(sa.split("_")[0])
		return samples

#Give a specific file format to be given as input, may be tab delimited like in nature pancan supplement
#File format expected and the program is quite strict about:
#uniprotid, gene_name, patient_id, mutation (ex; P238T meaning at position 238 P->T)
class ProteinMutationMapper(object):
	def __init__(self, infile=None, seqfile=None): #Seqfile is an tab file with id sequence
		if not infile:
			print("Input file in specified format is not provided! Terminating the program.")
			sys.exit(0)
		if not seqfile:
			print("Sequence file with uniprot id and sequence not provided. Cannot proceed further.")
			sys.exit(0)
			#check it exists and has the exact number of columns; also the headers
			#UniprotId GeneName SampleId Mutation
			#seqfile should have format:
			#uniprotid sequence
		#{ uniprot->Mutations[[ ],[ ], ... ] }
		self.mutation_map_={}
		uniseq={}
		fhs=open(seqfile)
		for s in fhs:
			ss=s.strip().split("\t")
			if len(ss)!=2:
				raise FormatError("Problem in Uniprot sequence file!!")
			else:
				uniseq[ss[0].strip()]=ss[1].strip()
		fhs.close()
		fhm=open(infile)
		for li in fhm:
			sli=li.strip().split("\t")
			if len(sli)!=4:
				raise FormatError("Mutation file is not consistent with the defined format. Please check the documentation!!")
			else:
				if sli[0].strip() in self.mutation_map_: #sli[0] ~ uniprotid
					#since the mutation array has been already initialized
					try: 
						maa=MutationAA(sli[3])
						samchange=sli[2]+"_"+maa.change_
						try:
							#global MISCOUNTER
							#MISCOUNTER=MISCOUNTER+1
							if maa.residue_==self.mutation_map_[sli[0].strip()].sequence_[maa.position_-1] and not samchange in self.mutation_map_[sli[0].strip()][maa.position_-1]:
								self.mutation_map_[sli[0].strip()][maa.position_-1]=samchange
						except IndexError:
							pass #Write something for this, will ya
					except FormatError:
						pass
				else:
					#Mutations(genename, sequence, mutationarray) #sli[1] is the gene name
					try:
						maa=MutationAA(sli[3])
						samchange=sli[2]+"_"+maa.change_
						try:
							#global MISCOUNTER
							#MISCOUNTER=MISCOUNTER+1
							mutsobj=Mutations(sli[1], uniseq[sli[0].strip()] ,['']*len(uniseq[sli[0].strip()])) #uniseq[sli[0]] is expected to return sequence
							if maa.residue_==mutsobj.sequence_[maa.position_-1] and not samchange in mutsobj[maa.position_-1]:
								mutsobj[maa.position_-1]=samchange
							self.mutation_map_[sli[0].strip()]=mutsobj
						except KeyError:
							print("%s not found in the uniprot proteome." %(sli[0].strip()))
						except IndexError:
							print("Index error encountered. %d %d" %(maa.position_-1, len(mutsobj.sequence_)))
							pass
					except FormatError:
						pass
		fhm.close()
	def number_of_mutations(self, uid, pos):
		return len(self.mutation_map_[uid][pos-1])
	def __getitem__(self, uid):
		return self.mutation_map_[uid]
	def uniprot_to_gene(self, uid):
		return self.mutation_map_[uid].gene_name_
	def gene_to_uniprot(self, genm):
		uid=None
		for ud in self.mutation_map_:
			if self.mutation_map_[ud].gene_name_==genm:
				uid=str(ud)
				break
		return uid
	def get_mutations_by_gene(self, genenm):
		if self.gene_to_uniprot(genenm):
			return self.mutation_map_[self.gene_to_uniprot(genenm)]
		else:
			print("Gene not found!!")
			return None
	def get_samples_mutated(self, uid):
		return self.mutation_map_[uid].samples_mutated()
	def get_samples_mutated_by_gene(self, gene):
		return self.get_mutations_by_gene(gene).samples_mutated()
			
#argv[1]=mutation file
#argv[2]=sequence file
#argv[3]=mutation list for each gene, output filename 
#def main():
#	if len(sys.argv)!=4:
#		print("Bad input. Program terminated!")
#		sys.exit(0)
#	else:
#		pmm=ProteinMutationMapper(sys.argv[1], sys.argv[2])
#		outf=open(sys.argv[3],"w")
#		outf.write("Protein MutationCounts TotalMutations\n")
#		for uid in pmm.mutation_map_:
#			outf.write("%s\n" %(str(pmm.mutation_map_[uid])))
#			print("%s processed." %(str(uid)))	
#		print("Processed Missense mutations: %d" %(MISCOUNTER))

#if __name__=="__main__":
#	main()
