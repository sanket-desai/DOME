import sys
from Bio import SeqIO
import os, fnmatch
from uniprot import *


def int_to_str_list(li):
    strlist=[]
    for i in li:
        strlist.append(str(i))
    return strlist

class FormatException(Exception):
    def __init__(self, e):
        super().__init__(e)

class IDLineRecord(object): #technically id record line
    def __init__(self, il):
        if not il.startswith("#=ID"):
            raise FormatException("Error in the record line %s" %(il.strip()))
        sil=il.strip().split("\t")
        self.domainname1_=sil[1]
        self.domainname2_=sil[2]
        domid1=sil[3].replace("(","")
        domid1=domid1[:domid1.find(".")]
        self.domainid1_=domid1
        domid1=sil[4].replace(")","")
        domid1=domid1[:domid1.find(".")]
        self.domainid2_=domid1
    def get_id_key(self):
        return self.domainname1_+":"+self.domainname2_

class PDB3DLineRecord(object): #technically line starting with 3D
    def __init__(self, il):
        if not il.startswith("#=3D"):
            raise FormatException("Error in the record line %s" %(il.strip()))
        sil=il.strip().split("\t")
        self.pdbid_=sil[1]
        siltemp=sil[2].split(":")
        self.chain1_=siltemp[0]
        temp=siltemp[1].split("-")
        self.chain1domainend_=0
        self.chain1domainend_=0
        self.chain2domainend_=0
        self.chain2domainstart_=0
        try:
            if len(temp)==2:
                self.chain1domainstart_=temp[0]
                self.chain1domainend_=temp[1]
            if len(temp)==3 and temp[0]=="":
                self.chain1domainstart_=0-int(temp[1])
                self.chain1domainend_=int(temp[2])
            if len(temp)==3 and temp[1]=="":
                self.chain1domainstart_=int(temp[0])
                self.chain1domainend_=0-int(temp[2])
            if len(temp)==4:
                self.chain1domainstart_=0-int(temp[1])
                self.chain1domainend_=0-int(temp[3])
        except ValueError:
            pass
        siltemp=sil[3].split(":")
        self.chain2_=siltemp[0]
        try:
            temp=siltemp[1].split("-")
            if len(temp)==2:
                self.chain2domainstart_=temp[0]
                self.chain2domainend_=temp[1]
            if len(temp)==3 and temp[0]=="":
                self.chain2domainstart_=0-int(temp[1])
                self.chain3domainend_=int(temp[2])
            if len(temp)==3 and temp[1]=="":
                self.chain2domainstart_=int(temp[0])
                self.chain2domainend_=0-int(temp[2])
            if len(temp)==4:
                self.chain2domainstart_=0-int(temp[1])
                self.chain2domainend_=0-int(temp[3])
        except ValueError:
            pass
        self.score_=float(sil[4])
        self.zscore_=float(sil[5])
        try:
            self.topology_=sil[6]
        except IndexError:
            self.topology_=""

class InteractionRecord(object): #individual 3d record residue pair entries
    def __init__(self, il):
        sil=il.strip().split("\t")
        self.residue1_=sil[0]
        self.residue2_=sil[1]
        self.position1_=int(sil[2])
        self.position2_=int(sil[3])
        self.contacttype_=sil[4]

class Interface(object): #individual PDB interacting residues / #=3D block
    def __init__(self, pdbinteractionblockvec):
        pdb3dlinerecord=""
        self.interactionresidues_=[]
        for i in pdbinteractionblockvec:
            if i.startswith("#=3D"):
                pdb3dlinerecord=i
                #self.pdb3dlinerecord_=PDB3DLineRecord(i)
            else:
                try:
                    #print(i)
                    self.interactionresidues_.append(InteractionRecord(i))
                except ValueError:
                    pass
        self.pdb3dlinerecord_=PDB3DLineRecord(pdb3dlinerecord)
    def get_number_of_interacting_residues(self):
        return len(self.interactionresidues_)
    def get_pdb_id(self):
        return self.pdb3dlinerecord_.pdbid_
    def get_chain1(self):
        return self.pdb3dlinerecord_.chain1_
    def get_chain2(self):
        return self.pdb3dlinerecord_.chain2_
    def get_chain1_residues(self):
        residues=[]
        for ir in self.interactionresidues_:
            residues.append(ir.residue1_)
        return residues
    def get_chain2_residues(self):
        residues=[]
        for ir in self.interactionresidues_:
            residues.append(ir.residue2_)
        return residues
    def get_chain1_positions(self):
        positions=[]
        for ir in self.interactionresidues_:
            positions.append(ir.position1_)
        return positions
    def get_chain2_positions(self):
        positions=[]
        for ir in self.interactionresidues_:
            positions.append(ir.position2_)
        return positions

class DomainInteraction(object): #One interaction block / list of interfaces
    def __init__(self, interactionblockvec):
        self.did_="" #domain interaction id / key
        self.interactionidobj_=""
        self.interfaceobjects_=[] #list of interface objects
        pdbinterfaceblock=[]
        for i in interactionblockvec:
            if not i.startswith("//"):
                if i.startswith("#=ID"):
                    self.interactionidobj_=IDLineRecord(i)
                elif i.startswith("#=3D"):
                    if len(pdbinterfaceblock) > 1:
                        #print(pdbinterfaceblock)
                        self.interfaceobjects_.append(Interface(pdbinterfaceblock))
                        pdbinterfaceblock=[]
                    pdbinterfaceblock.append(i)
                else:
                    pdbinterfaceblock.append(i)
            else:
                break
        self.interfaceobjects_.append(Interface(pdbinterfaceblock))
        self.did_=self.interactionidobj_.get_id_key()
    def get_id_key(self):
        return self.did_

class DomainInteractionDb(object):
    def __init__(self, didfile):
        self.domaininteractions_=[]
        intblock=[]
        df=open(didfile)
        for d in df:
            d=d.strip()
            if d.startswith("#=ID"):
                if len(intblock) >1 :
                    self.domaininteractions_.append(DomainInteraction(intblock))
                    #print(len(self.domaininteractions_))
                    intblock=[]
            intblock.append(d)
        df.close()
        self.domaininteractions_.append(DomainInteraction(intblock))

#PDB directory
#contains pdb_chain->uniprot entry map
#atom section sequence mapping with the uniprot database
class PDBDbHandler(object):
    def __init__(self, pdbdir):
        self.pdb_entry_map_={} #entry_name -> [ pdbid_chain, ]
        self.pdb_dir_=pdbdir
        self.pdb_files_=[]
        if not pdbdir.endswith("/"):
            pdbdir=pdbdir+"/"
        pdbs=fnmatch.filter(os.listdir(pdbdir), '*.pdb')
        for p in pdbs:
            self.pdb_files_.append(pdbdir+p)
            sp=SeqIO.parse(pdbdir+p, "pdb-seqres")
            for rec in sp:
                #check if dbxref is uniprot
                try:
                    if rec.dbxrefs[3].startswith("UNP"):
                        pid=rec.dbxrefs[0].replace("PDB:","").lower()
                        pchain=rec.annotations['chain']
                        pdbidchain=pid+"_"+pchain
                        entryname=rec.dbxrefs[3].replace("UNP:","")
                        if not pdbidchain in self.pdb_entry_map_:
                            self.pdb_entry_map_[pdbidchain]=entryname
                except IndexError:
                    pass
    def has_pdb_file(self, pdbid):
        return self.pdb_dir_+pdbid+".pdb" in self.pdb_files_
    def get_entry_name(self, pdbid, chainid):
        return self.pdb_entry_map_[pdbid+"_"+chainid]

#didfile="3did_flat.txt"
#did3dfh=open(didfile)
def main():
    #outfile=sys.argv[1]
    outfile="/mnt/d/code/python/lab/domain_mutations/dome2021/data/ppi/interaction/human_interfaces_3did_uniprotannotated.tsv"
    didfile="/mnt/d/code/python/lab/domain_mutations/dome2021/data/ppi/interaction/3did_flat.txt"
    unip="/mnt/d/code/python/lab/domain_mutations/dome2021/data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
    did=DomainInteractionDb(didfile)
    print("Domain interaction Db created!")
    humproteome=UniprotProteomeParser(unip)
    print("Proteome parsed!")
    ofile=open(outfile,'w')
    p=PDBDbHandler("/mnt/d/code/python/lab/domain_mutations/dome2021/data/3d_structure/pdb_download/human/didpdb_matching/")
    print("PDBDbHandler object created!")
    counter=0
    for d in did.domaininteractions_:
        for i in d.interfaceobjects_:
            counter+=1
            if p.has_pdb_file(i.get_pdb_id()):
                try:
                    p1matching=0 #number of matching residues in protein 1
                    p2matching=0
                    en1=p.get_entry_name(i.get_pdb_id(),i.get_chain1())
                    en2=p.get_entry_name(i.get_pdb_id(),i.get_chain2())
                    #if humproteome.is_entryname_within_uniprot_proteome(en1) and humproteome.is_entryname_within_uniprot_proteome(en2):
                    #    pr1=humproteome.get_protein_object_by_entryname(en1)
                    #    pr2=humproteome.get_protein_object_by_entryname(en2)
                    for ind1 in range(0, len(i.get_chain1_positions())):
                        if humproteome.sequence_residue_check_by_entry( en1, i.get_chain1_positions()[ind1], i.get_chain1_residues()[ind1]):
                            p1matching+=1
                    for ind2 in range(0, len(i.get_chain2_positions())):
                        if humproteome.sequence_residue_check_by_entry( en2, i.get_chain2_positions()[ind2], i.get_chain2_residues()[ind2]):
                            p2matching+=1
                    if (p1matching/len(i.get_chain1_positions())) == 1 and (p2matching/len(i.get_chain2_positions())) == 1: #80 percent residue match
                        #print(d.get_id_key()+"\t"+en1+"\t"+en2+"\t"+i.get_pdb_id()+"\t"+i.get_chain1()+"\t"+i.get_chain2()+"\t"+",".join( int_to_str_list( i.get_chain1_positions() ) )+"\t"+",".join( int_to_str_list( i.get_chain2_positions()) )+"\t"+",".join(i.get_chain1_residues()) +"\t"+",".join(i.get_chain2_residues()) )
                        print("%d %s" %(counter, d.get_id_key()))
                        sstr=d.get_id_key()+"\t"+en1+"\t"+en2+"\t"+i.get_pdb_id()+"\t"+i.get_chain1()+"\t"+i.get_chain2()+"\t"+",".join( int_to_str_list( i.get_chain1_positions() ) )+"\t"+",".join( int_to_str_list( i.get_chain2_positions()) )+"\t"+",".join(i.get_chain1_residues()) +"\t"+",".join(i.get_chain2_residues())
                        ofile.write("%s" %(sstr))
                except KeyError as ke:
                    print(ke)
                    pass
                    #print("Entry not found %s" %(i.get_pdb_id()) )
                    #print(d.get_id_key()+"\t"+i.get_pdb_id()+"\t"+i.get_chain1()+"\t"+i.get_chain2()+"\t"+",".join( int_to_str_list( i.get_chain1_positions() ) )+"\t"+",".join( int_to_str_list( i.get_chain2_positions()) )+"\t"+",".join(i.get_chain1_residues()) +"\t"+",".join(i.get_chain2_residues()) )
    ofile.close()

if __name__=="__main__":
    main()
