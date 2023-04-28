import shutil
import Bio.PDB
from Bio.SeqUtils import IUPACData
import sys
import os
from uniprot import *
from Bio import SeqIO

def is_struct_sequence_matching( accid, aseq, uniprot ):
    try:
        uobj=uniprot.get_protein_object_by_entryname(accid)
        useq=uobj.sequence_
        print(useq)
        return useq == aseq
    except KeyError:
        return False

def main():
    uniprot=UniprotProteomeParser(sys.argv[3])
    idir=sys.argv[1]#input dir
    odir=sys.argv[2]
    if not idir.endswith("/"):
        idir=idir+"/"
    filist=os.listdir(idir)
    matchcount=0
    p=Bio.PDB.PDBParser()
    for f in filist:
        #pdbid=sys.argv[1]
        #pdbfile=sys.argv[2]
        if f.endswith(".pdb"):
            #pdbid= f.strip().split("-")[1] #ex: AF-Q9Y6X2-F1-model_v1.pdb
            pdbid=""
            pdbfile=idir+f
            for record in SeqIO.parse(pdbfile, "pdb-seqres"):
                #print(record.dbxrefs)
                for i in record.dbxrefs:
                    if "_HUMAN" in i:
                        i=i.replace("UNP:","")
                        if i.find(" ")>-1:
                            si=i.split(" ")
                            for s in si:
                                if s.find("_HUMAN") > -1:
                                    pdbid=s
                                    break
                        else:
                            pdbid=i
            print("Processing: %s" %pdbid)
            #gmlfile="/media/param/DATA/code/python/lab/domain_mutations/dome2021/data/3d_structure/afolddb/dome3d/F1gml/"+pdbid+".edge"
            if len(pdbid) > 0:
                s=p.get_structure(pdbid, pdbfile)
                c1=s.get_chains()
                c=next(c1)
                res=c.get_residues()
                seqres=""
                while True:
                    try:
                        r=next(res)
                        b=r.resname[0]
                        for i in r.resname[1:]:
                            b+=str.lower(i)
                        aa=IUPACData.protein_letters_3to1[b]
                        seqres+=aa
                    except:
                        break
                if not odir.endswith("/"):
                    odir+="/"
                if is_struct_sequence_matching(pdbid, seqres, uniprot):
                    shutil.copy(pdbfile, odir)
                    matchcount+=1
                    print("Match found : %d %s" %( matchcount,pdbid))

if __name__=="__main__":
    main()
