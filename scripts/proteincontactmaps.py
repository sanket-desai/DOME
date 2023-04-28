import matplotlib.pyplot as plt
import networkx as nx
import Bio.PDB
from Bio.SeqUtils import IUPACData
import numpy
import sys
import os
from Bio import SeqIO

def res_to_single_letter(tlet):
    slet=tlet[0]
    for i in tlet[1:]:
        slet += str.lower(i)
    return IUPACData.protein_letters_3to1[slet]

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

def tuple_in_list(tup, li):
    ispresent=False
    for t in li:
        if (tup[0] in t) and (tup[1] in t):
            ispresent=True
            break
    return ispresent

class ProteinContactMap(nx.Graph):
    def __init__(self, pname , edgefile):
        self.protein_=pname
        nx.Graph.__init__(nx.read_weighted_edgelist(edgefile))



def main():
    idir=sys.argv[1]#input dir
    if not idir.endswith("/"):
        idir=idir+"/"
    filist=os.listdir(idir)
    odir=sys.argv[2]
    if not odir.endswith("/"):
        odir+="/"
    p=Bio.PDB.PDBParser()
    for f in filist:
        #pdbid=sys.argv[1]
        #pdbfile=sys.argv[2]
        if f.endswith(".pdb"):
            #pdbid= f.strip().split("-")[1] #ex: AF-Q9Y6X2-F1-model_v1.pdb
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
            gmlfile=odir+pdbid+".edge"
            s=p.get_structure(pdbid, pdbfile)
            c1=s.get_chains()
            c=next(c1)
            res=c.get_residues()
            aalist=[]
            while True:
                try:
                    r=next(res)
                    b=r.resname[0]
                    for i in r.resname[1:]:
                        b+=str.lower(i)
                    aa=IUPACData.protein_letters_3to1[b]
                    aalist.append(aa)
                except:
                    break
            g=nx.Graph()
            #print("Length of protein: %d" %len(aalist))
            mat=calc_dist_matrix(c,c)
            print("Chain length: %d" %(len(c)))
            nrow, ncol=mat.shape
            print("Number of rows and columns: %d %d" %(nrow, ncol))
            #sys.exit(0)
            pairlist=[] #pair in
            for row in range(len(aalist)):
                for col in range(row):
                    if row!=col:
                        print("Row col %d %d" %(row, col))
                        #pairlist.append(t1)
                        #if not str(row)+"-"+str(col) in pairlist  or not str(col)+"-"+str(row) in pairlist:
                        if mat[row, col] < 10:
                            #print("%s\t%d\t%s\t%d\t%s\t%f" %( pdbid, row+1, aalist[row], col+1, aalist[col] , mat[row,col]) )
                            g.add_weighted_edges_from([(row+1,col+1, float("{:.3f}".format(mat[row, col])) )])
            nx.write_weighted_edgelist(g, gmlfile)
if __name__=="__main__":
    main()
