from uniprot import *
import sys

#up=UniprotProteomeParser("/media/param/DATA/code/python/lab/domain_mutations/dome2021/data/uniprot/uniprot_hs_proteome_12_4_2021.txt")
#fi=open('../validation/cgc_uniq.tsv')
if len(sys.argv) == 3:
    up=UniprotProteomeParser(sys.argv[1])
    fi=open(sys.argv[2])
    for i in fi:
     try:
      print("%s\t%s" %(i.strip(),up.genename_entryname_[i.strip()]))
     except KeyError:
      pass
else:
    print("Usage: python3 history.py <uniprot_proteome.txt> <inputgenelistfile.tsv>")
