import sys
import pysam
from uniprot import *
proteomefile="../data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
clinmissensebgzip="../validation/clinvar_annotation/clinvar_snpeff_missense_sorted.vcf.gz"
def main():
    up=UniprotProteomeParser(proteomefile)
    v=pysam.VariantFile(clinmissensebgzip)
    while True:
        try:
            r=next(v)
            ann=r.info['ANN']
            keyfound=[]
            for a in ann:
                a=a.split("|")
                #print(a[4])
                #uentry='-'
                try:
                    uentry=up.ensgid_entryname_[a[4]]
                except KeyError as ke:
                    uentry='-'
                kchange=a[4]+"_"+a[10]
                if not kchange in keyfound:
                    if a[1]=='missense_variant':
                        keyfound.append(kchange)
                        poss=a[12].split("/")
                        poss=poss[0]
                        try:
                            pline=uentry+"_"+poss+"\t"+r.info['GENEINFO']+"\t"+";".join(r.info['CLNDISDB'])+"\t"+";".join(r.info['CLNSIG'])+"\t"+a[2]+"\t"+a[3]+"\t"+a[4]+"\t"+a[6]+"\t"+a[10]+"\t"+a[12]
                            print(pline)
                        except KeyError as ee:
                            pline=uentry+"_"+poss+"\t"+r.info['GENEINFO']+"\t"+"-"+"\t"+";".join(r.info['CLNSIG'])+"\t"+a[2]+"\t"+a[3]+"\t"+a[4]+"\t"+a[6]+"\t"+a[10]+"\t"+a[12]
                            print(pline)
        except KeyError as te:
            pass
        except Exception as e:
            print(e)
            break
if __name__=="__main__":
    main()
