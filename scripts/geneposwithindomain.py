from uniprot import *
import sys
import pandas as pd


def main():
    gd=pd.read_csv(sys.argv[2], header=0, sep="\t") #gene domain map
    muts=pd.read_csv(sys.argv[1], header=0, sep=",") #statsig format file
    u=UniprotProteomeParser(sys.argv[3])
    print("gene\tpos\tdomainID\tFDR\tcancer\tentryname")
    for i in muts.index:
        if len(gd[ gd["Gene"] == muts["gene"][i]]) > 0:
            a=gd[ gd["Gene"] == muts["gene"][i] ]
            gn=muts["gene"][i]
            ue=u.genename_entryname_[gn]
            if len(ue) > 5: # if entry name is assigned
                domfound=False
                for gi in a.index:
                    st=a["Start"][gi]
                    en=a["End"][gi]
                    if not pd.isnull(st) and not pd.isnull(en):
                        st=int(st)
                        en=int(en)
                        if int(muts["aapos"][i]) in range(st,en+1):
                            pfid=a["PfamID"][gi]
                            if not pd.isnull(pfid):
                                print(muts["gene"][i]+"\t"+str(muts["aapos"][i])+"\t"+pfid+"\t"+ str(muts["qvalue"][i]) +"\t"+ muts["tumors"][i]+"\t"+ ue )
                                domfound=True
                                break
                            #else:
                            #    print(muts["gene"][i]+"\t"+str(muts["aapos"][i])+"\t-\t"+ str(muts["qvalue"][i]) +"\t"+ muts["tumors"][i]+"\t"+ ue )
                if not domfound:
                    print(muts["gene"][i]+"\t"+str(muts["aapos"][i])+"\t-\t"+ str(muts["qvalue"][i]) +"\t"+ muts["tumors"][i]+"\t"+ ue )
                        #else:
                        #   print(muts["gene"][i]+"\t"+str(muts["aapos"][i])+"\t-\t"+ str(muts["qvalue"][i]) +"\t"+ muts["tumors"][i]+"\t"+ ue )
                    #else:
                    #    print(muts["gene"][i]+"\t"+str(muts["aapos"][i])+"\t-\t"+ str(muts["qvalue"][i]) +"\t"+ muts["tumors"][i]+"\t"+ ue )
        #else:
        #    print(muts["gene"][i]+"\t"+str(muts["aapos"][i])+"\t-\t"+ str(muts["qvalue"][i]) +"\t"+ muts["tumors"][i]+"\t"+ ue )

if __name__=="__main__":
    main()
