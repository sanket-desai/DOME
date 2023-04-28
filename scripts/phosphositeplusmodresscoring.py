#!/usr/bin/python3
#extracts functional sites from phosphositeplus file and mapped to individual proteins in form of a map
#columns: GENE    PROTEIN ACC_ID  HU_CHR_LOC      MOD_RSD SITE_GRP_ID     ORGANISM       MW_kD    DOMAIN  SITE_+/-7_AA    LT_LIT  MS_LIT  MS_CST  CST_CAT#
import sys
import os
import re
from uniprot import *
class PhosphositeplusRecordFormatError(Exception):
    def __init__(self, msg):
        self.message_=msg
        super().__init__(self.message_)

class PhosphositeplusRecord(object):
    def __init__(self, rec):
        srec=rec.rstrip().split('\t')
        self.gene_name_=srec[0]
        self.protein_=srec[1]
        self.accession_=srec[2]
        if self.accession_.find('-') != -1:
            self.accession_=self.accession_[:self.accession_.find('-')]
        self.chr_loc_=srec[3]
        self.mod_res_=srec[4]
        self.site_group_id_=srec[5]
        self.organism_=srec[6]
        self.mol_weight_=srec[7]
        self.domain_=srec[8]
        self.site_sequence_=srec[9]
    def __str__(self):
        return "%s\t%s\t%d\t%s\t%s\t%s" %(self.gene_name_, self.accession_, self.get_position(), self.get_residue() , self.get_modification_type(), self.mod_res_)
    def get_modification_type(self):
        mtype=self.mod_res_.split('-')
        if len(mtype)!=2:
            raise PhosphositeplusRecordFormatError("Modification type cannot be determined : %s" %(self.mod_res_))
        return mtype[1]
    def get_residue(self):
        mtype=self.mod_res_.split('-')
        if len(mtype)!=2:
            raise PhosphositeplusRecordFormatError("Modification type cannot be determined : %s" %(self.mod_res_))
        #return mtype[0][0]
        temp=re.findall(r'\D+', mtype[0])
        if len(temp) != 1:
            raise PhosphositeplusRecordFormatError("Protein residue not found in the modification information: %s" %(self.mod_res_))
        return temp[0]
    def get_position(self):
        mtype=self.mod_res_.split('-')
        if len(mtype)!=2:
            raise PhosphositeplusRecordFormatError("Modification type cannot be determined : %s" %(self.mod_res_))
        temp=re.findall(r'\d+', mtype[0])
        if len(temp) != 1:
            raise PhosphositeplusRecordFormatError("Protein position not found in the modification information: %s" %(self.mod_res_))
        return int(temp[0])
    def is_human(self):
        return self.organism_=='human'

class PhosphositeplusParser(object):
    def __init__(self, ifile):
        self.filename_=ifile
        fi=open(ifile)
        self.header_=fi.readline().strip().split('\t')
        self.phosphositeplus_functional_site_map_={} # uniprot_accession -> [phosphositerecords]
        for l in fi:
            psprec=PhosphositeplusRecord(l)
            if not psprec.accession_ in self.phosphositeplus_functional_site_map_:
                self.phosphositeplus_functional_site_map_[psprec.accession_]=[ psprec ]
            else:
                precs=self.phosphositeplus_functional_site_map_[psprec.accession_]
                precs.append(psprec)
                self.phosphositeplus_functional_site_map_[psprec.accession_]=precs
    def get_number_of_modifications(self, uacc): #query based on accession
        precs=self.phosphositeplus_functional_site_map_[uacc]
        return len(precs)
    def get_phosphositeplus_sites(self, uacc):
        return self.phosphositeplus_functional_site_map_[uacc]
    def get_uniprot_accessions(self):
        return self.phosphositeplus_functional_site_map_.keys()

def main():
    uniprot= UniprotProteomeParser("/media/param/DATA/code/python/lab/domain_mutations/dome2021/data/uniprot/uniprot_hs_proteome_12_4_2021.txt")
    fn='/media/param/DATA/code/python/lab/domain_mutations/dome2021/data/phosphositeplus/phosphositeplus_all_modres_human.tsv'
    psp=PhosphositeplusParser(fn)
    for u in psp.get_uniprot_accessions():
        try:
            uentry=uniprot.accession_entryname_[u]
            x=psp.get_phosphositeplus_sites(u)
            xx=x[1:3]
            for i in xx:
                print("%s\t%s" %(uentry, str(i)))
        except KeyError:
            pass
    #x=psp.get_phosphositeplus_sites('P61981')
    #xx=x[1:3]
    #for i in xx:
    #    print(i.get_position())
    #    print(i.get_modification_type())
    #    print(str(i))
if __name__=='__main__':
    main()
