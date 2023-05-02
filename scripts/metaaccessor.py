#!/usr/bin/python3

#Access all the databases including long and short range 3D interaction with residues within protein to prioritize mutations
#Basic aim of this script is to list down criteria using which we may call a position / mutation important
#Algorithm:
'''
Position within domain
Position geometrically close to domain
Position is non-neutral (CADD)
Position is itself a
- known statistically significant or analogous mutation position (given a database of statistically significant mutations)
- known phosphositeplus site
- known uniprot functional site
- mutagenesis experiment known (uniprot)

Geometric closeness to
- known statistically significant mutations
- known phosphositeplus (PTM) site
- known uniprot functional site

Database annotations
- Clinvar
- COSMIC
'''
import os
import sys
import pysam
import pandas as pd
import numpy as np
from aacids import *
import networkx as nx
from domainmultialignedseq import *
#Global paths
class RefPath(object):
    def __init__(self, path):
        #domerefdata_="../data/domerefdata/final/"
        if not path.endswith("/"):
            path=path+"/"
        self.refpath_=path
        self.analogous_=self.refpath_+"analogous_positions.tsv"
        self.statsig_=self.refpath_+"statsig_tcga_mutations.tsv"
        self.cadddb_=self.refpath_+"cadd_scores_missense.tsv.gz"
        self.cosmicdb_=self.refpath_+"cosmic_missense_counts.tsv.gz"
        self.cosmicgenomewidedb_=self.refpath_+"CosmicMutantExport_missense_genomewide_cor_entry_pchangecount_ord.tsv.gz"
        self.clinvardb_=self.refpath_+"clinvar_snpeff_missense.tsv.gz"
        self.pspdb_=self.refpath_+"phosphositeplus_human.tsv.gz"
        self.umutagendb_=self.refpath_+"uniprot_mutagenesis.tsv.gz"
        self.upfamdb_=self.refpath_+"uniprot_pfam.tsv.gz"
        self.usitesdb_=self.refpath_+"uniprot_sites.tsv.gz"
        self.edgedir_=self.refpath_+"pdbedgelist/"
        self.hotspot3ds_=self.refpath_+"3d_hotspots_residues_entry.tsv"
        self.didproteininterfaces_=self.refpath_+"protein_interface_positions.tsv"
        self.unidentrymap_=self.refpath_+"unid_entry_map.tsv"
        self.domalignpath_=self.refpath_+"mafft_outaln/"
        self.domainindexproteomeentropydb_=self.refpath_+ "proteome_domain_entropy_cosmicgw.tsv.gz"
        self.proteinposdomainindexdb_=self.refpath_+"protein_position_domain_index.tsv.gz"

class Query(object):
    def __init__(self, entry, refpos, altaa=""): #changed AA may not be given
        self.entryname_=entry
        self.position_=int(refpos)
        self.altaa_=altaa
    def __str__(self):
        return self.entryname_ + "|"+ str(self.position_)+"|"+ self.altaa_
    def tabix_query(self):
        return self.entryname_+":"+self.position_+"-"+self.position_

class ProteinInterfaces(object):
    def __init__(self, refpathobj):
        fi=open(refpathobj.didproteininterfaces_)
        self.proteininterfaces_={}
        fii=open(refpathobj.unidentrymap_)
        self.uniprotid_entryname_map_={}
        for i in fii:
            si=i.strip().split("\t")
            if len(si)==2:
                self.uniprotid_entryname_map_[si[0]]=si[1]
        for i in fi:
            si=i.strip().split("\t")
            try:
                self.proteininterfaces_[self.get_uniprotid_by_entryname(si[0])]=si[1].strip().split(",")
            except KeyError:
                pass
        #print(self.proteininterfaces_)
        #self.proteininterfaces_=pd.read_csv(refpathobj.didproteininterfaces_, header=0, sep="\t")
        #piis=pd.read_csv(refpathobj.didproteininterfaces_, header=0, sep="\t")

    def get_uniprotid_by_entryname(self, uid):
        return self.uniprotid_entryname_map_[uid]
    def get_interface_positions_by_entryname(self, entry):
        ipos=[]
        if entry in self.proteininterfaces_:
            pos=self.proteininterfaces_[entry]
            spos=[]
            for s in pos:
                ipos.append(int(s))
        return ipos

class Hotspot3dPositions(object):
    def __init__(self, refpathobj):
        self.hotspot3dframe_=pd.read_csv(refpathobj.hotspot3ds_, header=0, sep="\t")
    def get_hotspot3d_position_records(self, q):
        selmutframe=self.hotspot3dframe_[self.hotspot3dframe_.entry==q.entryname_]
        return selmutframe[selmutframe.pos == q.position_]
    def get_hotspot3d_positions_for_protein(self, q):
        selmutframe=self.hotspot3dframe_[self.hotspot3dframe_.entry==q.entryname_]
        return list(selmutframe.pos)
    def is_hotspot3d_position(self, q):
        return self.get_hotspot3d_position_records(q).shape[1] > 0
    def protein_has_hotspot3d_positions(self, q):
        return len(self.get_hotspot3d_positions_for_protein(q)) > 0

class SignificantMutations(object):
    def __init__(self, refpathobj):
        self.sigmutframe_=pd.read_csv( refpathobj.statsig_ , header=0, sep="\t")
    def get_significant_position_records(self, q): #returns a pandas dataframe object
        selmutframe= self.sigmutframe_[self.sigmutframe_.entryname==q.entryname_]
        return selmutframe[ selmutframe.pos == q.position_ ]
    def get_significant_position_records_for_protein(self, q):
        return self.sigmutframe_[self.sigmutframe_.entryname==q.entryname_]
    def get_significant_positions_for_protein(self, q):
        return list(self.get_significant_position_records_for_protein(q).pos)
    def is_position_significant(self, q): #check if the query is singificant (q being the query object)
        issignificant=False
        if self.get_significant_position_records(q).shape[0]>0:
            issignificant=True
        return issignificant
    def get_significance_scores(self,q): #assumes mutation is stat sig
        sigscores= []
        smr=self.get_significant_position_records(q)
        if smr.shape[0]>0:
            #sigscores=list( -np.log(smr.FDR) )
            #in the newest version fdr is replaced by H, R (hotspot or resistant mutation)
            sigscores=list(smr.FDR)
        return sigscores
    def get_tumour_types(self, q): #assumes mutation is stat sig
        ttypes= []
        smr=self.get_significant_position_records(q)
        if smr.shape[0]>0:
            ttypes=list(smr.cancer)
        return ttypes
    def get_tumour_types_for_protein(self, q):
        return list(self.get_significant_position_records_for_protein(q).cancer)

class AnalogousMutations(object):
    def __init__(self, refpathobj):
        self.anamutframe_=pd.read_csv( refpathobj.analogous_ , header=0, sep="\t")
    def get_analogous_positions_for_protein(self, entry): #given protein uniprot entry name
        fr=self.anamutframe_[self.anamutframe_.entry==entry]
        return list(fr.pos)
    def get_analogous_position_records(self, q):
        selmutframe= self.anamutframe_[self.anamutframe_.entry==q.entryname_]
        return selmutframe[ selmutframe.pos == q.position_ ]
    def is_position_analogous(self, q):
        isanalogous=False
        if self.get_analogous_position_records(q).shape[0]>0:
            isanalogous=True
        return isanalogous
    def get_analogous_position_residue(self, q):
        anaposres=''
        amr=self.get_analogous_position_records(q)
        if amr.shape[0]>0:
            anaposres=list(amr.residue)[0]
        return anaposres
    def get_analogous_significance_score(self, q):
        sigscore=0
        amr=self.get_analogous_position_records(q)
        if amr.shape[0]>0:
            sigscore=list(amr.statscore)[0]
        return sigscore
    def get_analogous_position_domain_id(self, q):
        domid=''
        amr=self.get_analogous_position_records(q)
        if amr.shape[0]>0:
            domid=list(amr.domainid)[0]
        return domid
    def get_analogous_position_domain_name(self, q):
        domid=''
        amr=self.get_analogous_position_records(q)
        if amr.shape[0]>0:
            domid=list(amr.domainname)[0]
        return domid
    '''
    def get_analosous_position_cosmic_mutation_count(self, q):
        cosmcount=0
        amr=self.get_analogous_position_records(q)
        if amr.shape[0]>0:
            cosmcount=list(amr.mcount)[0]
        return cosmcount
    def is_position_mutation_neutral(self, q):
        diffmutcadd=0
        amr=self.get_analogous_position_records(q)
        if amr.shape[0]>0:
            diffmutcadd=list(amr.diff_aa_caddphred)[0]
        return diffmutcadd <= 10
    '''
    def get_significant_protein_to_analogous_position(self, q):
        amr=self.get_analogous_position_records(q)
        sentry=[]
        if amr.shape[0]>0:
            sentry=list(amr.statsigentry)
        return sentry
    def get_significant_protein_position_to_analogous_position(self, q):
        amr=self.get_analogous_position_records(q)
        spos=[0]
        if amr.shape[0]>0:
            spos=list(amr.statpos)
        return spos
    def get_significant_protein_residue_to_analogous_position(self, q):
        amr=self.get_analogous_position_records(q)
        sentry=[]
        if amr.shape[0]>0:
            sentry=list(amr.statsigres)
        return sentry

class UniprotIDEntryNameMap(object):
    def __init__(self, refpathobj):
        fi=open(refpathobj.unidentrymap_)
        self.uniprotid_entryname_map_={}
        for i in fi:
            si=i.strip().split("\t")
            if len(si)==2:
                self.uniprotid_entryname_map_[si[0]]=si[1]
        fi.close()
    def get_entryname_by_uniprotid(self, uid):
        return self.uniprotid_entryname_map_[uid]


class UniprotPfamDbRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[0]
        self.domstart_=int(srec[1])
        self.domend_=int(srec[2])
        self.accession_=srec[3]
        self.domainid_=srec[4]
        self.domainname_=srec[5]
class UniprotPfamDb(object):
    def __init__(self, refpathobj):
        self.upfamtabixfile_=pysam.TabixFile(refpathobj.upfamdb_)
    def get_pfam_domain_records(self, q):
        tabrec=[] #string list
        try:
            tbi=self.upfamtabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                ii=UniprotPfamDbRecord(i)
                if q.position_ in range(ii.domstart_, ii.domend_):
                    tabrec.append(ii)
        except ValueError:
            pass
        return tabrec
    def get_pfam_domain_records_for_protein_entry(self, e):
        tabrec=[] #uniprotpfamdbrecord list
        try:
            tbi=self.upfamtabixfile_.fetch(e)
            for i in tbi:
                ii=UniprotPfamDbRecord(i)
                tabrec.append(ii)
        except ValueError:
            pass
        return tabrec
    def get_pfam_domain_records_for_protein(self, q):
        tabrec=[] #uniprotpfamdbrecord list
        try:
            tbi=self.upfamtabixfile_.fetch(q.entryname_)
            for i in tbi:
                ii=UniprotPfamDbRecord(i)
                tabrec.append(ii)
        except ValueError:
            pass
        return tabrec
    def is_within_domain(self, q):
        return len(self.get_pfam_domain_records(q))>0
    def get_pfam_domain_name(self, q):
        recs=self.get_pfam_domain_records(q)
        dnames=[]
        dname=""
        for r in recs:
            dnames.append(r.domainname_)
        if len(dnames)>0:
            dnames=list(np.unique(np.array(dnames)))
            dname=";".join(dnames)
            #only a hack to prevent joining two names
            dname=dnames[0]
        return dname
    def get_pfam_ids(self, q):
        recs=self.get_pfam_domain_records(q)
        dnames=[]
        for r in recs:
            dnames.append(r.domainid_)
        return dnames
    def get_pfam_domain_start(self, q):
        recs=self.get_pfam_domain_records(q)
        dnames=[]
        for r in recs:
            dnames.append(r.domainstart_)
        return dnames
    def get_pfam_domain_start(self, q):
        recs=self.get_pfam_domain_records(q)
        dnames=[]
        for r in recs:
            dnames.append(r.domainend_)
        return dnames
    def get_pfam_domain_accessions(self, q):
        recs=self.get_pfam_domain_records(q)
        dnames=[]
        for r in recs:
            dnames.append(r.accession_)
        return dnames

#Upon canonical transcript parsing from the
#change code in CADD parser
#ABL1_HUMAN      315     T       S       3.704133        25.4    ENST00000318560.6       9       130872896 C       G
class CADDScoreRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        '''
        self.entryname_=srec[0]
        self.refaa_=srec[1]
        self.altaa_=srec[2]
        self.geneid_=srec[3]
        self.transcriptid_=srec[4]
        self.genename_=srec[5]
        self.posaa_=int(srec[6])
        self.caddrawscore_=float(srec[7])
        self.caddphredscore_=float(srec[8])
        '''
        self.entryname_=srec[0]
        self.refaa_=srec[2]
        self.altaa_=srec[3]
        #self.geneid_=srec[3]
        self.transcriptid_=srec[6]
        #self.genename_=srec[5]
        self.posaa_=int(srec[1])
        self.caddrawscore_=float(srec[4])
        self.caddphredscore_=float(srec[5])
        self.chrom_=srec[7]
        self.pos_=int(srec[8])
        self.ref_=srec[9]
        self.alt_=srec[10]

class CADDScoresDb(object):
    def __init__(self, refpathobj):
        self.caddtabixfile_=pysam.TabixFile(refpathobj.cadddb_)
    def get_cadd_score_records(self, q):
        tabrec=[] #string list
        try:
            tbi=self.caddtabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                ii=CADDScoreRecord(i)
                if ii.posaa_ == q.position_:
                    tabrec.append(ii)
        except ValueError:
            pass
        return tabrec
    def get_mutation_caddphredscore(self, q):
        crecs=self.get_cadd_score_records(q)
        csc=0
        for c in crecs:
            if c.altaa_ == q.altaa_:
                csc=c.caddphredscore_
        return csc
    def get_position_caddphredscore(self, q):
        crecs=self.get_cadd_score_records(q)
        a=AminoAcids()
        #extract all the mutations which are dissimilar aa types
        cscores=[]
        for c in crecs:
            if not a.is_similar( c.refaa_, c.altaa_ ):
                cscores.append(c.caddphredscore_)
        cscores=np.array(cscores)
        if len(cscores)==0:
            return 0
        else:
            return np.mean(cscores)
    def is_mutation_neutral(self, q):
        return self.get_mutation_caddphredscore(q) <= 10
    def is_position_neutral(self, q):
        return self.get_position_caddphredscore(q) <=10
    def get_reference_aminoacid(self, q):
        crecs=self.get_cadd_score_records(q)
        if len(crecs)>0:
            return crecs[0].refaa_
        else:
            return ""

class ClinvarSignificanceRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[0]
        self.posaa_=int(srec[1])
        self.refaa_=srec[2]
        self.altaa_=srec[3]
        self.referenceid_=srec[5]
        self.significance_=srec[6]

class ClinvarSignificaneDb(object):
    def __init__(self, refpathobj):
        self.clinvartabixfile_=pysam.TabixFile(refpathobj.clinvardb_)
    def get_clinvar_significance_records(self, q):
        srecords=[]
        try:
            tbi=self.clinvartabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                cs=ClinvarSignificanceRecord(i)
                if cs.posaa_ == q.position_:
                    if cs.altaa_ == q.altaa_:
                        srecords.append(cs)
        except ValueError:
            pass
        return srecords
    def has_clinvar_record(self, q):
        return len(self.get_clinvar_significance_records(q)) > 0
    def get_mutation_significance(self, q):
        if q.altaa_ == None:
            return "Unknown"
        srecords=self.get_clinvar_significance_records(q)
        if len(srecords)==0:
            return "Unknown"
        else:
            return srecords[0].significance_
    def get_reference_id(self, q):
        if q.altaa_ == None:
            return ""
        srecords=self.get_clinvar_significance_records(q)
        if len(srecords)==0:
            return ""
        else:
            return srecords[0].referenceid_

class DomainProteomeEntropyDb(object):
    def __init__(self, refpathobj):
        self.domainproteomentropytabixfile_= refpathobj.domainindexproteomeentropydbh_ #pysam.TabixFile(refpathobj.domainindexproteomeentropydb_)
        self.proteinposdomainindextabixfile_= refpathobj.proteinposdomainindexdbh_ #pysam.TabixFile(refpathobj.proteinposdomainindexdb_)
    def get_protein_position_entropy(self, entry, pos):
        entropy=0
        try:
            for p in self.proteinposdomainindextabixfile_.fetch(entry, pos-1, pos):
                sp=p.split("\t")
                if int(sp[1])==pos:
                    pfam=sp[2]
                    pind=int(sp[3])
                    for pe in self.domainproteomentropytabixfile_.fetch(pfam, pind-1, pind):
                        spe=pe.split("\t")
                        #print(spe)
                        if int(spe[1]) == pind:
                            entropy=float(spe[16])
                            break
                break
        except Exception as e:
            #print(e)
            pass
        return entropy

class COSMICCountRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[0]
        self.posaa_=int(srec[1])
        self.refaa_=srec[2]
        self.altaa_=srec[3]
        self.cosmicid_=srec[5]
        self.mcount_=int(srec[9])

class COSMICGenomeWideCountRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[1]
        self.posaa_=int(srec[7])
        self.refaa_=srec[6]
        self.altaa_=srec[8]
        self.cosmicid_=srec[4]
        self.mcount_=int(srec[0])

class COSMICGenomeWideDb(object):
    def __init__(self, refpathobj):
        self.cosmicgenometabixfile_=pysam.TabixFile(refpathobj.cosmicgenomewidedb_)
    def get_position_cosmic_count_records(self, q):
        ccr=[]
        try:
            tbi=self.cosmicgenometabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                cre=COSMICGenomeWideCountRecord(i)
                if cre.posaa_ == q.position_:
                    ccr.append( cre )
        except ValueError:
            pass
        return ccr
    def get_alternate_residue_frequency(self, q):
        p=self.get_position_cosmic_count_records(q)
        altaa_count={}
        tmcount=0
        if len(p)>0:
            for i in p:
                if i.posaa_ == q.position_:
                    altaa_count[i.altaa_]=i.mcount_
                    tmcount+=i.mcount_
        altaa_freq={}
        for a in altaa_count:
            altaa_freq[a]=altaa_count[a]/tmcount
        return altaa_freq
    def get_position_cosmic_count(self, q):
        p=self.get_position_cosmic_count_records(q)
        pc=0
        if len(p)>0:
            for i in p:
                if i.posaa_ == q.position_:
                    pc+=i.mcount_
        return pc
    def get_mutation_cosmic_count(self, q):
        p=self.get_position_cosmic_count_records(q)
        mc=0
        if len(p)>0:
            for i in p:
                if i.altaa_==q.altaa_:
                    mc=i.mcount_
        return mc
    def get_mutation_cosmicid(self, q):
        p=self.get_position_cosmic_count_records(q)
        cosmid=""
        if len(p)>0:
            for i in p:
                if i.altaa_==q.altaa_:
                    cosmid=i.cosmicid_
        return cosmid

class COSMICCountDb(object):
    def __init__(self, refpathobj):
        self.cosmictabixfile_=pysam.TabixFile(refpathobj.cosmicdb_)
    def get_position_cosmic_count_records(self, q):
        ccr=[]
        try:
            tbi=self.cosmictabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                cre=COSMICCountRecord(i)
                if cre.posaa_ == q.position_:
                    ccr.append( cre )
        except ValueError:
            pass
        return ccr
    def get_position_cosmic_count(self, q):
        p=self.get_position_cosmic_count_records(q)
        pc=0
        if len(p)>0:
            for i in p:
                if i.posaa_ == q.position_:
                    pc+=i.mcount_
        return pc
    def get_mutation_cosmic_count(self, q):
        p=self.get_position_cosmic_count_records(q)
        mc=0
        if len(p)>0:
            for i in p:
                if i.altaa_==q.altaa_:
                    mc=i.mcount_
        return mc
    def get_mutation_cosmicid(self, q):
        p=self.get_position_cosmic_count_records(q)
        cosmid=""
        if len(p)>0:
            for i in p:
                if i.altaa_==q.altaa_:
                    cosmid=i.cosmicid_
        return cosmid

class PhosphoSitePlusDbRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[0]
        self.pos_=int(srec[1])
        self.modtype_=srec[5]

class PhosphoSitePlusDb(object):
    def __init__(self, refpathobj):
        self.psptabixfile_=pysam.TabixFile(refpathobj.pspdb_)
    def get_psp_records_for_protein(self,q):
        precs=[]
        try:
            tbi=self.psptabixfile_.fetch(q.entryname_)
            for i in tbi:
                ps=PhosphoSitePlusDbRecord(i)
                precs.append(ps)
        except ValueError:
            pass
        return precs
    def protein_has_psp_records(self , q):
        return len(self.get_psp_records_for_protein(q)) > 0
    def get_psp_records(self, q):
        precs=[]
        try:
            tbi=self.psptabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                ps=PhosphoSitePlusDbRecord(i)
                if ps.pos_ == q.position_:
                    precs.append(ps)
        except ValueError:
            pass
        return precs
    def has_modification(self,q):
        return len(self.get_psp_records(q)) > 0

class UniprotMutagenesisDbRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[0]
        self.pos_=int(srec[1])
        self.altaa_=srec[3]
        self.desc_=srec[5]

class UniprotMutatgenesisDb(object):
    def __init__(self, refpathobj):
        self.uniprotmuttabixfile_=pysam.TabixFile(refpathobj.umutagendb_)
    def get_mutagenesis_records(self, q):
        precs=[]
        try:
            tbi=self.uniprotmuttabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                umd=UniprotMutagenesisDbRecord(i)
                if umd.pos_ == q.position_:
                    if umd.altaa_ == q.altaa_:
                        precs.append(umd.desc_)
        except ValueError:
            pass
        return precs
    def has_mutagenesis_record(self, q):
        return len(self.get_mutagenesis_records(q)) > 0

class UniprotSitesDbRecord(object):
    def __init__(self, rec):
        srec=rec.strip().split("\t")
        self.entryname_=srec[0]
        self.start_=int(srec[1])
        self.end_=int(srec[2])
        self.type_=srec[6]
        self.desc_=srec[8]

class UniprotSiteDb(object):
    def __init__(self, refpathobj):
        self.uniprotsitetabixfile_=pysam.TabixFile(refpathobj.usitesdb_)
    def get_site_records(self, q):
        precs=[]
        try:
            tbi=self.uniprotsitetabixfile_.fetch(q.entryname_, q.position_-1, q.position_)
            for i in tbi:
                us=UniprotSitesDbRecord(i)
                if us.type_ == "DISULFID": #exception to the DISULFID record as start and end are not truly what they are
                    if q.position_ == us.start_ or q.position_ == us.end_:
                        precs.append(us)
                else:
                    if us.start_ == us.end_:
                        if q.position_ == us.start_:
                            precs.append(us)
                    else:
                        if q.position_ in range(us.start_, us.end_):
                            precs.append(us)
        except ValueError:
            pass
        return precs
    def get_site_records_for_protein(self, q):
        precs=[]
        try:
            tbi=self.uniprotsitetabixfile_.fetch(q.entryname_)
            for i in tbi:
                us=UniprotSitesDbRecord(i)
                if us.type_ == "DISULFID": #exception to the DISULFID record as start and end are not truly what they are
                    if q.position_ == us.start_ or q.position_ == us.end_:
                        precs.append(us)
                else:
                    if us.start_ == us.end_:
                        if q.position_ == us.start_:
                            precs.append(us)
                    else:
                        if q.position_ in range(us.start_, us.end_):
                            precs.append(us)
                #if q.position_ in range(us.start_, us.end_):
                #    precs.append(us)
        except ValueError:
            pass
        return precs
    def protein_has_site_records(self, q):
        return len(self.get_site_records_for_protein(q)) > 0
    def has_site_records(self, q):
        return len(self.get_site_records(q)) > 0

#Write for 3d Structure parsing
'''
search if the query is near to the known analogous mutations, domains, significant mutations, uniprot sites, phosphositeplus sites
Two types of geometric proximity: short range, long range

Algorithm
- Load query 3d edge map
- Get the neighbors of the position in query
- search for the positions, neighbors in the statistically significant/ analogous mutation db, domains, PSP sites, uniprot sites
- define short and long proximal positions
'''

#It should be a series / list of querys against which the proximity search can be performed rather than individual queries
class ProteinNetworkException(Exception):
    def __init__(self, msg):
        super().__init__(msg)

class ProteinNetwork(object):
    def __init__(self, entry, refpathobj):
        self.proteingraph_=nx.Graph()
        self.entryname_=entry
        try:
            self.proteingraph_=nx.read_weighted_edgelist(refpathobj.edgedir_+self.entryname_+".edge")
        except Exception:
            raise ProteinNetworkException("%s protein structure could not be located!!" %(refpathobj.edgedir_+self.entryname_+".edge"))
    def get_neighbor_positions(self, q): #Given a query object return the neightbor positions as strings
        neighborpositions=[]
        if not q.entryname_ == self.entryname_:
            raise ProteinNetworkException("Different network loaded! Expected protein %s, loaded protein: %s" %(q.entryname_, self.entryname_))
        try:
            edges=self.proteingraph_.edges(str(q.position_))
            for e in edges:
                neighborpositions.append(int(e[1]))
        except Exception:
            raise ProteinNetworkException("Query position: %d , not found in graph!" %(q.position_))
        return neighborpositions
    def get_neighbor_positions_map(self, qlist):
        neighborpositions={} #position-> list of neighbors
        for q in qlist:
            if not q.entryname_ == self.entryname_:
                raise ProteinNetworkException("Different network loaded! Expected protein %s, loaded protein: %s" %(q.entryname_, self.entryname_))
            qneighbor=[]
            try:
                neighborpositions[str(q.position_)]=self.get_neighbor_positions(q)
                #edges=self.proteingraph_.edges(str(q.position_))
                #for e in edges:
                #    qneighbor.append(e[1])
                #neighborpositions[str(q.position_)]=qneighbor
            except Exception:
                raise ProteinNetworkException("Query position: %d , not found in graph!" %(q.position_))
        return neighborpositions #returns a map of neighbors for each query position provided
def initialize_reference_path_object(refpath):
    return RefPath(refpath)

#Complete proximity search operations on the structure
'''
Algorithm
-   Initialize reference object
-   load all the databases and store the objects
-   run through the queries and search neighbors in the protein graph
-   search the neighbors in the respective datasets for being any important residues
-   store and display the list as a table
'''
class ReferenceDatabases(object):
    def __init__(self, rpath):
        self.refobj=initialize_reference_path_object(rpath)
        self.uid_entryname_map_=UniprotIDEntryNameMap(self.refobj)
        self.proteininterfaces_=ProteinInterfaces(self.refobj)
        self.hotspot3ds_=Hotspot3dPositions(self.refobj)
        self.statsigs_=SignificantMutations(self.refobj)
        self.analogmuts_=AnalogousMutations(self.refobj)
        self.uniprotsitedbh_=UniprotSiteDb(self.refobj)
        self.uniprotpfamdbh_=UniprotPfamDb(self.refobj)
        self.uniprotmutagendbh_=UniprotMutatgenesisDb(self.refobj)
        self.clinvardbh_=ClinvarSignificaneDb(self.refobj)
        self.cosmiccountdbh_=COSMICCountDb(self.refobj)
        self.cosmicgenomewidedbh_=COSMICGenomeWideDb(self.refobj)
        self.phosphositeplusdbh_=PhosphoSitePlusDb(self.refobj)
        self.caddscoresdbh_=CADDScoresDb(self.refobj)
        self.edgedir_=self.refobj.edgedir_
        self.domaligndir_=self.refobj.domalignpath_
        self.domainindexproteomeentropydbh_=pysam.TabixFile(self.refobj.domainindexproteomeentropydb_)
        self.proteinposdomainindexdbh_=pysam.TabixFile(self.refobj.proteinposdomainindexdb_)

class DOMEOutputRecordParser(object):
    def __init__(self, line, sind=5, aind=6, cind=7, gsind=17, gaind=18 , siteind=[12,13,19,20,22] ):
        self.columns_=line.strip().split("\t")
        self.statind_=sind
        self.analogind_=aind
        self.caddind_=cind
        self.geostatind_=gsind
        self.geoanalogind_=gaind
        self.sitesind_=siteind
    def get_tier(self): #tier 1- statsig; tier 2: analogous, cadd > 20, sites; tier 3: analogous, cadd > 20,
        tier=0
        if self.columns_[self.statind_] != "-":
            tier=1
            #print(self.columns_[self.statind_])
        else:
            if self.columns_[self.analogind_] != "-" and float(self.columns_[self.caddind_]) > 20:
            #testing with CADD > 10 / non-neutral as potential candidate
            #if self.columns_[self.analogind_] != "-" and float(self.columns_[self.caddind_]) > 10:
                temp= np.array(self.columns_) #np.array(self.sitesind_)
                ###ERROR
                if len(list(np.where(temp[self.sitesind_] != "-")[0])) > 0:
                    tier=2
                else:
                    tier=3
            elif self.columns_[self.analogind_] == "-" and float(self.columns_[self.caddind_]) > 20:
            #elif self.columns_[self.analogind_] == "-" and float(self.columns_[self.caddind_]) > 10:
                if self.columns_[self.geostatind_] != "-" or self.columns_[self.geoanalogind_] != "-":
                    temp= np.array(self.columns_)#np.array(self.sitesind_)
                    if len(list(np.where(temp[self.sitesind_] != "-")[0])) > 0:
                        tier=4
        return tier

def main():
    #read and create a list of queries (query objects)
    #input csv file with entryname, pos, (altaa - if available)
    #'''
    infi=open(sys.argv[1])
    entry_queries={}
    for i in infi:
        if not i.startswith("#"):
            si=i.strip().split(",")
            qsi=Query(si[0], int(si[1]))
            if len(si)==3:
                qsi=Query(si[0], int(si[1]), si[2]) #If altered amino acid is provided then use it to create query
            if si[0] in entry_queries:
                temp=entry_queries[si[0]]
                temp.append(qsi)
                entry_queries[si[0]]=temp
            else:
                entry_queries[si[0]]=[ qsi ]
    refpath="../data/domerefdata/final/"
    oheader=["tier","protein","mutpos","refaa","altaa", "domain", "statscore","analogtomut", "caddscore", "cosmiccount", "colcoscount", "colnumofresidues" , "clinvarsig", "uniprotsite", "phosphosite", "uniprotmutagen" , "hotspot3dann", "3d_dom_close", "3d_statsig_close" , "3d_analog_close","3d_uniprotsite_close", "3d_pspsite_close", "3d_hotspot_close", "interfaceposition", "contextscore", "neutralityscore","entropyscore", "alterfrequencyscore","biochemcontext","meanscore"]
    refdbobj=ReferenceDatabases(refpath)
    domprotentropydb=DomainProteomeEntropyDb(refdbobj)
    cosmicgwtbx=refdbobj.cosmicgenomewidedbh_
    print("\t".join(oheader))
    #uncomment to print header
    AA_=AminoAcids()
    for ent in entry_queries:
        structnotfound=False
        try:
            entnet=ProteinNetwork(ent, refdbobj) #Load protein network
        except ProteinNetworkException:
            structnotfound=True
        #To be duplicated in proteomeaccessor
        upfamrec=refdbobj.uniprotpfamdbh_.get_pfam_domain_records_for_protein_entry(ent)
        domname_alignments={} #stores domain alignments for a particular domain name (in this case for a protein)
        for upfa in upfamrec:
            try:
                domname_alignments[upfa.domainname_]= DomainMultiAlignedSequence(upfa.domainid_, refdbobj.domaligndir_+upfa.domainid_+".aln.fa")
            except FileNotFoundError:
                pass
        ######
        if not structnotfound:
            for q in entry_queries[ent]:
                '''
                #semicolon sep - SS
                outvec is output vector which contains following elements
                entryname, position, refaa, altaa, pfname
                (statscores(SS), stattumours (SS)) | (statsig to analogous, statsig to analog position, statsig to analogpos residue),
                caddscore (neutrality),
                cosmiccount, clinvarsig,
                uniprotsite, pspsite, uniprotmutagenesis
                '''
                outvec=[]
                out3dvec=[]
                q1chemcontext=[]
                q2refchemcontext=[]
                qneighbors=entnet.get_neighbor_positions(q) #return int positions
                outvec.append(q.entryname_)
                outvec.append(str(q.position_))
                outvec.append(refdbobj.caddscoresdbh_.get_reference_aminoacid(q)) #loding sequeces might help
                outvec.append(q.altaa_)
                for qn1 in qneighbors:
                    raa=refdbobj.caddscoresdbh_.get_reference_aminoacid( Query(q.entryname_, qn1) )
                    if raa != '':
                        q1chemcontext.append( raa )
                totcoscount=0
                colresidues=0
                if refdbobj.uniprotpfamdbh_.is_within_domain(q):
                    domname=refdbobj.uniprotpfamdbh_.get_pfam_domain_name(q)
                    domaln=[]
                    try:
                        domaln=domname_alignments[domname]
                    except KeyError:
                        pass
                    daccs=refdbobj.uniprotpfamdbh_.get_pfam_domain_accessions(q)
                    if len(daccs) > 1:
                        #print(daccs)
                        daccs=list(np.unique(np.array(daccs)))
                        #cant handle this currently
                        #sys.exit(0)
                    if len(daccs)==1:
                        dacc=daccs[0]
                        try:
                            anposs=domaln.get_analogous_positions(dacc, q.position_)
                            for aaa in anposs:
                                aaentryname = refdbobj.uid_entryname_map_.get_entryname_by_uniprotid(aaa)
                                aaentrypos = anposs[aaa]
                                qq=Query(aaentryname, aaentrypos)
                                ccount=refdbobj.cosmiccountdbh_.get_position_cosmic_count(qq)
                                totcoscount+=ccount
                                colresidues+=1
                        except IndexError:
                            pass
                        except KeyError:
                            pass
                        except Exception:
                            pass
                    outvec.append(domname)
                    #upfamrec=refdbobj.uniprotpfamdbh_.get_pfam_domain_records_for_protein(q)
                    neighbordomains=[]
                    #in later versions you may want to remove the same domain name in 3d output vector, as it is redundant
                    for qn in qneighbors:
                        for upf in upfamrec:
                            if qn in range(upf.domstart_, upf.domend_+1):
                                if not upf.domainname_ in neighbordomains:
                                    neighbordomains.append(upf.domainname_)
                    if len(neighbordomains) > 0:
                        out3dvec.append("%s" %(";".join(neighbordomains)))
                    else:
                        out3dvec.append('-')
                else: #check if geometrically connected to domain
                    outvec.append('-') #since not within domain print -
                    upfamrec=refdbobj.uniprotpfamdbh_.get_pfam_domain_records_for_protein(q)
                    neighbordomains=[]
                    for qn in qneighbors:
                        for upf in upfamrec:
                            if qn in range(upf.domstart_, upf.domend_+1):
                                if not upf.domainname_ in neighbordomains:
                                    neighbordomains.append(upf.domainname_)
                    if len(neighbordomains) > 0:
                        out3dvec.append("%s" %(";".join(neighbordomains)))
                    else:
                        out3dvec.append('-')
                sigstatus=False
                if refdbobj.statsigs_.is_position_significant(q): #if stat significant
                    sigstatus=True
                    #sscores_=refdbobj.statsigs_.get_significance_scores(q)
                    sscoresstr=refdbobj.statsigs_.get_significance_scores(q)
                    sstumours=refdbobj.statsigs_.get_tumour_types(q)
                    #Since tumors are not stored in statsig file, it is commented in the current version
                    outvec.append(";".join(sscoresstr))
                    #closer to any other stat sig in 3d
                    pstatsigpos=refdbobj.statsigs_.get_significant_positions_for_protein(q)
                    pstattumtypes=refdbobj.statsigs_.get_tumour_types_for_protein(q)
                    siggeoclose=[]
                    sigtumorgeoclose=[]
                    for pstatind in range(len(pstatsigpos)):
                        if (pstatsigpos[pstatind] in qneighbors):
                            siggeoclose.append(str(pstatsigpos[pstatind]))
                            sigtumorgeoclose.append(pstattumtypes[pstatind])
                    if len(siggeoclose)>0:
                        #out3dvec.append("%s|%s" %(";".join(siggeoclose), ";".join(sigtumorgeoclose)) )
                        out3dvec.append("%s" %(";".join(siggeoclose)) )
                    else:
                        out3dvec.append('-')
                else: #if not statistically significant
                    outvec.append("-")
                #else: #check if analogous and #check if geometrically closer to stat sig position
                if not sigstatus:
                    pstatsigpos=refdbobj.statsigs_.get_significant_positions_for_protein(q)
                    pstattumtypes=refdbobj.statsigs_.get_tumour_types_for_protein(q)
                    siggeoclose=[]
                    sigtumorgeoclose=[]
                    for pstatind in range(len(pstatsigpos)):
                        if (pstatsigpos[pstatind] in qneighbors):
                            siggeoclose.append(str(pstatsigpos[pstatind]))
                            sigtumorgeoclose.append(pstattumtypes[pstatind])
                    if len(siggeoclose)>0:
                        #out3dvec.append("%s|%s" %(";".join(siggeoclose), ";".join(sigtumorgeoclose)) )
                        out3dvec.append("%s" %(";".join(siggeoclose)) )
                    else:
                        out3dvec.append('-')
                #check if geometrically closer to analogous positions
                aposs=refdbobj.analogmuts_.get_analogous_positions_for_protein(q.entryname_)
                anpos3d=[]
                #print(aposs)
                #print(qneighbors)
                for ap in aposs:
                    if ap in qneighbors:
                        akey=q.entryname_+";"+str(ap)
                        anpos3d.append(akey)
                if len(anpos3d) > 0:
                    out3dvec.append( "|".join(anpos3d) )
                else:
                    out3dvec.append("-")
                if refdbobj.analogmuts_.is_position_analogous(q):
                    #xxx
                    sigp=refdbobj.analogmuts_.get_significant_protein_to_analogous_position(q)
                    sigps=refdbobj.analogmuts_.get_significant_protein_position_to_analogous_position(q)
                    sigr=refdbobj.analogmuts_.get_significant_protein_residue_to_analogous_position(q)
                    stemp=[]
                    for ie in range(0,len(sigp)):
                        stemp.append(str(sigp[ie])+";"+str(sigps[ie])+";"+sigr[ie])
                    #outvec.append(str(sigp)+";"+str(sigps)+";"+sigr)
                    outvec.append("|".join(stemp))
                else: #if nothing put a blank
                    outvec.append("-")
                cdscore=0
                if q.altaa_ =='':
                    cdscore=refdbobj.caddscoresdbh_.get_position_caddphredscore(q)
                else:
                    cdscore=refdbobj.caddscoresdbh_.get_mutation_caddphredscore(q)
                outvec.append(str(format(cdscore, '.2f')))
                if q.altaa_ == '':
                    outvec.append( str(refdbobj.cosmiccountdbh_.get_position_cosmic_count(q)) )
                else:
                    outvec.append( str(refdbobj.cosmiccountdbh_.get_mutation_cosmic_count(q)) )
                #column cosmic count

                if totcoscount ==0:
                    outvec.append("-")
                    outvec.append("-")
                else:
                    outvec.append(str( totcoscount ))
                    outvec.append(str(colresidues))
                #Check with clinvar
                if refdbobj.clinvardbh_.has_clinvar_record(q):
                    outvec.append(refdbobj.clinvardbh_.get_mutation_significance(q))
                else:
                    outvec.append('-')
                #uniprot sites
                usites=[]
                usitegeoclose=[]
                usitedescgeoclose=[]
                if refdbobj.uniprotsitedbh_.protein_has_site_records(q):
                    for qn in qneighbors:
                        for s in refdbobj.uniprotsitedbh_.get_site_records_for_protein(q):
                            if s.type_ == "DISULFID":
                                if qn == s.start_ or qn == s.end_:
                                    usitegeoclose.append(str(qn))
                                    #usitedescgeoclose.append(s.desc_)
                                    usitedescgeoclose.append(s.type_)
                            else:
                                if qn in range(s.start_, s.end_+1) :
                                    usitegeoclose.append(str(qn))
                                    #usitedescgeoclose.append(s.desc_)
                                    usitedescgeoclose.append(s.type_)
                                    break
                        break
                    out3dvec.append("%s|%s" %(";".join(usitegeoclose),";".join(usitedescgeoclose)))
                else:
                    out3dvec.append('-')
                if refdbobj.uniprotsitedbh_.has_site_records(q):
                    for s in refdbobj.uniprotsitedbh_.get_site_records(q):
                        #usites(s.desc_)
                        usites.append(s.type_)
                if len(usites)>0:
                    outvec.append(";".join(usites))
                else:
                    outvec.append("-")
                #phosphositeplus
                psites=[]
                psitegeoclose=[]
                psitetypegeoclose=[]
                if refdbobj.phosphositeplusdbh_.protein_has_psp_records(q):
                    precs=refdbobj.phosphositeplusdbh_.get_psp_records_for_protein(q)
                    if len(precs) > 0:
                        for qn in qneighbors:
                            for pr in precs:
                                if qn == pr.pos_:
                                    psitegeoclose.append(str(qn))
                                    psitetypegeoclose.append(pr.modtype_)
                        if len(psitegeoclose) > 0:
                            out3dvec.append("%s|%s" %(";".join(psitegeoclose), ";".join(psitetypegeoclose)))
                        else:
                            out3dvec.append('-')
                    else:
                        out3dvec.append('-')
                else:
                    out3dvec.append('-')
                #Check for actual match
                if refdbobj.phosphositeplusdbh_.has_modification(q):
                    modtypes=[]
                    precs=refdbobj.phosphositeplusdbh_.get_psp_records(q)
                    for pr in precs:
                        modtypes.append(pr.modtype_)
                    outvec.append(";".join(modtypes))
                else:
                    outvec.append("-")
                #try: #for value error thrown
                if refdbobj.uniprotmutagendbh_.has_mutagenesis_record(q):
                    mrecs=[]
                    for mr in refdbobj.uniprotmutagendbh_.get_mutagenesis_records(q):
                        #mrecs.append(mr.desc_)
                        mrecs.append(mr)
                    outvec.append(";".join(mrecs))
                else:
                    outvec.append("-")
                #Check for actual match in hotspot 3d
                if refdbobj.hotspot3ds_.is_hotspot3d_position(q):
                    rh=refdbobj.hotspot3ds_.get_hotspot3d_position_records(q)
                    try:
                        hclist0=list(rh.hclass)[0]
                        outvec.append(list(rh.hclass)[0] )
                    except Exception as err:
                        outvec.append("-")
                else:
                    outvec.append("-") #not 3d hotspot linked
                #Check for geometrical closeness to hotspot 3d
                hotspotmutspos=[]
                if refdbobj.hotspot3ds_.protein_has_hotspot3d_positions(q):
                    phps=refdbobj.hotspot3ds_.get_hotspot3d_positions_for_protein(q)
                    for qn in qneighbors:
                        if qn in phps and not (qn in hotspotmutspos):
                            hotspotmutspos.append(str(qn))
                if len(hotspotmutspos) > 0:
                    out3dvec.append("%s" %(";".join(hotspotmutspos)))
                else:
                    out3dvec.append("-")
                if q.position_ in refdbobj.proteininterfaces_.get_interface_positions_by_entryname(q.entryname_):
                    out3dvec.append("I") #Is an interface
                else:
                    out3dvec.append('-')
                dorec="\t".join(outvec) + "\t" + "\t".join(out3dvec)
                dr=DOMEOutputRecordParser(dorec)
                #print("%d\t%s" %(dr.get_tier(), dorec))
                #Changed code below this-----
                #if analogous;
                    #check for similar context
                dorec=dorec.split("\t")
                scores=[0.0,0.0,0.0,0.0,0.0]
                if dorec[6] != "-":
                    #scores=[0,0,0,0] #context, neutrality, entropy, cosmic context freq, context
                    ssigprofilematch={} #every profile contains; 1. PTM, 2. geometrically closer to PTM, AA, interface / sum is stored (4 is max)
                    statsigqueries=[] #list of query objects, representing statsig
                    #print(dorec)
                    #print("$$")
                    ssigs=dorec[6].split('|')
                    #print(ssigs)
                    sprofiles=[]
                    sdists=[]
                    for s in ssigs:
                        ss=s.split(";")
                        sprofile=[0,0,0,0]
                        #sprofile=[0,0,0,0,0]
                        #print(ss)
                        qs=Query(ss[0], int(ss[1]))
                        qstructnotfound=False
                        try:
                            qsentnet=ProteinNetwork(qs.entryname_, refdbobj) #Load protein network
                        except ProteinNetworkException:
                            qstructnotfound=True
                        if not qstructnotfound:
                            if ( (dorec[12]!="-") == (refdbobj.uniprotsitedbh_.has_site_records(qs)) ) or ( (dorec[13]!="-") == (refdbobj.phosphositeplusdbh_.has_modification(qs)) ):
                                sprofile[0]=1
                            qsunigeo=False
                            qspspgeo=False
                            qsusitegeoclose=[]
                            qsneighbors=qsentnet.get_neighbor_positions(qs)
                            #print(qsneighbors)
                            for qn2 in qsneighbors:
                                raa2=refdbobj.caddscoresdbh_.get_reference_aminoacid( Query( qs.entryname_ , qn2) )
                                if raa2 != '':
                                    q2refchemcontext.append( raa2 )
                            if refdbobj.uniprotsitedbh_.protein_has_site_records(qs):
                                for qsn in qsneighbors:
                                    for s in refdbobj.uniprotsitedbh_.get_site_records_for_protein(qs):
                                        if s.type_ == "DISULFID":
                                            if qsn == s.start_ or qsn == s.end_:
                                                qsusitegeoclose.append(str(qsn))
                                        else:
                                            if qsn in range(s.start_, s.end_+1) :
                                                qsusitegeoclose.append(str(qsn))
                                                break
                                    break
                                if len(qsusitegeoclose)>0:
                                    qsunigeo=True
                            qspsitegeoclose=[]
                            if refdbobj.phosphositeplusdbh_.protein_has_psp_records(qs):
                                precs=refdbobj.phosphositeplusdbh_.get_psp_records_for_protein(qs)
                                if len(precs) > 0:
                                    for qsn in qsneighbors:
                                        for pr in precs:
                                            if qsn == pr.pos_:
                                                qspsitegeoclose.append(str(qsn))
                                    if len(qspsitegeoclose) > 0:
                                        qspspgeo=True
                            if (qsunigeo or qspspgeo) == (dorec[19] != "-" or dorec[20] != "-"):
                                sprofile[1]=1
                            if AA_.is_similar(ss[2], dorec[2]):
                                sprofile[2]=1
                            if (qs.position_ in refdbobj.proteininterfaces_.get_interface_positions_by_entryname(qs.entryname_) ) == ( int(dorec[1]) in refdbobj.proteininterfaces_.get_interface_positions_by_entryname(dorec[0]) ) :
                                sprofile[3]=1
                            statsigqueries.append(sum(sprofile))
                        if len(q1chemcontext) != 0 and len(q2refchemcontext) != 0:
                            q1fr=AA_.get_amino_frequency(q1chemcontext)
                            q2fr=AA_.get_amino_frequency(q2refchemcontext)
                            q1aafreq=[]
                            q2aafreq=[]
                            for q1a in q1fr:
                                q1aafreq.append(q1fr[q1a])
                                q2aafreq.append(q2fr[q1a])
                            sdists.append(np.linalg.norm( np.array(q1aafreq) - np.array(q2aafreq)))

                        sprofiles.append(sum(sprofile)/len(sprofile))

                    #print("#sprofile")
                    #print(sprofile)
                    scores[0]= max(sprofiles)
                    #if sum(sprofile) == 4: #profile match
                    #    scores[0]=1
                    if float(dorec[7]) > 20:
                        scores[1]=1
                    entropy=domprotentropydb.get_protein_position_entropy(dorec[0], int(dorec[1]))
                    scores[2]=entropy
                    aachangefreq=0
                    try:
                        ppfamid=""
                        ppindex=0
                        for ppdi in domprotentropydb.proteinposdomainindextabixfile_.fetch(dorec[0], int(dorec[1])-1, int(dorec[1])):
                            sppdi=ppdi.split("\t")
                            if sppdi[1] == dorec[1]:
                                ppfamid=sppdi[2]
                                ppindex=int(sppdi[3])
                                break
                        if ppfamid != "":
                            if ppindex > 0:
                                for enrec in domprotentropydb.domainproteomentropytabixfile_.fetch(ppfamid, ppindex-1, ppindex):
                                    senrec=enrec.split("\t")
                                    #print(senrec)
                                    if senrec[0]==ppfamid+";"+str(ppindex):
                                        x=senrec[6]+"|"+senrec[10]
                                        x=x.split("|")
                                        if len(x) != 0:
                                            #print(x)
                                            modx=[]
                                            for ix in x:
                                                if ix != "":
                                                    modx.append(ix)
                                            x=modx
                                            aa=AminoAcidCosmicChangeFrequency(x, cosmicgwtbx.cosmicgenometabixfile_ , AA_)
                                            #print(aa.get_typechange_frequency(dorec[2],dorec[3], AA_))
                                            #print("----")
                                            aachangefreq = float("%.4f" %(aa.get_typechange_frequency(dorec[2],dorec[3], AA_)))
                                            #print(aachangefreq)
                                        else:
                                            aachangefreq=0
                            else:
                                for enrec in domprotentropydb.domainproteomentropytabixfile_.fetch(ppfamid, ppindex, ppindex):
                                    senrec=enrec.split("\t")
                                    #print(senrec)
                                    if senrec[0]==ppfamid+";"+str(ppindex):
                                        x=senrec[6]+"|"+senrec[10]
                                        x=x.split("|")
                                        if len(x) != 0:
                                            #print(x)
                                            modx=[]
                                            for ix in x:
                                                if ix != "":
                                                    modx.append(ix)
                                            x=modx
                                            aa=AminoAcidCosmicChangeFrequency(x, cosmicgwtbx.cosmicgenometabixfile_, AA_)
                                            aachangefreq = float("%.4f" %(aa.get_typechange_frequency(dorec[2],dorec[3], AA_)) )
                                        else:
                                            aachangefreq=0
                        '''
                        else:
                            print(ppfamid)
                            print(ppindex)
                            print(dorec)
                        '''
                    except Exception as e:
                        pass
                        #print(e)
                        #print(dorec)
                        #print("nothing of aa change frequency worked!")
                        #sys.exit(0)
                    scores[3]=aachangefreq
                #print(q1chemcontext)
                #print(q2refchemcontext)
                    if len(sdists) != 0:
                        scores[4]=1-min(sdists)
                #print(scores)
                meanscore=(scores[0]+scores[2]+scores[3]+scores[4])/4
                print("%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" %(dr.get_tier(), "\t".join(dorec), scores[0],scores[1],scores[2],scores[3], scores[4], meanscore))
                #print("%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" %(dr.get_tier(), "\t".join(dorec), scores[0],scores[1],scores[2],scores[3], scores[4], sum(scores)/len(scores)))
                #uncomment
                #dr=DOMEOutputRecordParser(dorec, statind, analogind, caddind, geostatind, geoanalogind, sitesind)
                #print("%s\t%s" %("\t".join(outvec), "\t".join(out3dvec) ) )
                #Needs to be uncommented
                #print( outvec )
                #print( out3dvec )


if __name__=="__main__":
    main()
