from domainconservationscore import *
from statisticallysignificantmutations import *
from uniprot import *
from aacids import *
import sys
import os

class AnalogousDomainPosition(object):
    def __init__(self, ape, apos, ares, asc, di, dn, ent, epo, ere ):
        self.analogprotein_entryname_=ape
        self.analogprotein_pos_=apos
        self.analogprotein_res_=ares
        self.analogstatscore_=asc
        self.domainid_=di
        self.domainname_=dn
        self.statsig_entryname_=ent
        self.statsig_pos_=epo
        self.statsig_res_=ere
    def __str__(self):
        return self.analogprotein_entryname_+"\t"+str(self.analogprotein_pos_)+"\t"+self.analogprotein_res_+"\t"+str(self.analogstatscore_)+"\t"+self.domainid_+"\t"+self.domainname_+"\t"+self.statsig_entryname_+"\t"+str(self.statsig_pos_)+"\t"+self.statsig_res_
def nonconserved_analogouspositions_from_hotspots(tcgastatsfile, genedomainmap, uniprotobj, domalndir):
    anadomposobjs=[]
    phm=ProteinHotspotMutations(tcgastatsfile)
    pdomrec=DomainRecords(genedomainmap)
    counter=0
    entrykey_hotspots={}
    aa=AminoAcids()
    for ename in uniprotobj.get_protein_entrynames():
        eprot=uniprotobj.get_protein_object_by_entryname(ename)
        plen=len(eprot)
        #print("processing %d %s " %( counter, ename))
        counter+=1
        #scoretrack=[0]*plen
        accs=uniprotobj.get_accessions_of_protein(ename)
        #analogmap={}
        #sigmut_analogmap={} #statistially significant mutation key -> anlog map
        tcgaacc=''
        hmuts=[]
        for a in accs:
            if phm.has_hotspot_mutations_by_accession(a):
                tcgaacc=a
                hmuts=phm.get_hotspot_mutations_by_accession(tcgaacc)
                break
        for h in hmuts:
            kk=ename+"_"+str(h.position_) #entryname_position key
            if not kk in entrykey_hotspots:
                entrykey_hotspots[kk]=h
        '''
            try:
                scoretrack[h.position_-1]=1
            except IndexError as ie:
                pass
        protein_mutsig_score_track[ename]=scoretrack
        protein_analogsig_score_track[ename]=[0]*plen
        '''
        #significant positions added
    #Going for analogous positions
    #for each of the position if inside domain check for analogous positions
    for ek in entrykey_hotspots:
        ekhot=entrykey_hotspots[ek]
        sek=ek.split("_")
        enm=sek[0]+"_"+sek[1] #stat sig entryname
        eaccs=uniprotobj.get_accessions_of_protein(enm)
        ppos=int(sek[2]) #stat sig pos
        try:
            eres=uniprotobj.get_residue_by_entry(enm, ppos)
            usefulacc='' #accession which is part of the domain file should be stored here
            analogmap={}
            if ekhot.is_within_domain():
                domainid=ekhot.domain_id_
                domalignment=domalndir+domainid+".aln.fa"
                try:
                    doma=DomainMultiAlignedSequence(domainid, domalignment)
                    for ea in eaccs:
                        if doma.has_sequence(ea):
                            analogmap=doma.get_analogous_positions(ea, ppos)
                except FileNotFoundError:
                    pass
                except IndexError: #when position is not found within sequence
                    pass
                #print(analogmap)
                if len(analogmap) > 0:
                    for an in analogmap:
                        if analogmap[an]>=0:
                            #print(analogmap[an])
                            try: #get protein entry name using accession
                                enn=uniprotobj.get_protein_object_by_accession(an)
                                #get residue of analogous position and position here and create an object / AnalogousDomainPosition
                                anares=enn.sequence_[ analogmap[an] -1]
                                anapos= analogmap[an]
                                #domname=pdomrec.get_domain_name_from_id(domainid)
                                #print(domname)
                                domname="-"
                                #statsig and analog should not be same
                                if not enn.entry_name_ == enm:
                                    #and the residue is similar
                                    if not aa.is_similar(anares, eres):
                                        #print( "%s %d %s %s %s %s %s %d %s" %(enn.entry_name_, anapos, anares, "0.6", domainid, domname, enm, ppos, eres))
                                        dpobj=AnalogousDomainPosition(enn.entry_name_, anapos, anares, 0.3, domainid, domname, enm, ppos, eres)
                                        anadomposobjs.append(dpobj)
                                #ennstrack=protein_analogsig_score_track[enn.entry_name_]
                                #ennstrack[analogmap[an] - 1]=ennstrack[ analogmap[an] -1]+0.6
                                #statanmapfo.write("%s\t%s\n" %(ek, enn.entry_name_+"_"+str(analogmap[an])))
                                #protein_analogsig_score_track[enn.entry_name_]=ennstrack
                            except KeyError: #accession not found - ignore
                                pass
                            except IndexError: #domain analogous position goes beyond protein
                                pass
        except IndexError:
            pass
    #statanmapfo.close()
    return anadomposobjs

def analogouspositions_from_hotspots(tcgastatsfile, genedomainmap, uniprotobj, domalndir):#returns an array of output
    anadomposobjs=[]
    phm=ProteinHotspotMutations(tcgastatsfile)
    pdomrec=DomainRecords(genedomainmap)
    counter=0
    entrykey_hotspots={}
    aa=AminoAcids()
    for ename in uniprotobj.get_protein_entrynames():
        eprot=uniprotobj.get_protein_object_by_entryname(ename)
        plen=len(eprot)
        #print("processing %d %s " %( counter, ename))
        counter+=1
        #scoretrack=[0]*plen
        accs=uniprotobj.get_accessions_of_protein(ename)
        #analogmap={}
        #sigmut_analogmap={} #statistially significant mutation key -> anlog map
        tcgaacc=''
        hmuts=[]
        for a in accs:
            if phm.has_hotspot_mutations_by_accession(a):
                tcgaacc=a
                hmuts=phm.get_hotspot_mutations_by_accession(tcgaacc)
                break
        for h in hmuts:
            kk=ename+"_"+str(h.position_) #entryname_position key
            if not kk in entrykey_hotspots:
                entrykey_hotspots[kk]=h
        '''
            try:
                scoretrack[h.position_-1]=1
            except IndexError as ie:
                pass
        protein_mutsig_score_track[ename]=scoretrack
        protein_analogsig_score_track[ename]=[0]*plen
        '''
        #significant positions added
    #Going for analogous positions
    #for each of the position if inside domain check for analogous positions
    for ek in entrykey_hotspots:
        ekhot=entrykey_hotspots[ek]
        sek=ek.split("_")
        enm=sek[0]+"_"+sek[1] #stat sig entryname
        eaccs=uniprotobj.get_accessions_of_protein(enm)
        ppos=int(sek[2]) #stat sig pos
        try:
            eres=uniprotobj.get_residue_by_entry(enm, ppos)
            usefulacc='' #accession which is part of the domain file should be stored here
            analogmap={}
            if ekhot.is_within_domain():
                domainid=ekhot.domain_id_
                domalignment=domalndir+domainid+".aln.fa"
                try:
                    doma=DomainMultiAlignedSequence(domainid, domalignment)
                    for ea in eaccs:
                        if doma.has_sequence(ea):
                            analogmap=doma.get_analogous_positions(ea, ppos)
                except FileNotFoundError:
                    pass
                except IndexError: #when position is not found within sequence
                    pass
                #print(analogmap)
                if len(analogmap) > 0:
                    for an in analogmap:
                        if analogmap[an]>=0:
                            #print(analogmap[an])
                            try: #get protein entry name using accession
                                enn=uniprotobj.get_protein_object_by_accession(an)
                                #get residue of analogous position and position here and create an object / AnalogousDomainPosition
                                anares=enn.sequence_[ analogmap[an] -1]
                                anapos= analogmap[an]
                                #domname=pdomrec.get_domain_name_from_id(domainid)
                                #print(domname)
                                domname="-"
                                #statsig and analog should not be same
                                if not enn.entry_name_ == enm:
                                    #and the residue is similar
                                    if aa.is_similar(anares, eres):
                                        #print( "%s %d %s %s %s %s %s %d %s" %(enn.entry_name_, anapos, anares, "0.6", domainid, domname, enm, ppos, eres))
                                        dpobj=AnalogousDomainPosition(enn.entry_name_, anapos, anares, 0.6, domainid, domname, enm, ppos, eres)
                                        anadomposobjs.append(dpobj)
                                #ennstrack=protein_analogsig_score_track[enn.entry_name_]
                                #ennstrack[analogmap[an] - 1]=ennstrack[ analogmap[an] -1]+0.6
                                #statanmapfo.write("%s\t%s\n" %(ek, enn.entry_name_+"_"+str(analogmap[an])))
                                #protein_analogsig_score_track[enn.entry_name_]=ennstrack
                            except KeyError: #accession not found - ignore
                                pass
                            except IndexError: #domain analogous position goes beyond protein
                                pass
        except IndexError:
            pass
    #statanmapfo.close()
    return anadomposobjs

def main():
    proteomeswissprotfile="../data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
    dalndir="../data/domerefdata/final/mafft_outaln/"
    #genedomainmap="../data/mappings/gene_domain_map.txt"
    genedomainmap="../data/uniprot/uniacc_pfam.tsv"
    uprot=UniprotProteomeParser(proteomeswissprotfile)
    #apfh=analogouspositions_from_hotspots(sys.argv[1], genedomainmap, uprot, dalndir)
    napfh=nonconserved_analogouspositions_from_hotspots(sys.argv[1], genedomainmap, uprot, dalndir)
    h="entry\tpos\tresidue\tstatscore\tdomainid\tdomainname\tstatsigentry\tstatpos\tstatsigres"
    print(h)
    #for a in apfh:
    #    print(str(a))
    for na in napfh:
        print(str(na))
if __name__=="__main__":
    main()
