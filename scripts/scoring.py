#!/usr/bin/python3
#Script to create track scores for indivual mutation residue scoring
#1 domain conservation score for individual position
#2 cadd score for individual mutation
#3 mutated or non-mutated statistically significant (analogous to one of such mutations)
#4 Functionally important positions to be scored as 1 or 0 based on presence or abscence

import sys
import os
from functionalsitescore import *
from domainconservationscore import *
from statisticallysignificantmutations import *
from uniprot import *
from caddrecordparser import *
import pysam
from phosphositeplusmodresscoring import *

proteomeswissprotfile="../data/uniprot/uniprot_hs_proteome_12_4_2021.txt"
#genedomainmap="../data/mappings/gene_domain_map.txt"
genedomainmap="../data/uniprot/uniacc_pfam.tsv"
domalnscoredir="../data/domerefdata/final/mafft_cons_score/"
domalndir="../data/domerefdata/final/mafft_outaln/"
functionalsitefile="../data/uniprot/uniprot_hs_proteome_12_4_2021_positions_impsites_withensg_uniq.txt"
hgncuniprotmapfile='../data/mappings/hgnc_uniprot_mapping2.tsv'
tcgastatsfile="../data/tcga_stats/Stats_Sig_Protein_db_all_cancer.tsv"
caddmissense="../data/cadd_download/cadd_missense_sorted.tsv"
phosphositeplusfile="../data/phosphositeplus/phosphositeplus_all_modres_human.tsv"
# all the scoring tracks would be of following file format:
# proteinid scores,comma,separated  totalscore length

def print_scoring_track(trackmap, outputfilehandle):
    for t in trackmap:
        totalscore=0
        tlen=len(trackmap[t])
        scorearr=[]
        for score in trackmap[t]:
            totalscore+=score
            scorearr.append(str(score))
        outputfilehandle.write("%s\t%s\t%d\t%d\n" %(t, ','.join(scorearr), totalscore, tlen ) )

def phoshositeplus_functional_scoring(uniprotobj):
    entryname_scoretrack={}
    phosphositeobj=PhosphositeplusParser(phosphositeplusfile)
    for ename in uniprotobj.get_protein_entrynames():
        #print(ename)
        plen=uniprotobj.get_protein_length(ename)
        scoretrack=[0]*plen
        entryname_scoretrack[ename]=scoretrack
    for uacc in phosphositeobj.phosphositeplus_functional_site_map_:
        try: #check if you can find entry name for accession
            pr=uniprotobj.get_protein_object_by_phosphositeplusid(uacc)
            scoretrack=entryname_scoretrack[pr.entry_name_]
            psplussites=phosphositeobj.phosphositeplus_functional_site_map_[uacc]
            for p in psplussites:
                try:
                    scoretrack[p.get_position()-1]=1
                except IndexError:
                    pass
            entryname_scoretrack[pr.entry_name_]=scoretrack
            #print("Processed %s" %(pr.entry_name_))
        except Exception as e:
            #print("Accession %s not found in the human proteome" %(uacc))
            for i in phosphositeobj.phosphositeplus_functional_site_map_[uacc]:
                print(str(i))
            #sys.exit(0)
            pass
    return entryname_scoretrack

def cadd_scoring(uniprotobj):
    #caddmissensefi=open(caddmissense)
    ensgsnotfound=[]
    entryname_scoretrack={} #first store addition of score then mean score for a position
    entryname_scorecounts={} #how many mutations at a protein position
    ensgstartendfi=open("../data/cadd_download/cadd_gene_start_end.tsv")
    ensg_tabquery_map={} #all string based
    for ese in ensgstartendfi:
        ese=ese.strip().split('\t')
        if not ese[0] in ensg_tabquery_map:
            ensg_tabquery_map[ese[0]]= (ese[1],int(ese[2]),int(ese[3]))
    #write a tabix based record extractor
    tb=pysam.TabixFile("../data/cadd_download/cadd_missense_sorted.tsv.gz")
    for en in ensg_tabquery_map:
        try:
            x=ensg_tabquery_map[en]
            recs=tb.fetch(x[0],x[1],x[2])
            upro=uniprotobj.get_protein_object_by_ensgid(en)
            scoretrack=[0]*len(upro)
            counttrack=[0]*len(upro)
            entryname_scoretrack[upro.entry_name_]=scoretrack
            entryname_scorecounts[upro.entry_name_]=counttrack
            for r in recs:
                c=CADDRecord(r)
                try:
                    scoretrack[int(c.protpos)-1]=scoretrack[int(c.protpos)-1]+float(c.rawscore)
                    counttrack[int(c.protpos)-1]=counttrack[int(c.protpos)-1]+1
                except IndexError: #position out of bounds for score track, ignore it
                    pass
            entryname_scoretrack[upro.entry_name_]=scoretrack
            entryname_scorecounts[upro.entry_name_]=counttrack
            print("Processed gene: %s, uniprot entry: %s" %(en, upro.entry_name_))
        except KeyError:
            print("%s not found in the uniprot proteome!" %(en))
    '''
    for crec in caddmissensefi:
        c=CADDRecord(crec)
        try:
            entryname=uniprotobj.ensgid_entryname_[c.geneid]
            if not entryname in entryname_scoretrack:
                print("Creating scoretrack for %s" %(c.geneid))
                upro=uniprotobj.get_protein_object_by_ensgid(c.geneid)
                scoretrack=[0]*len(upro)
                counttrack=[0]*len(upro)
                try:
                    scoretrack[int(c.protpos)-1]=float(c.rawscore)
                    counttrack[int(c.protpos)-1]=1
                except IndexError: #position out of bounds for score track, ignore it
                    pass
                entryname_scoretrack[entryname]=scoretrack
                entryname_scorecounts[entryname]=counttrack
            else:
                scoretrack=entryname_scoretrack[entryname]
                counttrack=entryname_scorecounts[entryname]
                try:
                    scoretrack[int(c.protpos)-1]=scoretrack[int(c.protpos)-1]+float(c.rawscore)
                    counttrack[int(c.protpos)-1]=counttrack[int(c.protpos)-1]+1
                except IndexError: #position out of bounds for score track, ignore it
                    pass
                entryname_scoretrack[entryname]=scoretrack
                entryname_scorecounts[entryname]=counttrack
        except KeyError:
            #print("%s ENSG ID not found in the uniprot proteome!" %(c.geneid))
            ensgsnotfound.append(c.geneid)
            pass
    '''
    entryname_meanscoretrack={}
    for e in entryname_scoretrack:
        ccountrack=entryname_scorecounts[e]
        sscoretrack=entryname_scoretrack[e]
        meanscoretrack=[0]*len(sscoretrack)
        for i in range(0,len(sscoretrack)):
            if counttrack[i]>0:
                meanscoretrack[i]=sscoretrack[i]/float(ccountrack[i])
        entryname_meanscoretrack[e]=meanscoretrack
    print(",".join(ensgsnotfound))
    print(len(ensgsnotfound))
    return entryname_meanscoretrack

def functional_site_scoring(uniprotobj):
    fsites=ProteinFunctionalSites(functionalsitefile)
    protein_functional_score_track={} #uniprot->[functional score track, binary scoring for each protein]
    for ename in uniprotobj.get_protein_entrynames():
        #print(ename)
        plen=uniprotobj.get_protein_length(ename)
        scoretrack=[0]*plen
        #protein_functional_score_track[ename]=[0]*plen #creates an empty track
        #condition wherein the tracks do not match
        #sites found
        if ename in fsites.protein_functionalsite_map_:
            sites=fsites.get_functional_sites(ename)
            for s in sites:
                startindex=s.feature_start_-1
                endindex=s.feature_end_-1
                for i in range(startindex, endindex+1):
                    try:
                        scoretrack[i]=1
                    except IndexError:
                        print("%d %d %s %s" %(startindex, endindex, s.feature_id_, ename))
            protein_functional_score_track[ename]=scoretrack
        #sites not found in protein
        else:
            protein_functional_score_track[ename]=scoretrack
    return protein_functional_score_track

def domain_conservation_scoring(uniprotobj):
    protein_domaincons_score_track={}
    pdr=ProteinDomainRecords(genedomainmap)
    print("domain records extracted")
    counter=0
    for ename in uniprotobj.get_protein_entrynames():
        plen=uniprotobj.get_protein_length(ename)
        print("processing %d %s " %( counter, ename))
        counter+=1
        scoretrack=[0]*plen
        #fetch accessions, check for domains in gene domain class and store the array - then extract scores for those
        accs=uniprotobj.get_accessions_of_protein(ename)
        for a in accs:
            if pdr.has_domains(a):
                domains=pdr.get_domains_for_protein(a)
                for dr in domains: #returns domain record
                    try:
                        dcons=DomainConservation(dr.domainid_, domalnscoredir+dr.domainid_+".aln.fa.score", domalndir+dr.domainid_+".aln.fa" )
                        #zscore=dcons.protein_position_residue_conservation_score_list(a)
                        zscore=dcons.domain_conservationscore_zscore_list(a)
                        z=0
                        try:
                            for i in range(dr.domainstart_-1, dr.domainend_):
                                scoretrack[i]=zscore[z]
                                z+=1
                        except IndexError as ie:
                            pass
                    except FileNotFoundError as fe:
                        pass
                        #print(ie)
                        #print("Error for domain: %s in protein: %s (start %d, end %d)" %(dr.domainid_, a, dr.domainstart_, dr.domainend_))
                break
        protein_domaincons_score_track[ename]=scoretrack
    return protein_domaincons_score_track
def mutational_statisticalsignificance_scoring(uniprotobj):
    protein_mutsig_score_track={}
    protein_analogsig_score_track={}
    phm=ProteinHotspotMutations(tcgastatsfile)
    #dev mode
    #krashot=phm.get_hotspot_mutations('KRAS')
    #for k in krashot:
    #    print("%s %s" %(k.position_, k.domain_id_))
    #sys.exit(0)
    pdomrec=ProteinDomainRecords(genedomainmap)
    counter=0
    statanmapfo=open("../data/statsig_analogpos.tsv",'w')
    sigtrackfo=open("../data/statsig_scores1.tsv",'w')
    anatrackfo=open("../data/analog_scores1.tsv",'w')
    statanmapfo.write("statsigpos\tanalogpos\n")
    entrykey_hotspots={}
    for ename in uniprotobj.get_protein_entrynames():
        eprot=uniprotobj.get_protein_object_by_entryname(ename)
        plen=len(eprot)
        #print("processing %d %s " %( counter, ename))
        counter+=1
        scoretrack=[0]*plen
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
            try:
                scoretrack[h.position_-1]=1
            except IndexError as ie:
                pass
        protein_mutsig_score_track[ename]=scoretrack
        protein_analogsig_score_track[ename]=[0]*plen
        #significant positions added
    print("SIgnificance position scores added")
    #Going for analogous positions
    #for each of the position if inside domain check for analogous positions
    for ek in entrykey_hotspots:
        ekhot=entrykey_hotspots[ek]
        sek=ek.split("_")
        enm=sek[0]+"_"+sek[1]
        eaccs=uniprotobj.get_accessions_of_protein(enm)
        ppos=int(sek[2])
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
            if len(analogmap) > 0:
                for an in analogmap:
                    if analogmap[an]>=0:
                        try: #get protein entry name using accession
                            enn=uniprotobj.get_protein_object_by_accession(an)
                            ennstrack=protein_analogsig_score_track[enn.entry_name_]
                            ennstrack[analogmap[an] - 1]=ennstrack[ analogmap[an] -1]+0.6
                            statanmapfo.write("%s\t%s\n" %(ek, enn.entry_name_+"_"+str(analogmap[an])))
                            protein_analogsig_score_track[enn.entry_name_]=ennstrack
                        except KeyError: #accession not found - ignore
                            pass
                        except IndexError: #domain analogous position goes beyond protein
                            pass
    statanmapfo.close()
    print_scoring_track(protein_mutsig_score_track, sigtrackfo)
    print_scoring_track(protein_analogsig_score_track, anatrackfo)
    sigtrackfo.close()
    anatrackfo.close()
    #Merge two tracks
    merged_mutsig_score_track={}
    for ee in protein_mutsig_score_track:
        domstrack=protein_analogsig_score_track[ee]
        sigtrack=protein_mutsig_score_track[ee]
        tottrack=[0]*len(sigtrack)
        for m in range(0,len(sigtrack)):
            tottrack[m]=sigtrack[m]+domstrack[m]
        merged_mutsig_score_track[ee]=tottrack
    return merged_mutsig_score_track

def main():
    #statout=open("../data/phosphositeplus_functional_scores.tsv",'w')
    statout=open("../data/statsig_analogous_scores5.tsv",'w')
    #domout=open("../data/domain_cons_zscores2.tsv",'w')
    #caddout=open("../data/caddmean_scores.tsv",'w')
    uprot=UniprotProteomeParser(proteomeswissprotfile)
    #caddsc=cadd_scoring(uprot)
    #for e in uprot.ensgid_entryname_:
    #    caddout.write("%s %s\n" %(e, uprot.ensgid_entryname_[e]))
    #caddout.close()
    print("Uniprot parsed")
    #print_scoring_track(caddsc, caddout)
    #caddout.close()
    #caddsc=cadd_scoring(uprot)
    mssc=mutational_statisticalsignificance_scoring(uprot)
    print_scoring_track(mssc, statout)
    #pspscores=phoshositeplus_functional_scoring(uprot)
    #print_scoring_track(pspscores, statout)
    statout.close()
    #pdoms=domain_conservation_scoring(uprot)
    #print(pdoms)
    #pfst=functional_site_scoring(uprot)
    #print_scoring_track(pdoms, domout)

if __name__=='__main__':
    main()
