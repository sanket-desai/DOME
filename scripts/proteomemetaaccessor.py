import os
import sys
import pysam
import pandas as pd
import numpy as np
from aacids import *
import networkx as nx
from metaaccessor import *


def main():
    infi=open(sys.argv[1])
    ofo=open(sys.argv[2],'w')
    #infi=open("j")
    entry_queries={}
    pindex=1
    for i in infi:
        si=i.strip().split(",")
        for l in list(range(1,int(si[1]))):
            qsi=Query(si[0], l)
            if si[0] in entry_queries:
                temp=entry_queries[si[0]]
                temp.append(qsi)
                entry_queries[si[0]]=temp
            else:
                entry_queries[si[0]]=[ qsi ]
    refpath="../data/domerefdata/final/"
    oheader=["tier","protein","mutpos","refaa","altaa", "domain", "statscore","analogtomut", "caddscore", "cosmiccount", "colcoscount", "colnumofresidues" , "clinvarsig", "uniprotsite", "phosphosite", "uniprotmutagen" , "hotspot3dann", "3d_dom_close", "3d_statsig_close" , "3d_analog_close","3d_uniprotsite_close", "3d_pspsite_close", "3d_hotspot_close", "interfaceposition"]
    refdbobj=ReferenceDatabases(refpath)
    ofo.write( "%s\n" %("\t".join(oheader)))
    for ent in entry_queries:
        structnotfound=False
        print("%d\t%s" %(pindex, ent))
        pindex+=1
        try:
            entnet=ProteinNetwork(ent, refdbobj) #Load protein network
        except ProteinNetworkException:
            structnotfound=True
        #To be duplicated with metaaccessor
        upfamrec=refdbobj.uniprotpfamdbh_.get_pfam_domain_records_for_protein_entry(ent)
        domname_alignments={} #stores domain alignments for a particular domain name (in this case for a protein)
        for upfa in upfamrec:
            try:
                domname_alignments[upfa.domainname_]= DomainMultiAlignedSequence(upfa.domainid_, refdbobj.domaligndir_+upfa.domainid_+".aln.fa")
            except FileNotFoundError as fnferr:
                pass
        if not structnotfound:
            for q in entry_queries[ent]:
                outvec=[]
                out3dvec=[]
                qneighbors=entnet.get_neighbor_positions(q) #return int positions
                outvec.append(q.entryname_)
                outvec.append(str(q.position_))
                outvec.append(refdbobj.caddscoresdbh_.get_reference_aminoacid(q)) #loding sequeces might help
                outvec.append(q.altaa_)
                totcoscount=0
                colresidues=0
                if refdbobj.uniprotpfamdbh_.is_within_domain(q):
                    domname=refdbobj.uniprotpfamdbh_.get_pfam_domain_name(q)
                    domaln=[]
                    try:
                        domaln=domname_alignments[domname]
                    except FileNotFoundError:
                        pass
                    except KeyError:
                        pass
                    daccs=refdbobj.uniprotpfamdbh_.get_pfam_domain_accessions(q)
                    if len(daccs) > 1:
                        #print(daccs)
                        #cant handle this currently
                        daccs=list(np.unique(np.array(daccs)))
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
                    sscoresstr=refdbobj.statsigs_.get_significance_scores(q)
                    sstumours=refdbobj.statsigs_.get_tumour_types(q)
                    outvec.append(";".join(sscoresstr))
                    #outvec.append(";".join(sscoresstr)+"|"+';'.join(sstumours))
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
                #except ValueError:
                #    print(q.entryname_)
                #    print(q.position_)
                #    sys.exit(0)
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
                #print(q.entryname_)
                #print(q.position_)
                #print(refdbobj.proteininterfaces_)
                #sys.exit(0)
                #print(refdbobj.proteininterfaces_.get_interface_positions_by_entryname(q.entryname_))
                if q.position_ in refdbobj.proteininterfaces_.get_interface_positions_by_entryname(q.entryname_):
                    out3dvec.append("I") #Is an interface
                else:
                    out3dvec.append('-')
                dorec="\t".join(outvec) + "\t" + "\t".join(out3dvec)
                dr=DOMEOutputRecordParser(dorec)
                ofo.write("%d\t%s\n" %(dr.get_tier(), dorec))
                #dr=DOMEOutputRecordParser(dorec, statind, analogind, caddind, geostatind, geoanalogind, sitesind)
                #print("%s\t%s" %("\t".join(outvec), "\t".join(out3dvec) ) )
                #Needs to be uncommented
                #print( outvec )
                #print( out3dvec )
    infi.close()
    ofo.close()

if __name__=="__main__":
    main()
