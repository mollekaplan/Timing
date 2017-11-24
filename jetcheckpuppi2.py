#! /usr/bin/env python

import ROOT
from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
import os
import math

o2pi = 1./(2.*math.pi)

def deltaphi(a,b):
    x = b-a
    if abs(x) <= math.pi:
        return x
    
    n = round(x*o2pi)
    return x - n*2.*math.pi

def deltar(a,b):
    deta = b.eta()-a.eta()
    dphi = deltaphi(a.phi(),b.phi())
    return math.sqrt(deta*deta+dphi*dphi)


def matchgenjet(jet, genjets):
    drmin = float("inf")
    matched = None
    for genjet in genjets:
        if genjet.pt() < 4.0:
            continue
        dr = deltar(jet,genjet) 
        if dr<0.2 and dr<drmin:
        #if dr<drmin:
            drmin = dr
            matched = genjet
            
    return matched

def matchgenjetloose(jet, genjets):
    drmin = float("inf")
    matched = None
    for genjet in genjets:
        dr = deltar(jet,genjet) 
        if dr<0.4 and dr<drmin:
        #if dr<drmin:
            drmin = dr
            matched = genjet
            
    return matched

def ispujet(jet,genjets):
    for genjet in genjets:
        if genjet.pt() < 4.0:
            continue
        dr = deltar(jet,genjet) 
        if dr<0.6:
            return False
    return True

#inputfiles = ['root://eoscms//eos/cms/store/group/upgrade/timing/pfintegration/Jan24jme/dymm200timing/MINIAODSIM/step3_dymm200timing_MINIAODSIM_55_1.root']

#inputfiles = ['root://eoscms//eos/cms/store/group/upgrade/timing/pfintegration/Jan24jme/dymm200timing/AODSIM/step3_dymm200timing_AODSIM_71_1.root']
#inputfiles = ['file:/eos/cms/store/group/upgrade/timing/pfintegration/Jan24jme/dymm200timing/AODSIM/step3_dymm200timing_AODSIM_71_1.root']

#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Jan24jme/dymm200notiming/AODSIM/'
#indir = '/eos/cms/store/relval/CMSSW_9_0_0_pre2/RelValZMM_14/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v1_2023D4TimingPU200-v1/10000/'

nfiles = -1

#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore3_alldensitylocal/dymm200timinggeom/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore2_alldensity/dymm200timinggeom/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore_v2/dymm200timinggeom/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore_v2/dymm200timinggeom/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod2_puppimod/dymm200timing/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar29jme_reminiaod_recombeta/dymm200timingmaxeta15/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar30jme/dymm200timingmaxeta15/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme/dymm200timingmaxeta15/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme/dymm200timing/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme/dymm200timingmaxeta100/MINIAODSIM/'
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme_reminiaod_vtxassign/dymm200timing/MINIAODSIM/'
#indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/PU200_timing/171110_104722/0000/'
#indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/NoPU_timing_/171109_155602/0000/'
#indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1sigma_PU200_timing/171116_091208/0000/'
#indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_PU200_timing/171116_095856/0000/'
#indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1sigma_NoPU_timing/171116_091101/0000/'
#indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_NoPU_timing/171116_095748/0000/'
indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/NoPU_timing_/171109_155602/0000/'
inputfiles = [indir + filename for filename in os.listdir(indir)[0:nfiles]]

#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore3_alldensitylocal/dymm200notiminggeom/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore2_alldensity/dymm200notiminggeom/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore_v2/dymm200notiminggeom/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod_puppimodmore_v2/dymm200notiminggeom/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme_reminiaod2_puppimod/dymm200timing/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar22jme/dymm200notiming/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar29jme_reminiaod_recombeta/dymm200notiming/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme/dymm200notiming/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme/dymm200timing/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme_reminiaod_vtxassign/dymm200notiming/MINIAODSIM/'
#indir2 = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/PU200_notiming/171110_104953/0000/'
#indir2 = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/NoPU_notiming_/171109_155445/0000/'
#indir2 = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1sigma_PU200_notiming/171116_091300/0000/'
#indir2 = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_PU200_notiming/171116_095955/0000/'
#indir2 = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_NoPU_notiming/171116_095645/0000/'
indir2 = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/NoPU_notiming_/171109_155445/0000/'
inputfiles2 = [indir2 + filename for filename in os.listdir(indir2)[0:nfiles]]



#print(inputfiles)


# Make VarParsing object
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile#VarParsing_Example
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

# Events takes either
# - single file name
# - list of file names
# - VarParsing options

# use Varparsing object
#events = Events (options)

events = Events (inputFiles=inputfiles)
events2 = Events (inputFiles=inputfiles2)
#events3 = Events (inputFiles=inputfiles3)



#vector<reco::Vertex>                  "inclusiveSecondaryVertices"   ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVertices"    ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVertices1D"   ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVertices1DWithBS"   ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVerticesLegacy"   ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVerticesLegacyWithBS"   ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVerticesWithBS"   ""                "RECO"    


# create handle outside of loop
handlevtx  = Handle ("std::vector<reco::Vertex>")
labelvtx = ("offlineSlimmedPrimaryVertices")

handlevtxalt  = Handle ("std::vector<reco::Vertex>")
labelvtxalt = ("offlinePrimaryVerticesLegacy")

handlexyz0 = Handle("ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>")
labelxyz0 = ("genParticles","xyz0")

handlet0 = Handle("float")
labelt0 = ("genParticles","t0")

handlepfcand = Handle("vector<pat::PackedCandidate>")
labelpfcand = ("packedPFCandidates")

handlejet = Handle("std::vector<reco::PFJet>")
labeljet = "AK4PFJets"

handlejetchs = Handle("std::vector<pat::Jet>")
labeljetchs = "slimmedJets"

handlejetpuppi = Handle("std::vector<pat::Jet>")
labeljetpuppi = "slimmedJetsPuppi"

handlegenjet = Handle("std::vector<reco::GenJet>")
labelgenjet = "slimmedGenJets"

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
#label = ("offlineSlimmedPrimaryVertices")


# Create histograms, etc.
#ROOT.gROOT.SetBatch()        # don't pop up canvases
#ROOT.gROOT.SetStyle('Plain') # white background

jetetachstiming = ROOT.TH1D("jetetachstiming","",50,-5.,5.)
jetetachsnotiming = ROOT.TH1D("jetetachsnotiming","",50,-5.,5.)
jetetachsnotiming.SetLineColor(ROOT.kRed)

jetetapuppitiming = ROOT.TH1D("jetetapuppitiming","",50,-5.,5.)
jetetapuppinotiming = ROOT.TH1D("jetetapuppinotiming","",50,-5.,5.)
jetetapuppinotiming.SetLineColor(ROOT.kRed)

matchedjetetapuppitiming = ROOT.TH1D("matchedjetetapuppitiming","",50,-5.,5.)
matchedjetetapuppinotiming = ROOT.TH1D("matchedjetetapuppinotiming","",50,-5.,5.)
matchedjetetapuppinotiming.SetLineColor(ROOT.kRed)

unmatchedjetetapuppitiming = ROOT.TH1D("unmatchedjetetapuppitiming","",50,-3.,3.)
unmatchedjetetapuppinotiming = ROOT.TH1D("unmatchedjetetapuppinotiming","",50,-3.,3.)
unmatchedjetetapuppinotiming.SetLineColor(ROOT.kRed)

unmatchedjetptpuppitiming = ROOT.TH1D("unmatchedjetptpuppitiming","",60,-0.,60.)
unmatchedjetptpuppinotiming = ROOT.TH1D("unmatchedjetptpuppinotiming","",60,-0.,60.)
unmatchedjetptpuppinotiming.SetLineColor(ROOT.kRed)

unmatchedjetmetapuppitiming = ROOT.TH1D("unmatchedjetmetapuppitiming","",50,-5.,5.)
unmatchedjetmetapuppinotiming = ROOT.TH1D("unmatchedjetmetapuppinotiming","",50,-5.,5.)
unmatchedjetmetapuppinotiming.SetLineColor(ROOT.kRed)

candetatiming = ROOT.TH1D("candetatiming","",100,-5.,5.)
candetanotiming = ROOT.TH1D("candetanotiming","",100,-5.,5.)
candetanotiming.SetLineColor(ROOT.kRed)

candetatimingsimple = ROOT.TH1D("candetatimingsimple","",100,-5.,5.)
candetanotimingsimple = ROOT.TH1D("candetanotimingsimple","",100,-5.,5.)
candetanotimingsimple.SetLineColor(ROOT.kRed)

candetapuppitiming = ROOT.TH1D("candetapuppitiming","",100,-5.,5.)
candetapuppinotiming = ROOT.TH1D("candetapuppinotiming","",100,-5.,5.)
candetapuppinotiming.SetLineColor(ROOT.kRed)

jetresponsetiming = ROOT.TH1D("jetresponsetiming","",100,0.,2.)
jetresponsenotiming = ROOT.TH1D("jetresponsenotiming","",100,0.,2.)
jetresponsenotiming.SetLineColor(ROOT.kRed)

notmatchedjetetapuppitiming = ROOT.TH1D("notmatchedjetetapuppitiming","",50,-3.,3.)
notmatchedjetetapuppinotiming = ROOT.TH1D("notmatchedjetetapuppinotiming","",50,-3.,3.)
notmatchedjetetapuppinotiming.SetLineColor(ROOT.kRed)


genjettiming = ROOT.TH2D("genjettiming","",60,0.,60.,50,-5.,5.)
genjetnotiming = ROOT.TH2D("genjetnotiming","",60,0.,60.,50,-5.,5.)
genjetnotiming.SetLineColor(ROOT.kRed)

matchtiming = ROOT.TH2D("matchtiming","",60,0.,60.,50,-5.,5.)
matchnotiming = ROOT.TH2D("matchnotiming","",60,0.,60.,50,-5.,5.)
matchnotiming.SetLineColor(ROOT.kRed)


#Create fake versus eta data set
#fake = RooRealVar("fake","",0,2)
#faketiming = RooDataSet("faketiming","",RooArgSet(eta,fake))
#fakenotiming = RooDataSet("fakenotiming","",RooArgSet(eta,fake))


maxevt = 10000

minjetpt = 30.

# loop over events
ievt=0
for event in events:

    event.getByLabel (labelvtx, handlevtx)
    event.getByLabel (labelxyz0, handlexyz0)
    event.getByLabel (labelt0, handlet0)
    event.getByLabel (labelpfcand, handlepfcand)
    #event.getByLabel (labeljet, handlejet)
    event.getByLabel (labeljetchs, handlejetchs)
    event.getByLabel (labeljetpuppi, handlejetpuppi)
    event.getByLabel (labelgenjet, handlegenjet)
    
    # get the product
    vtxs = handlevtx.product()
    #jets = handlejet.product()
    chsjets = handlejetchs.product()
    puppijets = handlejetpuppi.product()
    genjets = handlegenjet.product()
    t0 = handlet0.product()
    pfcands = handlepfcand.product()
 
    pvtx = vtxs[0]
            
    nmuons = 0
    for cand in pfcands:
        if abs(cand.pdgId()) != 13:
            continue
        if cand.pt()<20.:
            continue
        if not (cand.pvAssociationQuality()>=1 and cand.vertexRef().key()==0):
            continue
        nmuons += 1
        
    if nmuons < 2:
        continue
            
    mugenjets = []
    #print("first genjet loop")
    for genjet in genjets:
	genjettiming.Fill(genjet.pt(),genjet.eta())
	for genpart in genjet.getJetConstituentsQuick():
            #print("genpart")
            #genpart = genjet.getGenConstituent(ipart)
            if abs(genpart.pdgId())==13 and genpart.fromHardProcessFinalState():
                #printf("match")
                mugenjets.append(genjet)
                break
        if len(mugenjets)==2:
            break
        
    for jet in chsjets:
        if jet.pt()>minjetpt:
            mugenjet = matchgenjetloose(jet,mugenjets)
            if mugenjet:
                continue
              
            jetetachstiming.Fill(jet.eta())
              

    for jet in puppijets:
	#print(jet.jecFactor("Uncorrected")/jet.jecFactor("L2Relative"))
      
        if jet.pt()>minjetpt:
            mugenjet = matchgenjetloose(jet,mugenjets)
            if mugenjet:
                continue

          
            genjet = matchgenjet(jet,genjets)
            genjet2 = matchgenjetloose(jet,genjets)
	    ispu = ispujet(jet,genjets)
            
            if not genjet in mugenjets:
                jetetapuppitiming.Fill(jet.eta())
                            
            if genjet:
		matchtiming.Fill(jet.pt(),jet.eta())

	    if genjet:
		matchedjetetapuppitiming.Fill(jet.eta())
                
	    if genjet and abs(jet.eta())<3.:
                    jetresponsetiming.Fill(jet.energy()/genjet.energy())
            
	    if not genjet:
		notmatchedjetetapuppitiming.Fill(jet.eta())

	    if ispu:
                unmatchedjetetapuppitiming.Fill(jet.eta())
                unmatchedjetptpuppitiming.Fill(jet.pt())
                absetamin = float("inf")
                meta = -5.
                for cand in jet.getJetConstituentsQuick():
                    if abs(cand.charge())>0.:
                        abseta = abs(cand.eta())
                        if abseta<absetamin:
                            absetamin = abseta
                            meta = cand.eta()
                unmatchedjetmetapuppitiming.Fill(meta)

    for cand in pfcands:
        #if abs(cand.charge())>0. and cand.pt()>0.7:
        if cand.charge()==0.:
            candetapuppitiming.Fill(cand.eta(),cand.puppiWeight())
        if abs(cand.charge())>0.:
            if cand.pvAssociationQuality()>=1 and cand.vertexRef().key()==0:
                candetatiming.Fill(cand.eta())
            #usetime = cand.timeError()>0. and pvtx.tError()>0.
	    #usetime = pvtx.tError()>0.
	    usetime = False
            #if abs(cand.dz(0))<0.1 and (not usetime or abs(cand.dtime(0))<0.09):
	    #if abs(cand.dz(0))<0.1 and (not usetime):
	    if abs(cand.dz(0))<0.1:
                candetatimingsimple.Fill(cand.eta())

 
    ievt += 1
    #if ievt==maxevt:
    #    break;

print(ievt)

jetetachstiming.Scale(1./float(ievt))
jetetapuppitiming.Scale(1./float(ievt))
matchedjetetapuppitiming.Scale(1./float(ievt))
unmatchedjetetapuppitiming.Scale(1./float(ievt))
unmatchedjetptpuppitiming.Scale(1./float(ievt))
unmatchedjetmetapuppitiming.Scale(1./float(ievt))
candetatiming.Scale(1./float(ievt))
candetatimingsimple.Scale(1./float(ievt))
candetapuppitiming.Scale(1./float(ievt))
jetresponsetiming.Scale(1./float(ievt))
matchtiming.Scale(1./float(ievt))
genjettiming.Scale(1./float(ievt))
notmatchedjetetapuppitiming.Scale(1./float(ievt))


ievt = 0
for event in events2:

    event.getByLabel (labelvtx, handlevtx)
    event.getByLabel (labelxyz0, handlexyz0)
    event.getByLabel (labelt0, handlet0)
    event.getByLabel (labelpfcand, handlepfcand)
    #event.getByLabel (labeljet, handlejet)
    event.getByLabel (labeljetchs, handlejetchs)
    event.getByLabel (labeljetpuppi, handlejetpuppi)
    event.getByLabel (labelgenjet, handlegenjet)
    
    # get the product
    vtxs = handlevtx.product()
    #jets = handlejet.product()
    chsjets = handlejetchs.product()
    puppijets = handlejetpuppi.product()
    genjets = handlegenjet.product()
    t0 = handlet0.product()
    pfcands = handlepfcand.product()
 
    pvtx = vtxs[0]
            
    nmuons = 0
    for cand in pfcands:
        if abs(cand.pdgId()) != 13:
            continue
        if cand.pt()<20.:
            continue
        if not (cand.pvAssociationQuality()>=1 and cand.vertexRef().key()==0):
            continue
        nmuons += 1
        
    if nmuons < 2:
        continue
            
    mugenjets = []
    #print("first genjet loop")
    for genjet in genjets:
	genjetnotiming.Fill(genjet.pt(),genjet.eta())
	for genpart in genjet.getJetConstituentsQuick():
            #print("genpart")
            #genpart = genjet.getGenConstituent(ipart)
            if abs(genpart.pdgId())==13 and genpart.fromHardProcessFinalState():
                #printf("match")
                mugenjets.append(genjet)
                break
        if len(mugenjets)==2:
            break
              
    for jet in chsjets:
        if jet.pt()>minjetpt:
            mugenjet = matchgenjetloose(jet,mugenjets)
            if mugenjet:
                continue
              
            jetetachsnotiming.Fill(jet.eta())
              
    for jet in puppijets:
        if jet.pt()>minjetpt:
            mugenjet = matchgenjetloose(jet,mugenjets)
            if mugenjet:
                continue

          
            genjet = matchgenjet(jet,genjets)
            genjet2 = matchgenjetloose(jet,genjets)
	    ispu = ispujet(jet,genjets)
            
            if not genjet in mugenjets:
                jetetapuppinotiming.Fill(jet.eta())
                            
	    if genjet:
		matchnotiming.Fill(jet.pt(),jet.eta())

            if genjet:
                matchedjetetapuppinotiming.Fill(jet.eta())
                if abs(jet.eta())<3.:
                    jetresponsenotiming.Fill(jet.energy()/genjet.energy())
          
            if not genjet:
                notmatchedjetetapuppinotiming.Fill(jet.eta())

	    if ispu:
                unmatchedjetetapuppinotiming.Fill(jet.eta())
                unmatchedjetptpuppinotiming.Fill(jet.pt())
                absetamin = float("inf")
                meta = -5.
                for cand in jet.getJetConstituentsQuick():
                    if abs(cand.charge())>0.:
                        abseta = abs(cand.eta())
                        if abseta<absetamin:
                            absetamin = abseta
                            meta = cand.eta()
                unmatchedjetmetapuppinotiming.Fill(meta)
                  
    for cand in pfcands:
        #if abs(cand.charge())>0. and cand.pt()>0.7:
        if cand.charge()==0.:
            candetapuppinotiming.Fill(cand.eta(),cand.puppiWeight())
        if abs(cand.charge())>0.:
            if cand.pvAssociationQuality()>=1 and cand.vertexRef().key()==0:
                candetanotiming.Fill(cand.eta())
            #usetime = cand.timeError()>0. and pvtx.tError()>0.
            #usetime = pvtx.tError()>0.
            usetime = False
            #if abs(cand.dz(0))<0.1 and (not usetime or abs(cand.dtime(0))<0.09):
            #if abs(cand.dz(0))<0.1 and (not usetime):
            if abs(cand.dz(0))<0.1:
                candetatimingsimple.Fill(cand.eta())

                
    ievt += 1
    if ievt==maxevt:
        break;
   

jetetachsnotiming.Scale(1./float(ievt))
jetetapuppinotiming.Scale(1./float(ievt))
matchedjetetapuppinotiming.Scale(1./float(ievt))
unmatchedjetetapuppinotiming.Scale(1./float(ievt))
unmatchedjetptpuppinotiming.Scale(1./float(ievt))
unmatchedjetmetapuppinotiming.Scale(1./float(ievt))
candetanotiming.Scale(1./float(ievt))
candetanotimingsimple.Scale(1./float(ievt))
candetapuppinotiming.Scale(1./float(ievt))
jetresponsenotiming.Scale(1./float(ievt))
matchnotiming.Scale(1./float(ievt))
genjetnotiming.Scale(1./float(ievt))
notmatchedjetetapuppinotiming.Scale(1./float(ievt))

#c1 = ROOT.TCanvas()
#jetetachsnotiming.Draw("HIST")
#jetetachstiming.Draw("HISTSAME")
#c1.SaveAs("jetetachs.pdf")

#c2 = ROOT.TCanvas()
#jetetapuppinotiming.Draw("HIST")
#jetetapuppitiming.Draw("HISTSAME")
#c2.SaveAs("jetetapuppi.pdf")    

c3 = ROOT.TCanvas()
matchedjetetapuppinotiming.Draw("HIST")
matchedjetetapuppitiming.Draw("HISTSAME")    

leg3 = ROOT.TLegend(0.4,.8,.6,.9)
leg3.AddEntry(matchedjetetapuppitiming,"timing","L")
leg3.AddEntry(matchedjetetapuppinotiming,"no-timing","L")
leg3.Draw()

c3.SaveAs("matchedjetetapuppi.pdf")
c3.SaveAs("matchedjetetapuppi.root")

c4 = ROOT.TCanvas()
unmatchedjetetapuppinotiming.Draw("HIST")
unmatchedjetetapuppitiming.Draw("HISTSAME")

leg4 = ROOT.TLegend(0.4,.8,.6,.9)
leg4.AddEntry(unmatchedjetetapuppitiming,"timing","L")
leg4.AddEntry(unmatchedjetetapuppinotiming,"no-timing","L")
leg4.Draw()

c4.SaveAs("unmatchedjetetapuppi.pdf")
c4.SaveAs("unmatchedjetetapuppi.root")

#c4a = ROOT.TCanvas()
#unmatchedjetptpuppinotiming.Draw("HIST")
#unmatchedjetptpuppitiming.Draw("HISTSAME")
#c4a.SaveAs("unmatchedjetptpuppi.pdf")
#c4a.SaveAs("unmatchedjetptpuppi.root")

#c4a = ROOT.TCanvas()
#unmatchedjetmetapuppinotiming.Draw("HIST")
#unmatchedjetmetapuppitiming.Draw("HISTSAME")

#matchedpuppiratio = matchedjetetapuppitiming.Clone()
#matchedpuppiratio.Divide(matchedjetetapuppinotiming)

#c5 = ROOT.TCanvas()
#matchedpuppiratio.Draw("HIST")
#c5.SaveAs("matchedpuppiratio.pdf")
#c5.SaveAs("matchedpuppiratio.root")

#unmatchedpuppiratio = unmatchedjetetapuppitiming.Clone()
#unmatchedpuppiratio.Divide(unmatchedjetetapuppinotiming)

#c6 = ROOT.TCanvas()
#unmatchedpuppiratio.Draw()
#c6.SaveAs("unmatchedpuppiratio.pdf")
#c6.SaveAs("unmatchedpuppiratio.root")

#c7 = ROOT.TCanvas()
#candetanotiming.Draw("HIST")
#candetatiming.Draw("HISTSAME")
#c7.SaveAs("candeta.pdf")
#c7.SaveAs("candeta.root")

#c7a = ROOT.TCanvas()
#candetanotimingsimple.Draw("HIST")
#candetatimingsimple.Draw("HISTSAME")
#c7a.SaveAs("candetasimple.pdf")
#c7a.SaveAs("candetasimple.root")

#c7b = ROOT.TCanvas()
#candetanotimingsimplemod = candetanotimingsimple.Clone()
#candetanotimingsimplemod.SetLineColor(ROOT.kBlue)
#candetanotiming.Draw("HIST")
#candetanotimingsimplemod.Draw("HISTSAME")

#c7d = ROOT.TCanvas()
#candetapuppinotiming.Draw("HIST")
#candetapuppitiming.Draw("HISTSAME")
#c7d.SaveAs("candetapuppi.pdf")
#c7d.SaveAs("candetapuppi.root")

#c8 = ROOT.TCanvas()
#jetresponsetiming.Draw("HIST")
#jetresponsenotiming.Draw("HISTSAME")
#c8.SaveAs("jetresponse.pdf")
#c8.SaveAs("jetresponse.root")

#c9 = ROOT.TCanvas()
#genjettiming.Draw("HIST")
#genjetnotiming.Draw("HISTSAME")
#c9.SaveAs("denom.pdf")
#c9.SaveAs("denom.root")

#c10 = ROOT.TCanvas()
#matchtiming.Draw("HIST")
#matchnotiming.Draw("HISTSAME")
#c10.SaveAs("num.pdf")
#c10.SaveAs("num.root")

#efficiencytiming = matchtiming.Clone()
#efficiencytiming.Divide(genjettiming)

#efficiencynotiming = matchnotiming.Clone()
#efficiencynotiming.Divide(genjetnotiming)

#c11 = ROOT.TCanvas()
#efficiencytiming.Draw("HIST")
#efficiencynotiming.Draw("HISTSAME")
#c11.SaveAs("efficiency.pdf")
#c11.SaveAs("efficiency.root")

#c12 = ROOT.TCanvas()
#notmatchedjetetapuppinotiming.Draw("HIST")
#notmatchedjetetapuppitiming.Draw("HISTSAME")
#c12.SaveAs("notmatched.pdf")
#c12.SaveAs("notmatched.root")

#notmatchedratio = notmatchedjetetapuppitiming.Clone()
#notmatchedratio.Divide(notmatchedjetetapuppinotiming)

#c13 = ROOT.TCanvas()
#notmatchedratio.Draw()
#c13.SaveAs("notmatchedratio.pdf")
#c13.SaveAs("notmatchedratio.root")

input("Press Enter to continue...")
