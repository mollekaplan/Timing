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
        if dr<0.4 and dr<drmin:
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

def matchjetloose(genjet, jets):
    drmin = float("inf")
    matched = None
    for jet in jets:
        dr = deltar(genjet,jet)
        if dr<0.4 and dr<drmin:
            #if dr<drmin:
            drmin = dr
            matched = jet
            
    return matched

def matchjet(genjet, jets):
    drmin = float("inf")
    matched = None
    for jet in jets:
        if jet.pt() < 4.0:
            continue
	dr = deltar(genjet,jet)
        if dr<0.4 and dr<drmin:
            #if dr<drmin:
            drmin = dr
            matched = jet

    return matched

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
indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Aug26jme/dymmnoputiming/MINIAODSIM/'
inputfiles = [indir + filename for filename in os.listdir(indir)[0:nfiles]]
#for filename in os.listdir(indir)[0:nfiles]:
#    if filename != indir+'step3_dymm200timing_MINIAODSIM_119.root':
#	inputfiles = indir + filename


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
indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Aug26jme/dymmnopunotiming/MINIAODSIM/'
inputfiles2 = [indir2 + filename for filename in os.listdir(indir2)[0:nfiles]]
#for filename in os.listdir(indir2)[0:nfiles]:
#    if filename != indir2+'step3_dymm200notiming_MINIAODSIM_38.root':
#        inputfiles2 = indir2 + filename


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
handlepfcand = Handle("vector<pat::PackedCandidate>")
labelpfcand = ("packedPFCandidates")

handlejetpuppi = Handle("std::vector<pat::Jet>")
labeljetpuppi = "slimmedJetsPuppi"

handlechsjet = Handle("std::vector<pat::Jet>")
labelchsjet = "slimmedJets"

handlejetAK8 = Handle("std::vector<pat::Jet>")
labeljetAK8 = "slimmedJetsAK8"
#labeljetAK8 = ("slimmedJetsAK8PFPuppiSoftDropPacked","SubJets")

handlegenjet = Handle("std::vector<reco::GenJet>")
labelgenjet = "slimmedGenJets"

handlegenjetAK8 = Handle("std::vector<reco::GenJet>")
labelgenjetAK8 = "slimmedGenJetsAK8"

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
#label = ("offlineSlimmedPrimaryVertices")


# Create histograms, etc.
#ROOT.gROOT.SetBatch()        # don't pop up canvases
#ROOT.gROOT.SetStyle('Plain') # white background
#ROOT.gStyle.SetOptStat(0);

npuppijetstiming = ROOT.TH1D("npuppijetstiming","Number of Puppi Jets, No Pileup",20,0,20)
npuppijetsnotiming = ROOT.TH1D("npuppijetsnotiming","Number of Puppi Jets, No Pileup",20,0,20)
npuppijetsnotiming.SetLineColor(ROOT.kRed)

nchsjetstiming = ROOT.TH1D("nchsjetstiming","Number of CHS Jets, No Pileup",50,0,100)
nchsjetsnotiming = ROOT.TH1D("nchsjetsnotiming","Number of CHS Jets, No Pileup",50,0,100)
nchsjetsnotiming.SetLineColor(ROOT.kRed)

nAK8jetstiming = ROOT.TH1D("nAK8jetstiming","Number of AK8 Jets, No Pileup",10,0,10)
nAK8jetsnotiming = ROOT.TH1D("nAK8jetsnotiming","Number of AK8 Jets, No Pileup",10,0,10)
nAK8jetsnotiming.SetLineColor(ROOT.kRed)

npuppijetstiming2 = ROOT.TH1D("npuppijetstiming2","Number of Puppi Jets, No Pileup",20,0,20)
npuppijetsnotiming2 = ROOT.TH1D("npuppijetsnotiming2","Number of Puppi Jets, No Pileup",20,0,20)
npuppijetsnotiming2.SetLineColor(ROOT.kRed)

nchsjetstiming2 = ROOT.TH1D("nchsjetstiming2","Number of CHS Jets, No Pileup",50,0,100)
nchsjetsnotiming2 = ROOT.TH1D("nchsjetsnotiming2","Number of CHS Jets, No Pileup",50,0,100)
nchsjetsnotiming2.SetLineColor(ROOT.kRed)

nAK8jetstiming2 = ROOT.TH1D("nAK8jetstiming2","Number of AK8 Jets, No Pileup",10,0,10)
nAK8jetsnotiming2 = ROOT.TH1D("nAK8jetsnotiming2","Number of AK8 Jets, No Pileup",10,0,10)
nAK8jetsnotiming2.SetLineColor(ROOT.kRed)

maxevt = 5000

minjetpt = 30.

################ TIMING #################
ievt=0
for event in events:

    event.getByLabel (labelpfcand, handlepfcand)
    event.getByLabel (labeljetpuppi, handlejetpuppi)
    event.getByLabel (labelchsjet, handlechsjet)
    event.getByLabel (labelgenjet, handlegenjet)
    event.getByLabel (labeljetAK8, handlejetAK8)    
    event.getByLabel (labelgenjetAK8, handlegenjetAK8)	

    # get the product
    puppijets = handlejetpuppi.product()
    chsjets = handlechsjet.product()
    genjets = handlegenjet.product()
    pfcands = handlepfcand.product()
    jetsAK8 = handlejetAK8.product()
    genjetsAK8 = handlegenjetAK8.product() 
            
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
        for genpart in genjet.getJetConstituentsQuick():
            #print("genpart")
            #genpart = genjet.getGenConstituent(ipart)
            if abs(genpart.pdgId())==13 and genpart.fromHardProcessFinalState():
                #printf("match")
                mugenjets.append(genjet)
                break
        if len(mugenjets)==2:
            break


    npuppi=0
    npuppi2=0
    for jet in puppijets:
        mugenjet = matchgenjetloose(jet,mugenjets)
	if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            npuppi2 += 1
            if abs(jet.eta()) < 2.5:
                npuppi += 1
    npuppijetstiming.Fill(npuppi)
    npuppijetstiming2.Fill(npuppi2)

    nchs=0
    nchs2=0
    for jet in chsjets:
        mugenjet = matchgenjetloose(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
	    nchs2 += 1
            if abs(jet.eta()) < 2.5:
                nchs += 1
    nchsjetstiming.Fill(nchs)
    nchsjetstiming2.Fill(nchs2)

    nAK8=0
    nAK82=0
    for jet in jetsAK8:
        mugenjet = matchgenjetloose(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            nAK82 += 1
            if abs(jet.eta()) < 2.5:
                nAK8 += 1
    nAK8jetstiming.Fill(nAK8)
    nAK8jetstiming2.Fill(nAK82)

    ievt += 1
    if ievt==maxevt:
        break;

#npuppijetstiming.Scale(1./float(ievt))
#npuppijetstiming2.Scale(1./float(ievt))
#nchsjetstiming.Scale(1./float(ievt))
#nchsjetstiming2.Scale(1./float(ievt))
#nAK8jetstiming.Scale(1./float(ievt))
#nAK8jetstiming2.Scale(1./float(ievt))

################## NO TIMING #####################

ievt = 0
for event in events2:

    event.getByLabel (labelpfcand, handlepfcand)
    event.getByLabel (labeljetpuppi, handlejetpuppi)
    event.getByLabel (labelgenjet, handlegenjet)
    event.getByLabel (labeljetAK8, handlejetAK8)
    event.getByLabel (labelgenjetAK8, handlegenjetAK8)
    event.getByLabel (labelchsjet, handlechsjet)

    # get the product
    puppijets = handlejetpuppi.product()
    chsjets = handlechsjet.product()
    genjets = handlegenjet.product()
    pfcands = handlepfcand.product()
    jetsAK8 = handlejetAK8.product()
    genjetsAK8 = handlegenjetAK8.product()

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
        for genpart in genjet.getJetConstituentsQuick():
            #print("genpart")
            #genpart = genjet.getGenConstituent(ipart)
            if abs(genpart.pdgId())==13 and genpart.fromHardProcessFinalState():
                #printf("match")
                mugenjets.append(genjet)
                break
        if len(mugenjets)==2:
            break

    
    npuppi=0
    npuppi2=0
    for jet in puppijets:
        mugenjet = matchgenjetloose(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            npuppi2 += 1
            if abs(jet.eta()) < 2.5:
                npuppi += 1
    npuppijetsnotiming.Fill(npuppi)
    npuppijetsnotiming2.Fill(npuppi2)

    nchs=0
    nchs2=0
    for jet in chsjets:
        mugenjet = matchgenjetloose(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            nchs2 += 1
            if abs(jet.eta()) < 2.5:
                nchs += 1
    nchsjetsnotiming.Fill(nchs)
    nchsjetsnotiming2.Fill(nchs2)
    
    nAK8=0
    nAK82=0
    for jet in jetsAK8:
        mugenjet = matchgenjetloose(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            nAK82 += 1
            if abs(jet.eta()) < 2.5:
                nAK8 += 1
    nAK8jetsnotiming.Fill(nAK8)
    nAK8jetsnotiming2.Fill(nAK82)



    ievt += 1
    if ievt==maxevt:
        break;


#npuppijetsnotiming.Scale(1./float(ievt))
#npuppijetsnotiming2.Scale(1./float(ievt))
#nchsjetsnotiming.Scale(1./float(ievt))
#nchsjetsnotiming2.Scale(1./float(ievt))
#nAK8jetsnotiming.Scale(1./float(ievt))
#nAK8jetsnotiming2.Scale(1./float(ievt))


####FIX PLOTS for jet counting
c1 = ROOT.TCanvas()
npuppijetsnotiming.Draw("HIST")
npuppijetstiming.Draw("HISTSAME")
c1.SaveAs("puppinopu_tight.pdf")
c1.SaveAs("puppinopu_tight.root")

c2 = ROOT.TCanvas()
nchsjetsnotiming.Draw("HIST")
nchsjetstiming.Draw("HISTSAME")
c2.SaveAs("chsnopu_tight.pdf")
c2.SaveAs("chsnopu_tight.root")

c3 = ROOT.TCanvas()
nAK8jetsnotiming.Draw("HIST")
nAK8jetstiming.Draw("HISTSAME")
c3.SaveAs("AK8nopu_tight.pdf")
c3.SaveAs("AK8nopu_tight.root")

c4 = ROOT.TCanvas()
npuppijetsnotiming2.Draw("HIST")
npuppijetstiming2.Draw("HISTSAME")
c4.SaveAs("puppinopu_loose.pdf")
c4.SaveAs("puppinopu_loose.root")

c5 = ROOT.TCanvas()
nchsjetsnotiming2.Draw("HIST")
nchsjetstiming2.Draw("HISTSAME")
c5.SaveAs("chsnopu_loose.pdf")
c5.SaveAs("chsnopu_loose.root")

c6 = ROOT.TCanvas()
nAK8jetsnotiming2.Draw("HIST")
nAK8jetstiming2.Draw("HISTSAME")
c6.SaveAs("AK8nopu_loose.pdf")
c6.SaveAs("AK8nopu_loose.root")


input("Press Enter to continue...")
