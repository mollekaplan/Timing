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

######No Pileup Files######

nfiles = -1

##Timing##
#indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme_reminiaod_vtxassign/dymm200timing/MINIAODSIM/'
indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Aug26jme/dymmnoputiming/MINIAODSIM/'
inputfiles = [indir + filename for filename in os.listdir(indir)[0:nfiles]]
#for filename in os.listdir(indir)[0:nfiles]:
#    if filename != indir+'step3_dymm200timing_MINIAODSIM_119.root':
#	inputfiles = indir + filename

##No Timing##
#indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Mar31jme_reminiaod_vtxassign/dymm200notiming/MINIAODSIM/'
indir2 = '/eos/cms/store/group/upgrade/timing/pfintegration/Aug26jme/dymmnopunotiming/MINIAODSIM/'
inputfiles2 = [indir2 + filename for filename in os.listdir(indir2)[0:nfiles]]


######Pileup Files######

##Timing##
indir3 = '/eos/cms/store/group/upgrade/timing/pfintegration/Sep5jme/dymm200timing/MINIAODSIM/'
#inputfiles = [indir + filename for filename in os.listdir(indir)[0:nfiles]]
inputfiles3=[]
for filename in os.listdir(indir3)[0:nfiles]:
    if filename != 'step3_dymm200timing_MINIAODSIM_119.root' and filename != 'step3_dymm200timing_MINIAODSIM_93.root':
	inputfiles3.append(indir3 + filename)
	
##No Timing##
indir4 = '/eos/cms/store/group/upgrade/timing/pfintegration/Sep5jme/dymm200notiming/MINIAODSIM/'
inputfiles4 = [indir4 + filename for filename in os.listdir(indir4)[0:nfiles]]


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
events3 = Events (inputFiles=inputfiles3)
events4 = Events (inputFiles=inputfiles4) 


# create handle outside of loop
handlepfcand = Handle("vector<pat::PackedCandidate>")
labelpfcand = ("packedPFCandidates")

handlejetpuppi = Handle("std::vector<pat::Jet>")
labeljetpuppi = "slimmedJetsPuppi"

handlechsjet = Handle("std::vector<pat::Jet>")
labelchsjet = "slimmedJets"

handlejetAK8 = Handle("std::vector<pat::Jet>")
labeljetAK8 = "slimmedJetsAK8"

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
ROOT.gStyle.SetOptStat(0);

######Hists######
##No Pileup Graphs##
npuppijetstimingnp = ROOT.TH1D("npuppijetstimingnp","Number of Puppi Jets",20,0,20)
npuppijetsnotimingnp = ROOT.TH1D("npuppijetsnotimingnp","Number of Puppi Jets",20,0,20)
npuppijetsnotimingnp.SetLineColor(ROOT.kRed)

nAK8jetstimingnp = ROOT.TH1D("nAK8jetstimingnp","Number of AK8 Jets",10,0,10)
nAK8jetsnotimingnp = ROOT.TH1D("nAK8jetsnotimingnp","Number of AK8 Jets",10,0,10)
nAK8jetsnotimingnp.SetLineColor(ROOT.kRed)

npuppijetstimingnp2 = ROOT.TH1D("npuppijetstimingnp2","Number of Puppi Jets",20,0,20)
npuppijetsnotimingnp2 = ROOT.TH1D("npuppijetsnotimingnp2","Number of Puppi Jets",20,0,20)
npuppijetsnotimingnp2.SetLineColor(ROOT.kRed)

nAK8jetstimingnp2 = ROOT.TH1D("nAK8jetstimingnp2","Number of AK8 Jets",10,0,10)
nAK8jetsnotimingnp2 = ROOT.TH1D("nAK8jetsnotimingnp2","Number of AK8 Jets",10,0,10)
nAK8jetsnotimingnp2.SetLineColor(ROOT.kRed)


##Pileup Graphs##
npuppijetstiming = ROOT.TH1D("npuppijetstiming","Number of Puppi Jets",20,0,20)
npuppijetsnotiming = ROOT.TH1D("npuppijetsnotiming","Number of Puppi Jets",20,0,20)
npuppijetsnotiming.SetLineColor(ROOT.kRed)

nAK8jetstiming = ROOT.TH1D("nAK8jetstiming","Number of AK8 Jets",10,0,10)
nAK8jetsnotiming = ROOT.TH1D("nAK8jetsnotiming","Number of AK8 Jets",10,0,10)
nAK8jetsnotiming.SetLineColor(ROOT.kRed)

npuppijetstiming2 = ROOT.TH1D("npuppijetstiming2","Number of Puppi Jets",20,0,20)
npuppijetsnotiming2 = ROOT.TH1D("npuppijetsnotiming2","Number of Puppi Jets",20,0,20)
npuppijetsnotiming2.SetLineColor(ROOT.kRed)

nAK8jetstiming2 = ROOT.TH1D("nAK8jetstiming2","Number of AK8 Jets",10,0,10)
nAK8jetsnotiming2 = ROOT.TH1D("nAK8jetsnotiming2","Number of AK8 Jets",10,0,10)
nAK8jetsnotiming2.SetLineColor(ROOT.kRed)


######Running Events######
maxevt = 2500

minjetpt = 30.

################ TIMING, NO PILEUP #################
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
        mugenjet = matchgenjet(jet,mugenjets)
	if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            npuppi2 += 1
            if abs(jet.eta()) < 2.5:
                npuppi += 1
    npuppijetstimingnp.Fill(npuppi)
    npuppijetstimingnp2.Fill(npuppi2)

	
    nAK8=0
    nAK82=0
    for jet in jetsAK8:
        mugenjet = matchgenjet(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            nAK82 += 1
            if abs(jet.eta()) < 2.5:
                nAK8 += 1
    nAK8jetstimingnp.Fill(nAK8)
    nAK8jetstimingnp2.Fill(nAK82)


    ievt += 1
    if ievt==maxevt:
        break;

#npuppijetstiming.Scale(1./float(ievt))
#npuppijetstiming2.Scale(1./float(ievt))
#nchsjetstiming.Scale(1./float(ievt))
#nchsjetstiming2.Scale(1./float(ievt))
#nAK8jetstiming.Scale(1./float(ievt))
#nAK8jetstiming2.Scale(1./float(ievt))




################## NO TIMING, NO PILEUP #####################
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
        mugenjet = matchgenjet(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            npuppi2 += 1
            if abs(jet.eta()) < 2.5:
                npuppi += 1
    npuppijetsnotimingnp.Fill(npuppi)
    npuppijetsnotimingnp2.Fill(npuppi2)
	
    
    nAK8=0
    nAK82=0
    for jet in jetsAK8:
        mugenjet = matchgenjet(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            nAK82 += 1
            if abs(jet.eta()) < 2.5:
                nAK8 += 1
    nAK8jetsnotimingnp.Fill(nAK8)
    nAK8jetsnotimingnp2.Fill(nAK82)


    ievt += 1
    if ievt==maxevt:
        break;


#npuppijetsnotiming.Scale(1./float(ievt))
#npuppijetsnotiming2.Scale(1./float(ievt))
#nchsjetsnotiming.Scale(1./float(ievt))
#nchsjetsnotiming2.Scale(1./float(ievt))
#nAK8jetsnotiming.Scale(1./float(ievt))
#nAK8jetsnotiming2.Scale(1./float(ievt))



################ TIMING, PILEUP #################
ievt=0
for event in events3:

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
        mugenjet = matchgenjet(jet,mugenjets)
	if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            npuppi2 += 1
            if abs(jet.eta()) < 2.5:
                npuppi += 1
    npuppijetstiming.Fill(npuppi)
    npuppijetstiming2.Fill(npuppi2)

	
    nAK8=0
    nAK82=0
    for jet in jetsAK8:
        mugenjet = matchgenjet(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            nAK82 += 1
            if abs(jet.eta()) < 2.5:
                nAK8 += 1
    nAK8jetstiming.Fill(nAK8)
    nAK8jetstiming2.Fill(nAK82)


    ievt += 1
    if ievt==maxevt:
        break;
	
	
	
################## NO TIMING, PILEUP #####################
ievt = 0
for event in events4:

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
        mugenjet = matchgenjet(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            npuppi2 += 1
            if abs(jet.eta()) < 2.5:
                npuppi += 1
    npuppijetsnotiming.Fill(npuppi)
    npuppijetsnotiming2.Fill(npuppi2)
	
    
    nAK8=0
    nAK82=0
    for jet in jetsAK8:
        mugenjet = matchgenjet(jet,mugenjets)
        if abs(jet.eta()) < 4.7 and not mugenjet and jet.pt()>minjetpt:
            nAK82 += 1
            if abs(jet.eta()) < 2.5:
                nAK8 += 1
    nAK8jetsnotiming.Fill(nAK8)
    nAK8jetsnotiming2.Fill(nAK82)


    ievt += 1
    if ievt==maxevt:
        break;
	
	
	
	
	
	

######Make canvases and save the plots######

####No Pileup####
##c1##
c1 = ROOT.TCanvas()
npuppijetstimingnp.Draw("HIST")
npuppijetsnotimingnp.Draw("HISTSAME")

leg1 = ROOT.TLegend(0.7,.5,.89,.85)
leg1.AddEntry(npuppijetstimingnp,"timing","L")
leg1.AddEntry(npuppijetsnotimingnp,"no-timing","L")
leg1.Draw()

c1.SaveAs("puppinopu_tight.pdf")
c1.SaveAs("puppinopu_tight.root")

##c2##
c2 = ROOT.TCanvas()
nAK8jetstimingnp.Draw("HIST")
nAK8jetsnotimingnp.Draw("HISTSAME")

leg2 = ROOT.TLegend(0.7,.5,.89,.85)
leg2.AddEntry(nAK8jetstimingnp,"timing","L")
leg2.AddEntry(nAK8jetsnotimingnp,"no-timing","L")
leg2.Draw()

c2.SaveAs("AK8nopu_tight.pdf")
c2.SaveAs("AK8nopu_tight.root")

##c3##
c3 = ROOT.TCanvas()
npuppijetstimingnp2.Draw("HIST")
npuppijetsnotimingnp2.Draw("HISTSAME")

leg3 = ROOT.TLegend(0.7,.5,.89,.85)
leg3.AddEntry(npuppijetstimingnp2,"timing","L")
leg3.AddEntry(npuppijetsnotimingnp2,"no-timing","L")
leg3.Draw()

c3.SaveAs("puppinopu_loose.pdf")
c3.SaveAs("puppinopu_loose.root")

##c4##
c4 = ROOT.TCanvas()
nAK8jetstimingnp2.Draw("HIST")
nAK8jetsnotimingnp2.Draw("HISTSAME")

leg4 = ROOT.TLegend(0.7,.5,.89,.85)
leg4.AddEntry(nAK8jetstimingnp2,"timing","L")
leg4.AddEntry(nAK8jetsnotimingnp2,"no-timing","L")
leg4.Draw()

c4.SaveAs("AK8nopu_loose.pdf")
c4.SaveAs("AK8nopu_loose.root")

##c5##
c5 = c1.DrawClone()
c5.SetLogy()
c5.Draw()

leg5 = ROOT.TLegend(0.7,.5,.89,.85)
leg5.AddEntry(npuppijetstimingnp,"timing","L")
leg5.AddEntry(npuppijetsnotimingnp,"no-timing","L")
leg5.Draw()

c5.SaveAs("puppinopu_tight_log.pdf")
c5.SaveAs("puppinopu_tight_log.root")

##c6##
c6 = c3.DrawClone()
c6.SetLogy()
c6.Draw()

leg6 = ROOT.TLegend(0.7,.5,.89,.85)
leg6.AddEntry(npuppijetstimingnp2,"timing","L")
leg6.AddEntry(npuppijetsnotimingnp2,"no-timing","L")
leg6.Draw()

c6.SaveAs("puppinopu_loose_log.pdf")
c6.SaveAs("puppinopu_loose_log.root")



####Pileup####
##c7##
c7 = ROOT.TCanvas()
npuppijetsnotiming.Draw("HIST")
npuppijetstiming.Draw("HISTSAME")

leg7 = ROOT.TLegend(0.7,.5,.89,.85)
leg7.AddEntry(npuppijetstiming,"timing","L")
leg7.AddEntry(npuppijetsnotiming,"no-timing","L")
leg7.Draw()

c7.SaveAs("puppi_tight.pdf")
c7.SaveAs("puppi_tight.root")

##c8##
c8 = ROOT.TCanvas()
nAK8jetstiming.Draw("HIST")
nAK8jetsnotiming.Draw("HISTSAME")

leg8 = ROOT.TLegend(0.7,.5,.89,.85)
leg8.AddEntry(nAK8jetstiming,"timing","L")
leg8.AddEntry(nAK8jetsnotiming,"no-timing","L")
leg8.Draw()

c8.SaveAs("AK8_tight.pdf")
c8.SaveAs("AK8_tight.root")

##c9##
c9 = ROOT.TCanvas()
npuppijetsnotiming2.Draw("HIST")
npuppijetstiming2.Draw("HISTSAME")

leg9 = ROOT.TLegend(0.7,.5,.89,.85)
leg9.AddEntry(npuppijetstiming2,"timing","L")
leg9.AddEntry(npuppijetsnotiming2,"no-timing","L")
leg9.Draw()

c9.SaveAs("puppi_loose.pdf")
c9.SaveAs("puppi_loose.root")

##c10##
c10 = ROOT.TCanvas()
nAK8jetstiming2.Draw("HIST")
nAK8jetsnotiming2.Draw("HISTSAME")

leg10 = ROOT.TLegend(0.7,.5,.89,.85)
leg10.AddEntry(nAK8jetstiming2,"timing","L")
leg10.AddEntry(nAK8jetsnotiming2,"no-timing","L")
leg10.Draw()

c10.SaveAs("AK8_loose.pdf")
c10.SaveAs("AK8_loose.root")

##c11##
c11 = c7.DrawClone()
c11.SetLogy()
c11.Draw()

leg11 = ROOT.TLegend(0.7,.5,.89,.85)
leg11.AddEntry(npuppijetstiming,"timing","L")
leg11.AddEntry(npuppijetsnotiming,"no-timing","L")
leg11.Draw()

c11.SaveAs("puppi_tight_log.pdf")
c11.SaveAs("puppi_tight_log.root")

##c12##
c12 = c9.DrawClone()
c12.SetLogy()
c12.Draw()

leg12 = ROOT.TLegend(0.7,.5,.89,.85)
leg12.AddEntry(npuppijetstiming2,"timing","L")
leg12.AddEntry(npuppijetsnotiming2,"no-timing","L")
leg12.Draw()

c12.SaveAs("puppi_loose_log.pdf")
c12.SaveAs("puppi_loose_log.root")



####No pileup vs. pileup, timing####
##c13##
npuppijetstimingnp.SetLineColor(ROOT.kRed)

c13 = ROOT.TCanvas()
npuppijetstiming.Draw("HIST")
npuppijetstimingnp.Draw("HISTSAME")

leg13 = ROOT.TLegend(0.7,.5,.89,.85)
leg13.AddEntry(npuppijetstiming,"with pileup","L")
leg13.AddEntry(npuppijetstimingnp,"without pileup","L")
leg13.Draw()

c13.SaveAs("puppivs_tight.pdf")
c13.SaveAs("puppivs_tight.root")

##c14##
nAK8jetstimingnp.SetLineColor(ROOT.kRed)

c14 = ROOT.TCanvas()
nAK8jetstiming.Draw("HIST")
nAK8jetstimingnp.Draw("HISTSAME")

leg14 = ROOT.TLegend(0.7,.5,.89,.85)
leg14.AddEntry(nAK8jetstiming,"with pileup","L")
leg14.AddEntry(nAK8jetstimingnp,"without pileup","L")
leg14.Draw()

c14.SaveAs("AK8vs_tight.pdf")
c14.SaveAs("AK8vs_tight.root")

##c15##
npuppijetstimingnp2.SetLineColor(ROOT.kRed)

c15 = ROOT.TCanvas()
npuppijetstiming2.Draw("HIST")
npuppijetstimingnp2.Draw("HISTSAME")

leg15 = ROOT.TLegend(0.7,.5,.89,.85)
leg15.AddEntry(npuppijetstiming2,"with pileup","L")
leg15.AddEntry(npuppijetstimingnp2,"without pileup","L")
leg15.Draw()

c15.SaveAs("puppivs_loose.pdf")
c15.SaveAs("puppivs_loose.root")

##c16##
nAK8jetstimingnp2.SetLineColor(ROOT.kRed)

c16 = ROOT.TCanvas()
nAK8jetstiming2.Draw("HIST")
nAK8jetstimingnp2.Draw("HISTSAME")

leg16 = ROOT.TLegend(0.7,.5,.89,.85)
leg16.AddEntry(nAK8jetstiming2,"with pileup","L")
leg16.AddEntry(nAK8jetstimingnp2,"without pileup","L")
leg16.Draw()

c16.SaveAs("AK8vs_loose.pdf")
c16.SaveAs("AK8vs_loose.root")

##c17##
c17 = c13.DrawClone()
c17.SetLogy()
c17.Draw()

leg17 = ROOT.TLegend(0.7,.5,.89,.85)
leg17.AddEntry(npuppijetstiming,"with pileup","L")
leg17.AddEntry(npuppijetstimingnp,"without pileup","L")
leg17.Draw()

c17.SaveAs("puppivs_tight_log.pdf")
c17.SaveAs("puppivs_tight_log.root")

##c18##
c18 = c15.DrawClone()
c18.SetLogy()
c18.Draw()

leg18 = ROOT.TLegend(0.7,.5,.89,.85)
leg18.AddEntry(npuppijetstiming2,"with pileup","L")
leg18.AddEntry(npuppijetstimingnp2,"without pileup","L")
leg18.Draw()

c18.SaveAs("puppivs_loose_log.pdf")
c18.SaveAs("puppivs_loose_log.root")



####No pileup vs. pileup, no timing####
##c19##
npuppijetsnotiming.SetLineColor(ROOT.kBlue)
npuppijetsnotimingnp.SetLineColor(ROOT.kRed)

c19 = ROOT.TCanvas()
npuppijetsnotiming.Draw("HIST")
npuppijetsnotimingnp.Draw("HISTSAME")

leg19 = ROOT.TLegend(0.7,.5,.89,.85)
leg19.AddEntry(npuppijetsnotiming,"with pileup","L")
leg19.AddEntry(npuppijetsnotimingnp,"without pileup","L")
leg19.Draw()

c19.SaveAs("puppivsnt_tight.pdf")
c19.SaveAs("puppivsnt_tight.root")

##c20##
nAK8jetsnotiming.SetLineColor(ROOT.kBlue)
nAK8jetsnotimingnp.SetLineColor(ROOT.kRed)

c20 = ROOT.TCanvas()
nAK8jetsnotiming.Draw("HIST")
nAK8jetsnotimingnp.Draw("HISTSAME")

leg20 = ROOT.TLegend(0.7,.5,.89,.85)
leg20.AddEntry(nAK8jetsnotiming,"with pileup","L")
leg20.AddEntry(nAK8jetsnotimingnp,"without pileup","L")
leg20.Draw()

c20.SaveAs("AK8vsnt_tight.pdf")
c20.SaveAs("AK8vsnt_tight.root")

##c21##
npuppijetsnotiming2.SetLineColor(ROOT.kBlue)
npuppijetsnotimingnp2.SetLineColor(ROOT.kRed)

c21 = ROOT.TCanvas()
npuppijetsnotiming2.Draw("HIST")
npuppijetsnotimingnp2.Draw("HISTSAME")

leg21 = ROOT.TLegend(0.7,.5,.89,.85)
leg21.AddEntry(npuppijetsnotiming2,"with pileup","L")
leg21.AddEntry(npuppijetsnotimingnp2,"without pileup","L")
leg21.Draw()

c21.SaveAs("puppivsnt_loose.pdf")
c21.SaveAs("puppivsnt_loose.root")

##c22##
nAK8jetsnotiming2.SetLineColor(ROOT.kBlue)
nAK8jetsnotimingnp2.SetLineColor(ROOT.kRed)

c22 = ROOT.TCanvas()
nAK8jetsnotiming2.Draw("HIST")
nAK8jetsnotimingnp2.Draw("HISTSAME")

leg22 = ROOT.TLegend(0.7,.5,.89,.85)
leg22.AddEntry(nAK8jetsnotiming2,"with pileup","L")
leg22.AddEntry(nAK8jetsnotimingnp2,"without pileup","L")
leg22.Draw()

c22.SaveAs("AK8vsnt_loose.pdf")
c22.SaveAs("AK8vsnt_loose.root")

##c23##
c23 = c19.DrawClone()
c23.SetLogy()
c23.Draw()

leg23 = ROOT.TLegend(0.7,.5,.89,.85)
leg23.AddEntry(npuppijetsnotiming,"with pileup","L")
leg23.AddEntry(npuppijetsnotimingnp,"without pileup","L")
leg23.Draw()

c23.SaveAs("puppivsnt_tight_log.pdf")
c23.SaveAs("puppivsnt_tight_log.root")

##c24##
c24 = c21.DrawClone()
c24.SetLogy()
c24.Draw()

leg24 = ROOT.TLegend(0.7,.5,.89,.85)
leg24.AddEntry(npuppijetsnotiming2,"with pileup","L")
leg24.AddEntry(npuppijetsnotimingnp2,"without pileup","L")
leg24.Draw()

c24.SaveAs("puppivsnt_loose_log.pdf")
c24.SaveAs("puppivsnt_loose_log.root")


input("Press Enter to continue...")
