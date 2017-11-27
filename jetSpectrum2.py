#! /usr/bin/env python

import ROOT
from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
import os
import math
import gc
import array

gc.disable()

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


nfiles = -1

files={}
indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/PU200_timing/171110_104722/0000/'
files['PU200_timing'] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/NoPU_timing_/171109_155602/0000/'
files['NoPU_timing_'] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1sigma_PU200_timing/171116_091208/0000/'
files['1sigma_PU200_timing'] = [indir+filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1sigma_NoPU_timing/171116_091101/0000/'
files['1sigma_NoPU_timing'] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_PU200_timing/171116_095856/0000/'
files['1p5sigma_PU200_timing'] = [indir +filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_NoPU_timing/171116_095748/0000/'
files['1p5sigma_NoPU_timing'] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]


indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/PU200_notiming/171110_104953/0000/'
files['PU200_notiming'] = [indir +filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/NoPU_notiming_/171109_155445/0000/'
files['NoPU_notiming_'] = [indir +filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1sigma_PU200_notiming/171116_091300/0000/'
files['1sigma_PU200_notiming'] = [indir +filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1sigma_NoPU_notiming/171116_090950/0000/'
files['1sigma_NoPU_notiming'] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_PU200_notiming/171116_095955/0000/'
files['1p5sigma_PU200_notiming'] = [indir +filename for filename in os.listdir(indir)[0:nfiles]]

iindir = '/eos/cms/store/group/phys_higgs/cmshmm/amarini/RelValZMM_14/RelValZMM_14/1p5sigma_NoPU_notiming/171116_095645/0000/'
files['1p5sigma_NoPU_notiming'] = [indir +filename for filename in os.listdir(indir)[0:nfiles]]

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

histos={}
roc={}
def loop(inputfiles,label):
    global histos
    global roc
    print "-> Doing",label, inputfiles[0]
    events = Events (inputFiles=inputfiles)

    handlevtx  = Handle ("std::vector<reco::Vertex>")
    labelvtx = ("offlineSlimmedPrimaryVertices")

    handlepfcand = Handle("vector<pat::PackedCandidate>")
    labelpfcand = ("packedPFCandidates")

    handlejetchs = Handle("std::vector<pat::Jet>")
    labeljetchs = "slimmedJets"

    handlejetpuppi = Handle("std::vector<pat::Jet>")
    labeljetpuppi = "slimmedJetsPuppi"

    handlegenjet = Handle("std::vector<reco::GenJet>")
    labelgenjet = "slimmedGenJets"

    handleinfo = Handle("GenEventInfoProduct")
    labelinfo = "generator"


    # Create histograms, etc.
    #ROOT.gROOT.SetBatch()        # don't pop up canvases
    #ROOT.gROOT.SetStyle('Plain') # white background

    histos["matched_eta1_"+label]=ROOT.TH1D("matched_eta1_"+label,label,50,0.0,1.44)
    histos["unmatched_eta1_"+label]=ROOT.TH1D("unmatched_eta1_"+label,label,50,0.0,1.44)
	
    histos["matched_eta2_"+label]=ROOT.TH1D("matched_eta2_"+label,label,50,1.44,2.5)
    histos["unmatched_eta2_"+label]=ROOT.TH1D("unmatched_eta2_"+label,label,50,1.44,2.5)
	
    histos["matched_eta3_"+label]=ROOT.TH1D("matched_eta3_"+label,label,50,2.5,5.)
    histos["unmatched_eta3_"+label]=ROOT.TH1D("unmatched_eta3_"+label,label,50,2.5,5.)

	
    maxevt = 100

    minjetpt = 30.

    # loop over events
    ievt=0
    roc['eff_'+label]=0.
    roc['fake_'+label]=0.
    for event in events:

        event.getByLabel (labelvtx, handlevtx)
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
	    for genpart in genjet.getJetConstituentsQuick():
                #print("genpart")
                #genpart = genjet.getGenConstituent(ipart)
                if abs(genpart.pdgId())==13 and genpart.fromHardProcessFinalState():
                    #printf("match")
                    mugenjets.append(genjet)
                    break
            if len(mugenjets)==2:
                break
              
        for jet in puppijets:
	    #print(jet.jecFactor("Uncorrected")/jet.jecFactor("L2Relative"))
      
            if jet.pt()>minjetpt:
                mugenjet = matchgenjetloose(jet,mugenjets)
                if mugenjet:
                    continue
          
                genjet = matchgenjet(jet,genjets)
	        ispu = ispujet(jet,genjets)

		redo=0
	        if genjet:
		    try:
		        if abs(jet.eta())<1.44: 
			    histos['matched_eta1_'+label].Fill(abs(jet.eta()))
			    roc['eff_'+label]+=1.
		        if abs(jet.eta())>1.44 and abs(jet.eta())<2.5: histos['matched_eta2_'+label].Fill(abs(jet.eta()))
		        if abs(jet.eta())>2.5 and abs(jet.eta())<5: histos['matched_eta3_'+label].Fill(abs(jet.eta()))
		    except:
			print("matched broke :(")
			redo=1
			break

	        if ispu:
		    try:
		        if abs(jet.eta())<1.44: 
			    histos['unmatched_eta1_'+label].Fill(abs(jet.eta()))
			    roc['fake_'+label]+=1.
		        if abs(jet.eta())>1.44 and abs(jet.eta())<2.5: histos['unmatched_eta2_'+label].Fill(abs(jet.eta()))
		        if abs(jet.eta())>2.5 and abs(jet.eta())<5: histos['unmatched_eta3_'+label].Fill(abs(jet.eta()))
		    except:
			print("unmatched broke :(")
                        redo=1
                        break

        ievt += 1
        if ievt==maxevt:
            break;
 

    if redo==1:
	redo = 0	
	loop(inputfiles,label)

    print(ievt)

    try:
        histos['matched_eta1_'+label].Scale(1./float(ievt))
        histos['matched_eta2_'+label].Scale(1./float(ievt))
        histos['matched_eta3_'+label].Scale(1./float(ievt))
	
        histos['unmatched_eta1_'+label].Scale(1./float(ievt))
        histos['unmatched_eta2_'+label].Scale(1./float(ievt))
        histos['unmatched_eta3_'+label].Scale(1./float(ievt))
    except:
    	print("unmatched broke :(")
        redo=1

    if redo==1:
        redo = 0
        loop(inputfiles,label)

    print(roc)


for index in files:
    loop(files[index],index)


ratio={}
#loop for 1.5 and 1 sigma, efficiency and fake ratios
for eta in ["eta1","eta2","eta3"]:
    for sig in ['1sigma','1p5sigma']:
	for timing in ['timing','notiming']:
	    ratio['eff_'+sig+'_'+eta+'_'+timing] = histos['matched_'+eta+'_'+sig+'_'+'PU200'+'_'+timing].Clone()
	    ratio['eff_'+sig+'_'+eta+'_'+timing].Divide(histos['matched_'+eta+'_'+sig+'_'+'NoPU'+'_'+timing])
	    ratio['fake_'+sig+'_'+eta+'_'+timing] = histos['unmatched_'+eta+'_'+sig+'_'+'PU200'+'_'+timing].Clone()
	    ratio['fake_'+sig+'_'+eta+'_'+timing].Divide(histos['matched_'+eta+'_'+sig+'_'+'NoPU'+'_'+timing])

		
#loop for 2 sigma, efiiciency and fake ratios
for eta in ["eta1","eta2","eta3"]:
    for timing in ['timing','notiming']:
	ratio['eff_2sigma_'+eta+'_'+timing] = histos['matched_'+eta+'_'+'PU200'+'_'+timing].Clone()
	ratio['eff_2sigma_'+eta+'_'+timing].Divide(histos['matched_'+eta+'_'+'NoPU'+'_'+timing+'_'])
		
	ratio['fake_2sigma_'+eta+'_'+timing] = histos['unmatched_'+eta+'_'+'PU200'+'_'+timing].Clone()
	ratio['fake_2sigma_'+eta+'_'+timing].Divide(histos['matched_'+eta+'_'+'NoPU'+'_'+timing+'_'])
 

efflist=[]
fakelist=[]
for timing in ['timing','notiming']:
    for sig in ['1sigma','1p5sigma']:
        efflist.append(roc['eff_'+sig+'_PU200_'+timing]/roc['eff_'+sig+'_NoPU_'+timing])
	fakelist.append(roc['fake_'+sig+'_PU200_'+timing]/roc['eff_'+sig+'_NoPU_'+timing])

    efflist.append(roc['eff_PU200_'+timing]/roc['eff_NoPU_'+timing+'_'])
    fakelist.append(roc['fake_PU200_'+timing]/roc['eff_NoPU_'+timing+'_'])

efflist_t=[]
fakelist_t=[]
for i in range(0,3): 
    efflist_t.append(efflist[i])
    fakelist_t.append(fakelist[i])

efflist_nt=[]
fakelist_nt=[]
for i in range(3,6): 
    efflist_nt.append(efflist[i])
    fakelist_nt.append(fakelist[i])

et=array.array('d',efflist_t)
ft=array.array('d',fakelist_t)
ent=array.array('d',efflist_nt)
fnt=array.array('d',fakelist_nt)
##NEED TO MAKE CANVASES, then done (after ROC)
canvas={}
count=1
for eta in ["eta1","eta2","eta3"]:
    for sig in ['1sigma','1p5sigma','2sigma']:
	for ef in ['eff','fake']:
	    canvas['c'+str(count)] = ROOT.TCanvas()
	    ratio[ef+'_'+sig+'_'+eta+'_notiming'].SetLineColor(ROOT.kRed)
	    ratio[ef+'_'+sig+'_'+eta+'_notiming'].Draw("HIST")
	    ratio[ef+'_'+sig+'_'+eta+'_timing'].Draw("HISTSAME")
	    
            leg = ROOT.TLegend(0.4,.8,.6,.9)
            leg.AddEntry(ratio[ef+'_'+sig+'_'+eta+'_timing'],"timing","L")
            leg.AddEntry(ratio[ef+'_'+sig+'_'+eta+'_notiming'],"no-timing","L")
            leg.Draw()
	
	    ratio[ef+'_'+sig+'_'+eta+'_notiming'].GetXaxis().SetTitle("|#eta|")
	
	    canvas['c'+str(count)].SaveAs(ef+'_'+sig+'_'+eta+'.pdf')
	    canvas['c'+str(count)].SaveAs(ef+'_'+sig+'_'+eta+'.root')
	
	    count += 1

canvas['c'+str(count)] = ROOT.TCanvas()
roc_t = ROOT.TGraph(3,et,ft)
roc_nt = ROOT.TGraph(3,ent,fnt)
roc_nt.SetLineColor(2)
roc_t.SetLineColor(4)
roc_t.SetMarkerColor(1)
roc_nt.SetMarkerColor(1)
roc_t.SetMarkerSize(2)
roc_nt.SetMarkerSize(2)
roc_t.SetMarkerStyle(3)
roc_nt.SetMarkerStyle(3)

mg=ROOT.TMultiGraph()
mg.Add(roc_nt)
mg.Add(roc_t)
mg.Draw("ACP")

leg = ROOT.TLegend(0.8,.8,.9,.9)
leg.AddEntry(roc_t,"Timing","L")
leg.AddEntry(roc_nt,"No Timing","L")
leg.Draw()

mg.GetXaxis().SetTitle("Efficiency")
mg.GetYaxis().SetTitle("Rejection")

canvas['c'+str(count)].SaveAs('roc.pdf')
canvas['c'+str(count)].SaveAs('roc.root')

input("Press Enter to continue...")
