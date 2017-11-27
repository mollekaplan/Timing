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
    roc['eff_eta1_'+label]=0.
    roc['fake_eta1_'+label]=0.
    roc['eff_eta2_'+label]=0.
    roc['fake_eta2_'+label]=0.
    roc['eff_eta3_'+label]=0.
    roc['fake_eta3_'+label]=0.
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
		    if abs(jet.eta())<1.44: 
			histos['matched_eta1_'+label].Fill(abs(jet.eta()))
			roc['eff_eta1_'+label]+=1.
		    if abs(jet.eta())>1.44 and abs(jet.eta())<2.5: 
			histos['matched_eta2_'+label].Fill(abs(jet.eta()))
			roc['eff_eta2_'+label]+=1.
		    if abs(jet.eta())>2.5 and abs(jet.eta())<5: 
			histos['matched_eta3_'+label].Fill(abs(jet.eta()))
			roc['eff_eta3_'+label]+=1.
		    if abs(jet.eta())<5:
			roc['eff_tot_'+label]+=1.
			
	        if ispu:
		    if abs(jet.eta())<1.44: 
			histos['unmatched_eta1_'+label].Fill(abs(jet.eta()))
			roc['fake_eta1_'+label]+=1.
		    if abs(jet.eta())>1.44 and abs(jet.eta())<2.5: 
			histos['unmatched_eta2_'+label].Fill(abs(jet.eta()))
			roc['fake_eta2_'+label]+=1.
		    if abs(jet.eta())>2.5 and abs(jet.eta())<5: 
			histos['unmatched_eta3_'+label].Fill(abs(jet.eta()))
			roc['fake_eta3_'+label]+=1.
		    if abs(jet.eta())<5:
			roc['fake_tot_'+label]+=1.

        ievt += 1
        if ievt==maxevt:
            break;

    try:
    histos['matched_eta1_'+label].Scale(1./float(ievt))
    histos['matched_eta2_'+label].Scale(1./float(ievt))
    histos['matched_eta3_'+label].Scale(1./float(ievt))
	
    histos['unmatched_eta1_'+label].Scale(1./float(ievt))
    histos['unmatched_eta2_'+label].Scale(1./float(ievt))
    histos['unmatched_eta3_'+label].Scale(1./float(ievt))

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
 

for eta in ['eta1','eta2','eta3','tot']:
    for sig in ['1sigma','1p5sigma']:
        efflist_t.append(roc['eff_'+eta+'_'+sig+'_PU200_timing']/roc['eff_'+eta+'_'+sig+'_NoPU_timing'])
        fakelist_t.append(roc['fake_'+eta+'_'+sig+'_PU200_timing']/roc['eff_'+eta+'_'+sig+'_NoPU_timing'])
        efflist_nt.append(roc['eff_'+eta+'_'+sig+'_PU200_notiming']/roc['eff_'+eta+'_'+sig+'_NoPU_notiming'])
        fakelist_nt.append(roc['fake_'+eta+'_'+sig+'_PU200_notiming']/roc['eff_'+eta+'_'+sig+'_NoPU_notiming'])

    efflist_t.append(roc['eff_'+eta+'_PU200_timing']/roc['eff_'+eta+'_NoPU_timing_'])
    fakelist_t.append(roc['fake_'+eta+'_PU200_timing']/roc['eff_'+eta+'_NoPU_timing_'])
    efflist_nt.append(roc['eff_'+eta+'_PU200_notiming']/roc['eff_'+eta+'_NoPU_notiming_'])
    fakelist_nt.append(roc['fake_'+eta+'_PU200_notiming']/roc['eff_'+eta+'_NoPU_notiming_'])

et1=array.array('d',efflist_t[0:3])
ft1=array.array('d',fakelist_t[0:3])
ent1=array.array('d',efflist_nt[0:3])
fnt1=array.array('d',fakelist_nt[0:3])

et2=array.array('d',efflist_t[3:6])
ft2=array.array('d',fakelist_t[3:6])
ent2=array.array('d',efflist_nt[3:6])
fnt2=array.array('d',fakelist_nt[3:6])

et3=array.array('d',efflist_t[6:9])
ft3=array.array('d',fakelist_t[6:9])
ent3=array.array('d',efflist_nt[6:9])
fnt3=array.array('d',fakelist_nt[6:9])

ettot=array.array('d',efflist_t[9:12])
fttot=array.array('d',fakelist_t[9:12])
enttot=array.array('d',efflist_nt[9:12])
fnttot=array.array('d',fakelist_nt[9:12])


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

##Will clean up and make loop later
##eta1##
canvas['c'+str(count)] = ROOT.TCanvas()
roc_t1 = ROOT.TGraph(3,et1,ft1)
roc_nt1 = ROOT.TGraph(3,ent1,fnt1)
roc_nt1.SetLineColor(2)
roc_t1.SetLineColor(4)
roc_t1.SetMarkerColor(1)
roc_nt1.SetMarkerColor(1)
roc_t1.SetMarkerSize(2)
roc_nt1.SetMarkerSize(2)
roc_t1.SetMarkerStyle(3)
roc_nt1.SetMarkerStyle(3)

mg1=ROOT.TMultiGraph()
mg1.Add(roc_nt1)
mg1.Add(roc_t1)
mg1.Draw("ACP")

leg1 = ROOT.TLegend(0.7,.8,.9,.9)
leg1.AddEntry(roc_t1,"Timing","L")
leg1.AddEntry(roc_nt1,"No Timing","L")
leg1.Draw()

mg1.GetXaxis().SetTitle("Efficiency")
mg1.GetYaxis().SetTitle("Rejection")

canvas['c'+str(count)].SaveAs('roc1.pdf')
canvas['c'+str(count)].SaveAs('roc1.root')

count +=1

##eta2##
canvas['c'+str(count)] = ROOT.TCanvas()
roc_t2 = ROOT.TGraph(3,et2,ft2)
roc_nt2 = ROOT.TGraph(3,ent2,fnt2)
roc_nt2.SetLineColor(2)
roc_t2.SetLineColor(4)
roc_t2.SetMarkerColor(1)
roc_nt2.SetMarkerColor(1)
roc_t2.SetMarkerSize(2)
roc_nt2.SetMarkerSize(2)
roc_t2.SetMarkerStyle(3)
roc_nt2.SetMarkerStyle(3)

mg2=ROOT.TMultiGraph()
mg2.Add(roc_nt2)
mg2.Add(roc_t2)
mg2.Draw("ACP")

leg2 = ROOT.TLegend(0.7,.8,.9,.9)
leg2.AddEntry(roc_t2,"Timing","L")
leg2.AddEntry(roc_nt2,"No Timing","L")
leg2.Draw()

mg2.GetXaxis().SetTitle("Efficiency")
mg2.GetYaxis().SetTitle("Rejection")

canvas['c'+str(count)].SaveAs('roc2.pdf')
canvas['c'+str(count)].SaveAs('roc2.root')

count +=1


##eta3##
canvas['c'+str(count)] = ROOT.TCanvas()
roc_t3 = ROOT.TGraph(3,et3,ft3)
roc_nt3 = ROOT.TGraph(3,ent3,fnt3)
roc_nt3.SetLineColor(2)
roc_t3.SetLineColor(4)
roc_t3.SetMarkerColor(1)
roc_nt3.SetMarkerColor(1)
roc_t3.SetMarkerSize(2)
roc_nt3.SetMarkerSize(2)
roc_t3.SetMarkerStyle(3)
roc_nt3.SetMarkerStyle(3)

mg3=ROOT.TMultiGraph()
mg3.Add(roc_nt3)
mg3.Add(roc_t3)
mg3.Draw("ACP")

leg3 = ROOT.TLegend(0.7,.8,.9,.9)
leg3.AddEntry(roc_t3,"Timing","L")
leg3.AddEntry(roc_nt3,"No Timing","L")
leg3.Draw()

mg3.GetXaxis().SetTitle("Efficiency")
mg3.GetYaxis().SetTitle("Rejection")

canvas['c'+str(count)].SaveAs('roc3.pdf')
canvas['c'+str(count)].SaveAs('roc3.root')

count +=1


##tot##
canvas['c'+str(count)] = ROOT.TCanvas()
roc_ttot = ROOT.TGraph(3,ettot,fttot)
roc_nttot = ROOT.TGraph(3,enttot,fnttot)
roc_nttot.SetLineColor(2)
roc_ttot.SetLineColor(4)
roc_ttot.SetMarkerColor(1)
roc_nttot.SetMarkerColor(1)
roc_ttot.SetMarkerSize(2)
roc_nttot.SetMarkerSize(2)
roc_ttot.SetMarkerStyle(3)
roc_nttot.SetMarkerStyle(3)

mgtot=ROOT.TMultiGraph()
mgtot.Add(roc_nttot)
mgtot.Add(roc_ttot)
mgtot.Draw("ACP")

legtot = ROOT.TLegend(0.7,.8,.9,.9)
legtot.AddEntry(roc_ttot,"Timing","L")
legtot.AddEntry(roc_nttot,"No Timing","L")
legtot.Draw()

mgtot.GetXaxis().SetTitle("Efficiency")
mgtot.GetYaxis().SetTitle("Rejection")

canvas['c'+str(count)].SaveAs('roctot.pdf')
canvas['c'+str(count)].SaveAs('roctot.root')

count +=1
input("Press Enter to continue...")
