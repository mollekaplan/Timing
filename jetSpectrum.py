import ROOT 
import os,sys, math
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser

parser=OptionParser()
#parse
opts,args=parser.parse_args()
minjetpt = 30.
minjeteta=0.
maxjeteta = 2.5
maxevt = 10000

def deltaphi(a,b):
    x = b-a
    if abs(x) <= math.pi:
        return x
    while x> math.pi: x-= 2*math.pi
    while x< -math.pi: x+= 2*math.pi
    return x

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
            drmin = dr
            matched = genjet
            
    return matched

def matchgenjetloose(jet, genjets):
    drmin = float("inf")
    matched = None
    for genjet in genjets:
        dr = deltar(jet,genjet) 
        if dr<0.4 and dr<drmin:
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

def mjjcalc(jets,muons):
    mjj = 0.
    jetslist = []
    for jet1 in jets:
        if abs(jet1.eta()) > maxjeteta or abs(jet1.eta()) < minjeteta: continue
        if jet1.pt()< minjetpt : continue
        if deltar(jet1,muons[0])<0.4: continue
        if deltar(jet1,muons[1])<0.4: continue
	for jet2 in jetslist:
            j1 = ROOT.TLorentzVector()
            j1.SetPtEtaPhiM ( jet1.pt(), jet1.eta(), jet1.phi(), jet1.mass() )
            j2 = ROOT.TLorentzVector()
            j2.SetPtEtaPhiM ( jet2.pt(),jet2.eta(),jet2.phi(),jet2.mass() )
            jj = j1+j2
            mjj_i = jj.M()
            if mjj_i>mjj: mjj = mjj_i
	jetslist.append(jet1)
    return mjj

def mjjcalc2(jets1,jets2, muons,eta1,eta2,eta3,eta4):
    mjj = 0.
    jetslist = []
    for jet1 in jets1 or jets2:
	if jet1 in jets1:
            if abs(jet1.eta()) > eta1 or abs(jet1.eta()) < eta2: continue
	    if jet1.pt()< minjetpt : continue
            if deltar(jet1,muons[0])<0.4: continue
            if deltar(jet1,muons[1])<0.4: continue
	elif jet1 in jets2: 
	    if abs(jet1.eta()) > eta3 or abs(jet1.eta()) < eta4: continue
 	else: print "not in list!"
	
        for jet2 in jetslist:
            j1 = ROOT.TLorentzVector()
            j1.SetPtEtaPhiM ( jet1.pt(), jet1.eta(), jet1.phi(), jet1.mass() )
            j2 = ROOT.TLorentzVector()
            j2.SetPtEtaPhiM ( jet2.pt(),jet2.eta(),jet2.phi(),jet2.mass() )
            jj = j1+j2
            mjj_i = jj.M()
            if mjj_i>mjj: mjj = mjj_i
        jetslist.append(jet1)
    return mjj


## DY timing
nfiles=-1

files={}
indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Sep5jme/dymm200timing/MINIAODSIM/'
#dymm200timing/MINIAODSIM/step3_dymm200timing_MINIAODSIM_119.root
#/eos/cms/store/group/upgrade/timing/pfintegration/Sep5jme/dymm200timing/MINIAODSIM/step3_dymm200timing_MINIAODSIM_93.root
files["DY200_timing"] = [indir + filename for filename in os.listdir(indir)[0:nfiles] if (filename != "step3_dymm200timing_MINIAODSIM_119.root" and filename !="step3_dymm200timing_MINIAODSIM_93.root")]

indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Sep5jme/dymm200notiming/MINIAODSIM/'
files["DY200_notiming"] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Aug26jme/dymmnoputiming/MINIAODSIM/'
files["DYnopu_timing"] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]

indir = '/eos/cms/store/group/upgrade/timing/pfintegration/Aug26jme/dymmnopunotiming/MINIAODSIM/'
files["DYnopu_notiming"] = [indir + filename for filename in os.listdir(indir)[0:nfiles]]


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

#print "TODO"

histos={}
jetsdict={}
def loop(inputfiles,label="",timingbool=False):
    print "-> Doing",label, inputfiles[0]
    global histos
    global jetsdict
    events = Events (inputFiles=inputfiles)
    
    handlejetpuppi = Handle("std::vector<pat::Jet>")
    labeljetpuppi = "slimmedJetsPuppi"
    
    handlepfcand = Handle("vector<pat::PackedCandidate>")
    labelpfcand = ("packedPFCandidates")

    handlechsjet = Handle("std::vector<pat::Jet>")
    labelchsjet = "slimmedJets"

    handleinfo = Handle("GenEventInfoProduct")
    labelinfo = "generator"

    handlegenjet = Handle("std::vector<reco::GenJet>")
    labelgenjet = "slimmedGenJets"

    handleVtx = Handle("vector<reco::Vertex>")
    labelVtx = ("offlineSlimmedPrimaryVertices")

    histos["tot_"+ label]=ROOT.TH1D("tot_"+label,label,1,0,1)
    histos["tot_w_"+ label]=ROOT.TH1D("tot_w_"+label,label,1,0,1)

    histos["njets_gen_"+ label]=ROOT.TH1D("njets_gen_"+label,label,20,0,20)
    histos["njets_gen_w_"+ label]=ROOT.TH1D("njets_gen_w_"+label,label,20,0,20)

    histos["nvtx_"+ label]=ROOT.TH1D("nvtx_"+label,label,500,0,500)

    histos["njets_"+ label]=ROOT.TH1D("njets_"+label,label,20,0,20)
    histos["njets_w_"+ label]=ROOT.TH1D("njets_w_"+label,label,20,0,20)
    histos["njets_chs_"+ label]=ROOT.TH1D("njets_chs_"+label,label,20,0,20)
    histos["njets_chs_w_"+ label]=ROOT.TH1D("njets_chs_w_"+label,label,20,0,20)

    histos["njets_pu_"+ label]=ROOT.TH1D("njets_pu_"+label,label,20,0,20)
    histos["njets_pu_w_"+ label]=ROOT.TH1D("njets_pu_w_"+label,label,20,0,20)
    histos["njets_matched_"+ label]=ROOT.TH1D("njets_matched_"+label,label,20,0,20)
    histos["njets_matched_w_"+ label]=ROOT.TH1D("njets_matched_w_"+label,label,20,0,20)

    histos["njets_chs_pu_"+ label]=ROOT.TH1D("njets_chs_pu_"+label,label,20,0,20)
    histos["njets_chs_pu_w_"+ label]=ROOT.TH1D("njets_chs_pu_w_"+label,label,20,0,20)
    histos["njets_chs_matched_"+ label]=ROOT.TH1D("njets_chs_matched_"+label,label,20,0,20)
    histos["njets_chs_matched_w_"+ label]=ROOT.TH1D("njets_chs_matched_w_"+label,label,20,0,20)

    histos["njets2d_w_"+ label]=ROOT.TH2D("njets2d_w_"+label,label,20,0,20,20,0,20)

    histos["mjj_"+label]=ROOT.TH1D("mjj_"+label,label,10,0,200)
    histos["mjj25_"+label]=ROOT.TH1D("mjj25_"+label,label,10,0,200) 
    histos["mjj47_"+label]=ROOT.TH1D("mjj47_"+label,label,10,0,200) 
    histos["mjjgen_"+label]=ROOT.TH1D("mjjgen_"+label,label,10,0,200)

    jetsdict["puppi_"+label] = []

    print(histos)

    debug=0
    for ievt, event in enumerate(events):
    #for event in events:
        event.getByLabel (labeljetpuppi, handlejetpuppi)
        event.getByLabel (labelpfcand, handlepfcand)
        event.getByLabel (labelinfo, handleinfo)
        event.getByLabel (labelchsjet, handlechsjet)
        event.getByLabel (labelgenjet, handlegenjet)
        event.getByLabel (labelVtx, handleVtx)
        
        pfcands = handlepfcand.product()
        puppijets = handlejetpuppi.product()
        chsjets = handlechsjet.product()
        genjets = handlegenjet.product()
        info = handleinfo.product()

        histos["tot_"+label].Fill(0.5)
        histos["tot_w_"+label].Fill(0.5,info.weight())

        histos["nvtx_"+label].Fill(handleVtx.product().size())

        nmuons = 0
        muons = []
        for cand in pfcands:
            if abs(cand.pdgId()) != 13:
                continue
            if cand.pt()<20.:
                continue
            if not (cand.pvAssociationQuality()>=1 and cand.vertexRef().key()==0):
                continue
            nmuons += 1
            muons.append(cand)
            
        if nmuons < 2:
            continue
        
	njets_gen=0
        genjetsColl=[]
        for jet in genjets:
            if jet.pt()< minjetpt : continue
            if abs(jet.eta()) > maxjeteta or abs(jet.eta()) < minjeteta: continue
            if deltar(jet,muons[0])<0.4: continue
            if deltar(jet,muons[1])<0.4: continue
            njets_gen +=1
            genjetsColl.append( jet )
        ## FILL HERE
        histos["njets_gen_"+label].Fill(njets_gen)
        histos["njets_gen_w_"+label].Fill(njets_gen,info.weight())
    
        ## I have two muons
        njets=0
        njets_matched=0
        njets_pu=0
       
        puppijetsColl=[]
        for jet in puppijets:
	    if jet.pt()< minjetpt : continue
            if deltar(jet,muons[0])<0.4: continue
            if deltar(jet,muons[1])<0.4: continue
    	    jetsdict["puppi_"+label].append(jet)
	    if abs(jet.eta()) > maxjeteta or abs(jet.eta()) < minjeteta: continue
	    njets +=1
            if ispujet(jet,genjets): njets_pu += 1
            if matchgenjetloose(jet, genjets): njets_matched+=1
            puppijetsColl.append( (jet,ispujet(jet,genjets),matchgenjetloose(jet, genjets)!=None))

        if debug <10:
            print "---------------------------------"
            print "mu:    ", ';'.join(["(%.1f,%.2f,%.2f)"%(mu.pt(),mu.eta(),mu.phi()) for mu in muons])
            print "genj  :", ';'.join(["(%.1f,%.2f,%.2f)"%(j.pt(),j.eta(),j.phi()) for j in genjetsColl])
            print "genjA :", ';'.join(["(%.1f,%.2f,%.2f)"%(j.pt(),j.eta(),j.phi()) for j in genjets])
            print "puppi :", ';'.join(["(%.1f,%.2f,%.2f,%d,%d)"%(j.pt(),j.eta(),j.phi(),int(ispu),int(ismatched)) for j,ispu,ismatched in puppijetsColl])
            print "---------------------------------"
            debug +=1
            
        ## FILL HERE
	redo=0
	try:
            histos["njets_"+label].Fill(njets)
            histos["njets_w_"+label].Fill(njets,info.weight())
            histos["njets_pu_"+label].Fill(njets_pu)
            histos["njets_pu_w_"+label].Fill(njets_pu,info.weight())
            histos["njets_matched_"+label].Fill(njets_matched)
            histos["njets_matched_w_"+label].Fill(njets_matched,info.weight())
            histos["njets2d_w_"+ label].Fill(njets_gen,njets,info.weight())
	except:
	    print "puppi jets died :("
	    redo=1
	    break
	    

        njets_chs=0
        njets_chs_matched=0
        njets_chs_pu=0
        for jet in chsjets:
            if jet.pt()< minjetpt : continue
            if abs(jet.eta()) > maxjeteta or abs(jet.eta()) < minjeteta: continue
            if deltar(jet,muons[0])<0.4: continue
            if deltar(jet,muons[1])<0.4: continue
            njets_chs +=1
            if ispujet(jet,genjets): njets_chs_pu += 1
            if matchgenjetloose(jet, genjets): njets_chs_matched+=1

	redo = 0
	try:
	    histos["njets_chs_"+label].Fill(njets_chs)
            histos["njets_chs_w_"+label].Fill(njets_chs,info.weight())
            histos["njets_chs_pu_"+label].Fill(njets_chs_pu)
            histos["njets_chs_pu_w_"+label].Fill(njets_chs_pu,info.weight())
            histos["njets_chs_matched_"+label].Fill(njets_chs_matched)
            histos["njets_chs_matched_w_"+label].Fill(njets_chs_matched,info.weight())
	except: 
	    print "chs jets died :("
            redo=1
            break

	##invariant mass calculation##
	mjj_puppi = mjjcalc(puppijets,muons)
	if mjj_puppi > 0.: histos["mjj_"+label].Fill(mjj_puppi)

	if timingbool:
	    mjj_puppi25 = mjjcalc2(puppijets,jetsdict["puppi_pu200_notiming"],muons,2.5,0.,4.7,2.5)
	    if mjj_puppi25 > 0.: histos["mjj25_"+label].Fill(mjj_puppi25)
	    mjj_puppi47 = mjjcalc2(puppijets,jetsdict["puppi_pu200_notiming"],muons,4.7,2.5,2.5,0.)
            if mjj_puppi47 > 0.: histos["mjj47_"+label].Fill(mjj_puppi47)

	mjj_gen = mjjcalc(genjets,muons)
	if mjj_gen > 0.: histos["mjjgen_"+label].Fill(mjj_gen)

        if maxevt > 0 and ievt > maxevt: break

    if redo==1:
	redo = 0	
	loop(inputfiles,label)

if __name__ == "__main__": 
    out=ROOT.TFile.Open("njets.root","RECREATE")
    loop(files["DY200_notiming"],"pu200_notiming")
    #loop(files["DYnopu_timing"],"nopu_timing")
    loop(files["DY200_timing"],"pu200_timing",True)
    #loop(files["DYnopu_notiming"],"nopu_notiming")

    out.cd()
    for hStr in histos:
        histos[hStr].Write()

    print "-> DONE"
