from DataFormats.FWLite import Events, Handle
from ROOT import TH1F, TCanvas, TLorentzVector, TFile
import sys
import math
#import CommonTools.CandAlgos
#import PhysicsTools.HepMCCandAlgos
#import FWCore.ParameterSet.Config as cms
from array import array

inputfiles = ['/afs/cern.ch/user/r/rselvati/wpigamma/CMSSW_8_0_24/src/StandardModel/WPiGamma/test/WPiGamma_pythia8_MINIAOD_1.root']

events = Events(inputfiles)

nevents = 1
maxEvents = 250
PV = 0
pTpi = 0
pTgamma = 0
pTe = 0
pTmu = 0
deta = 0
dphi = 0
etapi = 0
phipi = 0
deltaRMax = 10000
deltapTMax = 10000
pTmuMin = -1000
pion = 0
pion1 = 0
cont = 0
cont2 = 0
mu_selection = 0
mu_selection_event = 0
checker = None
checker1 = None
checker2 = None
cheker3 = None
checker4 = None
checker5 = None
checker6 = None
mu_ID = 0
gen_ID = 0
gen_mother = 0
eTphMin = -1000
pTpiMin = -1000
ph_cont = 0
ph_cont1 = 0
ph_cont2 = 0
candidate_ph = TLorentzVector()
candidate_pi = TLorentzVector()
candidate_eta = 0
candidate_phi = 0
photon_eta = 0
photon_phi = 0
pxpi = 0
pypi = 0
pzpi = 0
pxph = 0
pyph = 0
pzph = 0
pTph = 0
pi_from_w = 0
photon_from_w = 0
pi_and_photon_from_w = 0
counter = 0
counter1 = 0
track_iso = 0
ecal_iso = 0
hcal_iso = 0
calo_iso = 0
iso_sum = 0

genParticles_h = Handle ("vector<reco::GenParticle>")
genParticles_l = ('prunedGenParticles')

PFCandidates_h = Handle ("vector<pat::PackedCandidate>")
PFCandidates_l = ("packedPFCandidates")

slimmedMuons_h = Handle ("vector<pat::Muon>")
slimmedMuons_l = ("slimmedMuons")

slimmedPhotons_h = Handle ("vector<pat::Photon>")
slimmedPhotons_l = ("slimmedPhotons")



inv_mass_1 = TH1F("Invariant mass1", "W mass1", 200,0,120)
inv_mass_2 = TH1F("Invariant mass2", "W mass2", 200,0,120)
track_iso_hist = TH1F("Track iso", "Track isolation", 200,0,60)
ecal_iso_hist = TH1F("Ecal iso", "Ecal isolation", 200,0,40)
hcal_iso_hist = TH1F("Hcal iso", "Hcal isolation", 200,0,25)
calo_iso_hist = TH1F("Calo iso", "Calo isolation", 200,0,140)
iso_sum_hist = TH1F("Sum iso", "Sum isolation", 150,0,3)


for event in events:

    event.getByLabel (genParticles_l,genParticles_h)

    genParticles = genParticles_h.product()

    event.getByLabel (PFCandidates_l,PFCandidates_h)

    PFCandidates = PFCandidates_h.product()
    
    event.getByLabel (slimmedMuons_l,slimmedMuons_h)

    slimmedMuons = slimmedMuons_h.product()

    event.getByLabel (slimmedPhotons_l,slimmedPhotons_h)
    
    slimmedPhotons = slimmedPhotons_h.product()
    

    print "Event n: ", nevents, "/", maxEvents

    pTmuMin = -1000
    checker = False
    mu_ID = 0
    for mu in slimmedMuons:
        if mu.pt()>24 and mu.pt() > pTmuMin and mu.isMediumMuon()==True:
            pTmuMin = mu.pt()
            #print "pT :", mu.pt(), "Eta: ", mu.eta(), "phi:", mu.phi()
            #print "dxy: ", mu.innerTrack().dxy(), "dz: ", mu.innerTrack().dz()

            #if mu.isMediumMuon()==True:# and (mu.pfIsolationR04().sumChargedHadronPt + max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt)/mu.pt()) < 0.15:

                #PV = mu.originalObjectRef().vertexRef()
            mu_ID = mu.pdgId()
            mu_selection += 1 #checking how many mu over total pass the selection
            if checker == False: #checking how many events contain mu passing the if selection
                mu_selection_event += 1
            checker = True
    
    #matchMap = cms.EDProducer("MCTruthDeltaRMatcher", src = cms.InputTag("PFCandidates"), matched = cms.InputTag("genParticles"), distMin = cms.double(0.3))
    #print PFCandidates.data_[matchMap.map_.first].pt()
    cont1 = 0
    pTpiMax = -1000
    checker1 = False
    checker3 = False #together with checker4, used to calculate how many times both reconstructed pi and gamma come from generated particles whose mother is a W
    checker5 = False
    checker6 = False
    counter = 0
    for cand in PFCandidates:
        if cand.pdgId()*mu_ID == 211*13 and checker == True and cand.pt()>=20 and cand.trackHighPurity()==True and cand.pt()>pTpiMax and cand.fromPV()==3:
            pTpiMax = cand.pt()
            candidate_pi = cand.p4()
            #print "candidate vertex: ",cand.vertexRef()
            #print "pions we like: ", cont1
            #print "pdgId: ", cand.pdgId(), "muId: ", mu_ID, "pTpi: ", cand.pt(), "pxpi: ", cand.px(), "pypi: ", cand.py(), "pzpi: ", cand.pz(), "pienergy: ", cand.energy()

            pxpi = cand.px()
            pypi = cand.py()
            pzpi = cand.pz()
            pTpi = cand.pt()
            candidate_eta = cand.eta()
            candidate_phi = cand.phi()
            
            deltapTMax = 10000
            deltaRMax = 0.3
            gen_mother = 0
            gen_px = 0
            gen_py = 0
            gen_pz = 0
            gen_pt = 0

            #gen_ID = cand.pdgId()#so if it doesn't enter the following loop-and-if, it doesn't display "reconstructed particle is different..." either 
            for gen in genParticles: #matching candidate for W reconstruction with MC truth
                deltaR = math.sqrt((candidate_eta-gen.eta())*(candidate_eta-gen.eta())+(candidate_phi-gen.phi())*(candidate_phi-gen.phi()))
                deltapT = math.fabs(cand.pt()-gen.pt())

                if deltaR <= deltaRMax and deltapT < deltapTMax:
                    deltapTMax = deltapT
                    gen_ID = gen.pdgId()
                    gen_mother = gen.mother().pdgId()
                    gen_px = gen.px()
                    gen_py = gen.py()
                    gen_pz = gen.pz()
                    gen_pt = gen.pt()

            if not gen_ID == 211 and not gen_ID == -211:
                print "|||corresponding generated particle is not a pion, but: ", gen_ID
                if counter == 0:
                    checker6 = True #avoiding to count the case in which only a single reco-pion, not matching with a gen-pion, passes the selection

            if gen_ID == 211 or gen_ID == -211:
                print "\\\identity of generated PION's mother: ", gen_mother
                print "reco px: ", pxpi, "reco py: ", pypi, "reco pz: ", pzpi, "reco pT: ", pTpi
                print "gen px: ", gen_px, "gen py: ", gen_py, "gen pz: ", gen_pz, "gen pT: ", gen_pt
                if gen_mother == 24 or gen_mother == -24:
                    pi_from_w += 1
                    checker3 = True
                    gen_mother = 0
                    gen_ID = 0
                    counter += 1
                    if cont1 >= 1:
                        pi_from_w += 1
                        checker5 = True #activates in case the first reco particle to fulfill the conditions is not a pion from a W, while the second is

            if counter >= 1 and cont1 >= 1:
                pi_from_w = pi_from_w-1 #without this, a case in which the first particle to pass the selection is a pion but the second and final is not, is counted as a good reconstruction 
                checker3 = False
            if checker5 == True:
                checker3 = True
                
            cont1 += 1
            checker1 = True

    if not cont1==0:
        cont2 += 1

    cont += cont1

    ph_cont1 = 0
    eTphMin = -1000
    checker2 = False
    checker4 = False #together with checker3, used to calculate how many times both reconstructed pi and gamma come from generated particles whose mother is a W
    for photon in slimmedPhotons:
        if photon.et()>20 and photon.et()>eTphMin and checker1 == True:# and  photon.caloIso()+photon.trackIso()<5:
            eTphMin = photon.et()
            print "px: ", photon.px(), "py: ", photon.py(), "pz: ", photon.pz(), "energy: ", photon.energy()
            candidate_ph = photon.p4(0)
            checker2 = True
            ph_cont1 +=1
            
            pxph = photon.px()
            pyph = photon.py()
            pzph = photon.pz()
            pTph = photon.pt()
            photon_eta = photon.eta()
            photon_phi = photon.phi()
            track_iso = photon.trackIso()
            ecal_iso = photon.ecalIso()
            hcal_iso = photon.hcalIso()
            calo_iso = photon.caloIso()
            iso_sum = (track_iso+calo_iso)/photon.et()

            track_iso_hist.Fill(track_iso)
            ecal_iso_hist.Fill(ecal_iso)
            hcal_iso_hist.Fill(hcal_iso)
            calo_iso_hist.Fill(calo_iso)
            iso_sum_hist.Fill(iso_sum)
            
            
            gen_px = 0
            gen_py = 0
            gen_pz = 0
            gen_pt = 0
            deltapTMax = 10000
            deltaRMax = 0.3

            #gen_ID = cand.pdgId()#so if it doesn't enter the following loop-and-if, it doesn't display "reconstructed particle is different..." either 
            for gen in genParticles:
                deltaR = math.sqrt((photon_eta-gen.eta())*(photon_eta-gen.eta())+(photon_phi-gen.phi())*(photon_phi-gen.phi()))
                deltapT = math.fabs(photon.pt()-gen.pt())

                if deltaR <= deltaRMax and deltapT < deltapTMax :
                    deltapTMax = deltapT
                    gen_ID = gen.pdgId()
                    gen_mother = gen.mother().pdgId()
                    gen_px = gen.px()
                    gen_py = gen.py()
                    gen_pz = gen.pz()
                    gen_pt = gen.pt()
                
            if not gen_ID == 22:
                print "|||corresponding generated particle is not a photon, but: ", gen_ID
                print "______reco px: ", pxph, "reco py: ", pyph, "reco pz: ", pzph, "reco pT: ", pTph
                print "gen px: ", gen_px, "gen py: ", gen_py, "gen pz: ", gen_pz, "gen pT: ", gen_pt
            
            if gen_ID == 22:
                print "\\\identity of generated photon's mother: ", gen_mother
                print "______reco px: ", pxph, "reco py: ", pyph, "reco pz: ", pzph, "reco pT: ", pTph
                print "gen px: ", gen_px, "gen py: ", gen_py, "gen pz: ", gen_pz, "gen pT: ", gen_pt
                if gen_mother == 24 or gen_mother == -24:
                    photon_from_w += 1
                    checker4 = True

    if checker1 == True and checker2 == False and checker6 == False:
        counter1 += 1

    if not ph_cont1 == 0:
        ph_cont2 += 1 #incrementing the number of events with at least one photon

    ph_cont += ph_cont1

    
    if checker1 == True and checker2 == True:
        print "INV MASS: ", (candidate_ph + candidate_pi).M(), "costheta: ", (pxpi*pxph+pypi*pyph+pzpi*pzph)/(math.sqrt(pxpi**2+pypi**2+pzpi**2)*math.sqrt(pxph**2+pyph**2+pzph**2))
        if checker3 == False or checker4 == False:
            inv_mass_1.SetLineColor(3)
            inv_mass_1.Fill((candidate_ph + candidate_pi).M())
        if checker3 == True and checker4 == True:
            pi_and_photon_from_w += 1
            inv_mass_2.SetLineColor(2)
            inv_mass_2.Fill((candidate_ph + candidate_pi).M())


    print " "

    if nevents>=maxEvents :
        break

    nevents += 1

print "n of mu passing selection: ", mu_selection, "|| n of events with mu passing selection: ", mu_selection_event
print "n of pi: ", cont, "|| n of events with at least one pi: ", cont2#, "|| n of pi per event: ", float(cont)/128
print "n of photons: ", ph_cont, "|| n of events with at least one photon: ", ph_cont2
print "number of matching reconstructed-generated pi, originating from W: ", pi_from_w-counter1
print "number of matching reconstructed-generated photon, originating from W: ", photon_from_w
print "both pi and photon come frome W: ", pi_and_photon_from_w

#newfile = TFile("mass_and_isolation.root","recreate")
#canvas1 = TCanvas()
#inv_mass_2.Draw("hist")
#inv_mass_1.Draw("SAME,hist")
#canvas1.Print("Wmass.png")
#canvas1.Write()

#canvas2 = TCanvas()
#track_iso_hist.Draw("hist")
#canvas2.Print("TrackIso.png")
#track_iso_hist.Write()

#canvas3 = TCanvas()
#ecal_iso_hist.Draw("hist")
#canvas3.Print("EcalIso.png")
#ecal_iso_hist.Write()

#canvas4 = TCanvas()
#hcal_iso_hist.Draw("hist")
#canvas4.Print("HcalIso.png")
#hcal_iso_hist.Write()

#canvas5 = TCanvas()
#calo_iso_hist.Draw("hist")
#canvas5.Print("CaloIso.png")
#calo_iso_hist.Write()

#canvas6 = TCanvas()
#iso_sum_hist.Draw("hist")
#canvas6.Print("IsoSum.png")
#iso_sum_hist.Write()




    



