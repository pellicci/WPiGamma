from DataFormats.FWLite import Events, Handle
from ROOT import TH1D, TCanvas, TLorentzVector
import sys
import math
from array import array

inputfiles = ['/afs/cern.ch/user/r/rselvati/wpigamma/CMSSW_8_0_24/test/WPiGamma_pythia8_MINIAOD_1.root']

events = Events(inputfiles)

nevents = 1
maxEvents = 250
wpi = 0
wgamma = 0
wpos = 0
wneg = 0
welec = 0
wmu = 0
wnue = 0
wnumu = 0
pTpi = 0
pTgamma = 0
pTe = 0
pTmu = 0
deta = 0
dphi = 0
etapi = 0
phipi = 0
deltaRMax = 10000
pTpiMax = 1000
pTmuMax = 1000
difference1 = 0
difference2 = 0
difference3 = 0
deltaRcomp = 0
identity = 0
pTcomp = 0
pion = 0
pion1 = 0
PV = 0
PV_list = []
cont = 0
cont2 = 0
muoni = 0
checker = None

genParticles_h = Handle ("vector<reco::GenParticle>")
genParticles_l = ('prunedGenParticles')

PFCandidates_h = Handle ("vector<pat::PackedCandidate>")
PFCandidates_l = ("packedPFCandidates")

slimmedMuons_h = Handle ("vector<pat::Muon>")
slimmedMuons_l = ("slimmedMuons")

histo1 = TH1D("pTpi","pT of pi", 100,0,200)
histo2 = TH1D("pTgamma","pT of gamma", 100,0,200)
histo3 = TH1D("pTe","pT of e", 70,0,200)
histo4 = TH1D("pTmu","pT of mu", 70,0,200)
histo5 = TH1D("deltaR pi-gamma", "deltaR",100,0,15)

for event in events:

    event.getByLabel (genParticles_l,genParticles_h)

    genParticles = genParticles_h.product()

    event.getByLabel (PFCandidates_l,PFCandidates_h)

    PFCandidates = PFCandidates_h.product()
    
    event.getByLabel (slimmedMuons_l,slimmedMuons_h)

    slimmedMuons = slimmedMuons_h.product()
    

    print "Event n: ", nevents, "/", maxEvents

    for gpart in genParticles:
         if gpart.pdgId()==24 and gpart.numberOfDaughters()==2 : #W+
             wpos += 1
             
             #print "W+) Daughter n1: ", gpart.daughter(0).pdgId(), "and daughter n2: ", gpart.daughter(1).pdgId()
             if gpart.daughter(0).pdgId()==211 or gpart.daughter(1).pdgId()==211 :
                 wpi+=1
                 if gpart.daughter(0).pdgId()==211 : 
                     pTpi = gpart.daughter(0).pt()
                     #print "pTpi generated: ", pTpi
                     etapi = gpart.daughter(0).eta()
                     phipi = gpart.daughter(0).phi()
                 if gpart.daughter(1).pdgId()==211 :
                     pTpi = gpart.daughter(1).pt()
                     #print "pTpi generated: ", pTpi
                     etapi = gpart.daughter(1).eta()
                     phipi = gpart.daughter(1).phi()


                 for cand in PFCandidates:
                     difference1 = cand.phi()-phipi
                     difference2 = cand.eta()-etapi
                     deltaRcomp = math.sqrt(difference1*difference1 + difference2*difference2) #deltaR between gen pi and compared reco particle
                     if deltaRcomp < deltaRMax :
                         deltaRMax = deltaRcomp
                         identity = cand.pdgId() #identity of the particle with closest eta to the generated pi
                         pTcomp = cand.pt() #pT of this particle, to be compared with the sample pi's pT
                 if identity == 211:
                     pion += 1 #number of particles which are pi closest to the gen pi
                 #print "closest phi particle ID: ", identity, "and its pT: ", pTcomp
                 deltaRMax = 10000

                 
                 #for cand in PFCandidates:
                 #    difference3 = math.fabs(cand.pt()-pTpi)
                 #    if difference3 < pTpiMax :
                 #        pTpiMax = difference3
                 #        identity = cand.pdgId() #identity of the particle with closest eta to the generated pi
                 #        pTcomp = cand.pt() #pT of this particle, to be compared with the sample pi's pT
                 #if identity == 211:
                 #    pion1 += 1 #number of particles which are pi closest to the gen pi
                 #print "closest phi particle ID/////: ", identity, "and its pT/////: ", pTcomp
                 #pTpiMax = 1000

             
             histo1.Fill(pTpi) #filling histo1 with pT of pi

             if gpart.daughter(0).pdgId()==22 or gpart.daughter(1).pdgId()==22 :
                 wgamma += 1
             if gpart.daughter(0).pdgId()==22 :
                 pTgamma = gpart.daughter(0).pt()
             if gpart.daughter(1).pdgId()==22 :
                 pTgamma = gpart.daughter(1).pt()
             histo2.Fill(pTgamma) #filling histo2 with pT of gamma
             if gpart.daughter(0).pdgId()==211 and gpart.daughter(1).pdgId()==22 :
                 deta = gpart.daughter(0).eta()-gpart.daughter(1).eta()
                 dphi = gpart.daughter(0).phi()-gpart.daughter(1).phi()
                 deltaR = math.sqrt(deta*deta + dphi*dphi)
             if gpart.daughter(0).pdgId()==22 and gpart.daughter(1).pdgId()==211 :
                 deta = gpart.daughter(1).eta()-gpart.daughter(0).eta()
                 dphi = gpart.daughter(1).phi()-gpart.daughter(0).phi()
                 deltaR = math.sqrt(deta*deta + dphi*dphi)
             histo5.Fill(deltaR)
                 


         if gpart.pdgId()==-24 and gpart.numberOfDaughters()==2 : #W-
             wneg += 1
           
             #print "W-) Daughter n1: ", gpart.daughter(0).pdgId(), "and daughter n2: ", gpart.daughter(1).pdgId()
             if gpart.daughter(0).pdgId()==11 or gpart.daughter(1).pdgId()==11 :
                 welec += 1
             if gpart.daughter(0).pdgId()==11 : #filling histo3 with pT of e-
                 pTe = gpart.daughter(0).pt()
             if gpart.daughter(1).pdgId()==11 :
                 pTe = gpart.daughter(1).pt()
             histo3.Fill(pTe)
             if gpart.daughter(0).pdgId()==-12 or gpart.daughter(1).pdgId()==-12 :
                 wnue += 1
             if gpart.daughter(0).pdgId()==13 or gpart.daughter(1).pdgId()==13 :
                 wmu += 1
             if gpart.daughter(0).pdgId()==13 : #filling histo4 with pT of mu-
                 pTmu = gpart.daughter(0).pt()
             if gpart.daughter(1).pdgId()==13 :
                 pTmu = gpart.daughter(1).pt()
             histo4.Fill(pTmu)
             if gpart.daughter(0).pdgId()==-14 or gpart.daughter(1).pdgId()==-14 :
                 wnumu += 1
                 
    pTmuMax = 1000
    checker = False
    for mu in slimmedMuons:
        if mu.pt()>24 and mu.pt() < pTmuMax:
            pTmuMax = mu.pt()
            #print "is global muon: ", mu.isGlobalMuon()
            #print mu.originalObjectRef().pvAssociationQuality()
            quality_list = [1,2,3,5,6,7]
            print "pT :", mu.pt(), "Eta: ", mu.eta(), "phi:", mu.phi()
            print "dxy: ", mu.innerTrack().dxy(), "dz: ", mu.innerTrack().dz()
            # mu.originalObjectRef().pvAssociationQuality() in quality_list and
            if mu.isMediumMuon()==True and (mu.pfIsolationR04().sumChargedHadronPt + max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt)/mu.pt()) < 0.15 and mu.dB() < 0.2 :
                PV = mu.originalObjectRef().vertexRef()
                checker = True
                
    cont1 = 0
    for cand in PFCandidates:
        if cand.vertexRef()==PV and cand.pt()>=20 and checker==True and cand.trackHighPurity()==True:
            cont1 += 1
            print "pions we like: ", cont1
    if not cont1==0:
        cont2 += 1

    cont += cont1

    #for cand in PFCandidates:
    #    if cand.pdgId()==13 and cand.pt() < pTmuMax:
    #        pTmuMax = cand.pt()
    #        PV = cand.vertexRef()
    #        PV_list.append(PV)
    #pTmuMax = 1000
    #for element in PV_list:
    #    for cand in PFCandidates:
    #        if cand.vertexRef()==element and cand.pdgId()==211:
    #            cont += 1
    #print "cont: ", cont
    #cont = 0

    if nevents>=maxEvents :
        break

    nevents += 1

print "pion number: ", pion, "pion number using pT matching: ", pion1
print "ntot charged particles: ", cont, "n of events with at least one pi: ", cont2, "n of pi per event: ", float(cont)/128
#if not welec==wnue :
#    print "Number of electrons =! from number of anti-nue"
#if not wmu==wnumu :
#    print "Number of muons =! from number of anti-numu"

#print "Number of W+ dacaying into 2 particles: ", wpos, "and number of W- decaying into 2 particles: ", wneg
#print "Number of pi+ from  W+: ", wpi, " and number of gamma from W+: ", wgamma
#print "Number of e- and anti-nu couples from W-: ", welec, "and number of mu and anti-numu from W-: ", wmu
    

#canvas1 = TCanvas()
#histo1.Draw("hist")
#canvas1.Print ("pTpi.png")
#canvas2 = TCanvas()
#histo2.Draw("hist")
#canvas2.Print("pTgamma.png")
#canvas3 = TCanvas()
#histo3.Draw("hist")
#canvas3.Print("pTe.png")
#canvas4 = TCanvas()
#histo4.Draw("hist")
#canvas4.Print("pTmu.png")
#canvas5 = TCanvas()
#histo5.Draw("hist")
#canvas5.Print("deltaRpigamma.png")



