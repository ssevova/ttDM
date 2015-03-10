//================================================================================================
//
// Perform preselection semi-leptonic ttbar events following B2G-14-004 and produce bacon bits
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => output bacon bits file name
//   argv[3] => dataset type: "mcsig", "mcttbar", "mcv0jets", "mcbkg", "e", "mu"
//   argv[4] => JSON file for run-lumi filtering of data, specify "none" for MC or no filtering
//   argv[5] => cross section (pb), ignored for data
//________________________________________________________________________________________________

// bacon object headers
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"

// JSON file parser
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"

#include "DMSAna/Utils/interface/mt2w_bisect.h"

// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TMath.h>

// Other C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

// Jet corrections
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"


//=== FUNCTION DECLARATIONS ======================================================================================

// Object selections
bool passJetSel(const baconhep::TJet *jet);
bool passEleSel(const baconhep::TElectron *electron, const double rho);
bool passMuonSel(const baconhep::TMuon *muon);
bool hasMuonOverlap(const baconhep::TElectron *electron, const TClonesArray *muonArr);

// check if MC event is V + 0-jets		   
bool isV0Jets(const TClonesArray* genParArr);

float computeMt(const float pt, const float phi, const float met, const float metphi);
float computeMT2W(const TLorentzVector &lep, 
                  const std::vector<TLorentzVector> &bjetsv,
                  const std::vector<TLorentzVector> &jetsv, 
                  const float met, const float metphi);
float computeTTbarCorr(const float pt);
double effArea(const double eta);


//=== MAIN =======================================================================================================

int main(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings and constants
  //==============================================================================================================  
  
  // handle input arguments
  const std::string infilename   = argv[1];
  const std::string outfilename  = argv[2];
  const std::string dstypename   = argv[3];
  const std::string jsonfilename = argv[4];
  const double      xsec         = atof(argv[5]);
    
  // determine dataset type
  enum {
    kMCSig=1,  // MC signal
    kMCTTbar,  // MC ttbar
    kMCV0Jets, // MC W/Z + Jets (0-jets events only extracted from inclusive sample)
    kMCBkg,    // Other MC background
    kEle,      // single electron data
    kMu        // single muon data
  };

  unsigned int dstype=0;
  if     (dstypename.compare("mcsig")==0)    { dstype = kMCSig; } 
  else if(dstypename.compare("mcttbar")==0)  { dstype = kMCTTbar; } 
  else if(dstypename.compare("mcv0jets")==0) { dstype = kMCV0Jets; }
  else if(dstypename.compare("mcbkg")==0)    { dstype = kMCBkg; } 
  else if(dstypename.compare("e")==0)        { dstype = kEle; }
  else if(dstypename.compare("mu")==0)       { dstype = kMu;  }
  assert(dstype>0);

  const std::string cmssw_base = getenv("CMSSW_BASE");
  
  // Trigger bits mapping file
  std::string trigfilename = cmssw_base + std::string("/src/BaconAna/DataFormats/data/HLTFile_v1");
  baconhep::TTrigger trigger(trigfilename);

  // JEC correction files
  std::vector<JetCorrectorParameters> corrParams;
  if(dstype==kEle || dstype==kMu) {
    corrParams.push_back(JetCorrectorParameters( (cmssw_base + std::string("/src/DMSAna/Utils/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).c_str() ));
    corrParams.push_back(JetCorrectorParameters( (cmssw_base + std::string("/src/DMSAna/Utils/data/Summer13_V1_DATA_L2Relative_AK5PF.txt")).c_str() ));
    corrParams.push_back(JetCorrectorParameters( (cmssw_base + std::string("/src/DMSAna/Utils/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).c_str() ));
    corrParams.push_back(JetCorrectorParameters( (cmssw_base + std::string("/src/DMSAna/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).c_str() ));

  } else {
    corrParams.push_back(JetCorrectorParameters( (cmssw_base + std::string("/src/DMSAna/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt")).c_str() ));
    corrParams.push_back(JetCorrectorParameters( (cmssw_base + std::string("/src/DMSAna/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt")).c_str() ));
    corrParams.push_back(JetCorrectorParameters( (cmssw_base + std::string("/src/DMSAna/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt")).c_str() ));
  }

  //
  // Cuts
  //
  const unsigned int NJETS_CUT  = 3;
  const unsigned int NBJETS_CUT = 0; 
  const double MET_CUT          = 120;
  const double ELE_PT_CUT       = 30;
  const double ELE_ETA_CUT      = 2.5;
  const double MUON_PT_CUT      = 30;
  const double MUON_ETA_CUT     = 2.1;
  const double JET_PT_CUT       = 30;
  const double JET_ETA_CUT      = 4.0;
  const double BTAG_WP          = 0.679;  // CSV medium working point
  
  //
  // constants
  //
  const double ELE_MASS  = 0.000511;
  const double MUON_MASS = 0.105658369;
  
  const int ELE_PDGID  = 11;  // e-
  const int MUON_PDGID = 13;  // mu-
  const int TOP_PDGID  = 6;   // t
  
  const double ECAL_GAP_LOW  = 1.4442;
  const double ECAL_GAP_HIGH = 1.566;
  
  
  // Print summary of selection cuts
  std::cout << " ===== Cuts ===== " << std::endl;
  std::cout << " -- MET > " << MET_CUT << std::endl;
  std::cout << " -- Electron:";
  std::cout << " pT > "     << ELE_PT_CUT;
  std::cout << ", |eta| < " << ELE_ETA_CUT << " (gap excl.)" << std::endl;
  std::cout << " -- Muon:    ";
  std::cout << " pT > "     << MUON_PT_CUT;
  std::cout << ", |eta| < " << MUON_ETA_CUT << std::endl;
  std::cout << " -- Jet:     ";
  std::cout << " pT > "     << JET_PT_CUT;
  std::cout << ", |eta| < " << JET_ETA_CUT << std::endl; 
  std::cout << " -- Number of jets   >= " << NJETS_CUT << std::endl;
  std::cout << " -- Number of b-tags >= " << NBJETS_CUT << std::endl; 
  std::cout << " -- b-tag def: CSV   >  " << BTAG_WP << std::endl;
  std::cout << std::endl;

  
  //--------------------------------------------------------------------------------------------------------------
  // Set up output file
  //==============================================================================================================
  unsigned int runNum, lumiSec, evtNum;          // event ID
  unsigned int metfilter;                        // MET filter bits
  unsigned int npv, npu;                         // number of PV / PU
  unsigned int njets, nbjets;                    // jet multiplicity
  float rho;                                     // event energy density
  float scale1fb;                                // cross section scale factor per 1/fb
  float evtWeight;                               // event weight (NOT from cross section normalization)
  float pfmet, pfmetphi, pfmt;                   // PF MET
  float mT2W;
  
  int lepId;                                     // lepton PDG ID
  TLorentzVector *lepton=0;                      // lepton 4-vector
  
  //
  // Jet ordering:
  // jet1 => leading jet
  // jet2 => sub-leading jet
  //
  float           jet1flavAlgo,  jet2flavAlgo;  // jet flavor, algorithmic definition (MC only)
  float           jet1flavPhys,  jet2flavPhys;  // jet flavor, physics definition (MC only)
  TLorentzVector *jet1vec=0,    *jet2vec=0;     // jet 4-vector

  TFile *outFile = new TFile(outfilename.c_str(),"RECREATE");
  TH1D hTotalEvents("TotalEvents","TotalEvents",1,-10,10);
  TTree *outTree = new TTree("Events","Events");
  
  outTree->Branch("runNum",    &runNum,    "runNum/i");
  outTree->Branch("lumiSec",   &lumiSec,   "lumiSec/i");
  outTree->Branch("evtNum",    &evtNum,    "evtNum/i");
  outTree->Branch("metfilter", &metfilter, "metfilter/i");
  outTree->Branch("npv",       &npv,       "npv/i");
  outTree->Branch("npu",       &npu,       "npu/i");
  outTree->Branch("njets",     &njets,     "njets/i");
  outTree->Branch("nbjets",    &nbjets,    "nbjets/i");
  outTree->Branch("rho",       &rho,       "rho/F");
  outTree->Branch("scale1fb",  &scale1fb,  "scale1fb/F");
  outTree->Branch("evtWeight", &evtWeight, "evtWeight/F");
  outTree->Branch("pfmet",     &pfmet,     "pfmet/F");
  outTree->Branch("pfmetphi",  &pfmetphi,  "pfmetphi/F");
  outTree->Branch("pfmt",      &pfmt,      "pfmt/F");
  outTree->Branch("mT2W",      &mT2W,      "mT2W/F");
  
  // lepton variables
  outTree->Branch("lepId", &lepId, "lepId/I");
  outTree->Branch("lepton", "TLorentzVector", &lepton);
  
  // jet 1 variables
  outTree->Branch("jet1flavAlgo", &jet1flavAlgo, "jet1flavAlgo/F");
  outTree->Branch("jet1flavPhys", &jet1flavPhys, "jet1flavPhys/F");
  outTree->Branch("jet1",     "TLorentzVector", &jet1vec);

  // jet 2 variables
  outTree->Branch("jet2flavAlgo", &jet2flavAlgo, "jet2flavAlgo/F");
  outTree->Branch("jet2flavPhys", &jet2flavPhys, "jet2flavPhys/F");
  outTree->Branch("jet2",     "TLorentzVector", &jet2vec);


  //--------------------------------------------------------------------------------------------------------------
  // Process input file
  //==============================================================================================================
  baconhep::TEventInfo *info = 0; TBranch *infoBr   = 0;
  TClonesArray *eleArr       = 0; TBranch *eleBr    = 0;
  TClonesArray *muonArr      = 0; TBranch *muonBr   = 0;
  TClonesArray *jetArr       = 0; TBranch *jetBr    = 0;
  TClonesArray *pvArr        = 0; TBranch *pvBr     = 0;
  TClonesArray *genParArr    = 0; TBranch *genParBr = 0;
    
  std::cout << "Processing " << infilename << "..." << std::endl;    
  TFile *infile    = TFile::Open(infilename.c_str()); assert(infile);
  TTree *eventTree = (TTree*)infile->Get("Events");   assert(eventTree);
  
  hTotalEvents.Add((TH1D*)infile->Get("TotalEvents"));

  eventTree->SetBranchAddress("Info",     &info,    &infoBr);
  eventTree->SetBranchAddress("Electron", &eleArr,  &eleBr);
  eventTree->SetBranchAddress("Muon",     &muonArr, &muonBr);
  eventTree->SetBranchAddress("Jet05",    &jetArr,  &jetBr);
  eventTree->SetBranchAddress("PV",       &pvArr,   &pvBr);
  if(dstype!=kEle && dstype!=kMu) { eventTree->SetBranchAddress("GenParticle", &genParArr, &genParBr); }
  
  // Set up object to handle good run-lumi filtering if necessary
  bool hasJSON = false;
  baconhep::RunLumiRangeMap rlrm;
  if(jsonfilename.compare("none")!=0) {
    rlrm.addJSONFile(jsonfilename);
    hasJSON = true;
  }

  // Set up JEC
  FactorizedJetCorrector jetcorr(corrParams);

  //
  // loop over events
  //
  for(unsigned int ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    infoBr->GetEntry(ientry);
    
    //
    // MET requirement
    //
    if(info->pfMETC < MET_CUT) continue;

    double xsWgt  = 1;  // event weight for cross section normalization
    double evtWgt = 1;  // event weight for other corrections
 
    if(dstype==kEle || dstype==kMu) {
      if(hasJSON) {
        // JSON filter
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(!rlrm.hasRunLumi(rl)) continue;
      }
    } else {
      xsWgt = 1000.*xsec/hTotalEvents.GetEntries();
    }
  
    //
    // primary vertex requirement
    //
    if(!(info->hasGoodPV)) continue;
  
    //
    // trigger requirement
    // Note: to avoid overlaps in data, only the electron (muon) trigger is required in the electron (muon) dataset 
    //
    const std::string eleTrigName = "HLT_Ele27_WP80_v*";
    const std::string muTrigName  = "HLT_IsoMu24_eta2p1_v*";
    bool passTrigger = false;
    if     (dstype==kEle) { passTrigger = trigger.pass(eleTrigName, info->triggerBits); }
    else if(dstype==kMu)  { passTrigger = trigger.pass(muTrigName,  info->triggerBits); }
    else                  { passTrigger = trigger.pass(eleTrigName, info->triggerBits) || trigger.pass(muTrigName, info->triggerBits); }
    if(!passTrigger) continue;
    
    
    if(dstype==kMCV0Jets) {
      //
      // check for V + 0-jet events
      //     
      genParArr->Clear();
      genParBr->GetEntry(ientry);
      if(!isV0Jets(genParArr)) continue;
    } 
    
    //
    // lepton selection: only one good electron or muon in the event
    //
    unsigned int nLeptons=0;
    int lid=0;
    TLorentzVector vLep;
    
    eleArr->Clear();
    eleBr->GetEntry(ientry);

    muonArr->Clear();
    muonBr->GetEntry(ientry);
    
    for(int i=0; i<muonArr->GetEntriesFast(); i++) {
      const baconhep::TMuon *muon = (baconhep::TMuon*)muonArr->At(i);
    
      if(muon->pt	 <= MUON_PT_CUT)  continue;
      if(fabs(muon->eta) >= MUON_ETA_CUT) continue;
      if(!passMuonSel(muon))	          continue;    

      ///// Good muon found! /////

      nLeptons++;
      if(nLeptons>1) break;
      
      if(dstype==kEle) continue;  // don't accept muons in the electron dataset
      
      // trigger object matching
      if(!trigger.passObj(muTrigName, muon->hltMatchBits)) continue;

      lid = (muon->q < 0) ? MUON_PDGID : -MUON_PDGID;
      vLep.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);
    }
    if(nLeptons>1) continue;
            
    for(int i=0; i<eleArr->GetEntriesFast(); i++) {
      const baconhep::TElectron *electron = (baconhep::TElectron*)eleArr->At(i);
  
      // exclude ECAL region around barrel-endcap transition
      if(fabs(electron->scEta)>=ECAL_GAP_LOW && fabs(electron->scEta)<=ECAL_GAP_HIGH) continue;
    
      if(electron->pt	     <= ELE_PT_CUT)   continue;
      if(fabs(electron->eta) >= ELE_ETA_CUT)  continue;
      if(!passEleSel(electron, info->rhoIso)) continue;
      
      // muon overlap cleaning
      if(hasMuonOverlap(electron, muonArr)) continue;

      ///// Good electron found! /////

      nLeptons++;
      if(nLeptons>1) break;
      
      if(dstype==kMu) continue;  // don't accept electrons in the muon dataset
          
      // trigger object matching
      if(!trigger.passObj(eleTrigName, electron->hltMatchBits)) continue;

      lid = (electron->q < 0) ? ELE_PDGID : -ELE_PDGID;
      vLep.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);
    }		 
    if(nLeptons>1) continue;
        
    if(lid==0) continue;  // no good leptons found
     
  
    //
    // jet selection
    //
    const baconhep::TJet *jet1=0, *jet2=0;
    TLorentzVector vJet1, vJet2;
    std::vector<TLorentzVector> goodJetsv, bJetsv;
    jetArr->Clear();
    jetBr->GetEntry(ientry);
    for(int i=0; i<jetArr->GetEntriesFast(); i++) {
      const baconhep::TJet *jet = (baconhep::TJet*)jetArr->At(i);

      TLorentzVector vRawJet;
      vRawJet.SetPtEtaPhiM(jet->ptRaw, jet->eta, jet->phi, (jet->mass)*(jet->ptRaw)/(jet->pt));

      jetcorr.setJetPt(vRawJet.Pt());
      jetcorr.setJetEta(vRawJet.Eta());
      jetcorr.setJetPhi(vRawJet.Phi());
      jetcorr.setJetE(vRawJet.E());
      jetcorr.setRho(info->rhoJet);
      jetcorr.setJetA(jet->area);
      jetcorr.setJetEMF(-99.0);
      double jec = jetcorr.getCorrection();

      TLorentzVector vJet;
      vJet.SetPtEtaPhiM(vRawJet.Pt()*jec, vRawJet.Eta(), vRawJet.Phi(), vRawJet.M()*jec);

      // lepton cleaning
      if(vJet.DeltaR(vLep) < 0.5) continue;

      if(vJet.Pt()        <= JET_PT_CUT)  continue;
      if(fabs(vJet.Eta()) >= JET_ETA_CUT) continue;
      if(!passJetSel(jet))                continue;
      
      if(fabs(jet->eta) < 2.4 && jet->csv > BTAG_WP) {
        bJetsv.push_back(vJet);
      } else {
        goodJetsv.push_back(vJet);
      }

      if(!jet1 || vJet.Pt() > vJet1.Pt()) {
        jet2  = jet1;
        jet1  = jet;
        vJet2 = vJet1;
        vJet1 = vJet;

      } else if(!jet2 || vJet.Pt() > vJet2.Pt()) {
        jet2  = jet;
        vJet2 = vJet;
      }
    }
    if(goodJetsv.size() + bJetsv.size() < NJETS_CUT) continue;
//    if(bJetsv.size() < NBJETS_CUT) continue;
  
  
    ////////// === Event passes preselection! HURRAY! === //////////
    
    // Read in PV branch to get number of primary vertices
    pvArr->Clear();
    pvBr->GetEntry(ientry);
    
    //
    // compute ttbar MC pT correction
    // Note: are cap at pT(top) = 400 GeV and the factor of 1.001 are standard prescriptions,
    //       or just for B2G-14-004? 
    //
    if(dstype==kMCTTbar) {
      genParArr->Clear();
      genParBr->GetEntry(ientry);
      
      double w1=1, w2=1;
      for(int i=0; i<genParArr->GetEntriesFast(); i++) {
        const baconhep::TGenParticle *p = (baconhep::TGenParticle*)genParArr->At(i);
        if(p->pdgId ==  TOP_PDGID && p->status == 3) { w1 = computeTTbarCorr(TMath::Min((float)400.,p->pt)); }
        if(p->pdgId == -TOP_PDGID && p->status == 3) { w2 = computeTTbarCorr(TMath::Min((float)400.,p->pt)); }
      }
      evtWgt *= 1.001*sqrt(w1*w2);
    }
  
    //
    // Fill output tree
    //
    runNum    = info->runNum;
    lumiSec   = info->lumiSec;
    evtNum    = info->evtNum;
    metfilter = info->metFilterFailBits;
    npv       = pvArr->GetEntriesFast();
    npu       = info->nPUmean;
    njets     = goodJetsv.size() + bJetsv.size();
    nbjets    = bJetsv.size();
    rho       = info->rhoJet;
    scale1fb  = xsWgt;
    evtWeight = evtWgt;
    pfmet     = info->pfMETC;
    pfmetphi  = info->pfMETCphi;
    pfmt      = computeMt(vLep.Pt(), vLep.Phi(), info->pfMET, info->pfMETphi);
    mT2W      = computeMT2W(vLep, bJetsv, goodJetsv, info->pfMET, info->pfMETphi);

    lepId  = lid;
    lepton = &vLep;
    
    jet1vec      = &vJet1;
    jet1flavAlgo = jet1 ? jet1->mcFlavor     : -999;
    jet1flavPhys = jet1 ? jet1->mcFlavorPhys : -999;
    
    jet2vec      = &vJet2;
    jet2flavAlgo = jet2 ? jet2->mcFlavor     : -999;
    jet2flavPhys = jet2 ? jet2->mcFlavorPhys : -999;

    outTree->Fill();
  }
  
  delete infile;
  infile=0, eventTree=0;
  
  outFile->Write();
  outFile->Close();
  
  delete info;
  delete eleArr;
  delete muonArr;
  delete jetArr;
  delete pvArr;
  delete genParArr;
  
  return 0;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
bool passJetSel(const baconhep::TJet *jet)
{
  // Loose PFJet ID
  if(jet->neuHadFrac >= 0.99) return false;
  if(jet->neuEmFrac  >= 0.99) return false;
  if(jet->nParticles <= 1)    return false;  
  if(fabs(jet->eta)<2.4) {
    if(jet->chHadFrac == 0)    return false;
    if(jet->nCharged  == 0)    return false;
    if(jet->chEmFrac  >= 0.99) return false;  
  }
  return true;
}

//--------------------------------------------------------------------------------------------------
bool passEleSel(const baconhep::TElectron *electron, const double rho)
{
  if(!(electron->typeBits & baconhep::kEcalDriven)) return false;
  //if(!(electron->typeBits & baconhep::kTrackerDriven)) return false;
  
  // cut-based tight ID (https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification#Electron_ID_Working_Points)
  // and PF isolation (https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation#Example_for_electrons)
  
  double iso = electron->chHadIso03 + TMath::Max( 0.0,(electron->gammaIso03 + electron->neuHadIso03 - rho*effArea(electron->scEta)) );
  if(iso >= 0.1*(electron->pt)) return false;
  
  if(fabs(electron->scEta) <= 1.479) {
    if(fabs(electron->dEtaIn)       >= 0.004)                       return false;
    if(fabs(electron->dPhiIn)       >= 0.03)                        return false;
    if(electron->sieie              >= 0.01)                        return false;
    if(electron->hovere             >= 0.12)                        return false;
    if(fabs(electron->d0)           >= 0.02)                        return false;
    if(fabs(electron->dz)           >= 0.1)                         return false;
    if(fabs(1.0 - electron->eoverp) >= 0.05*(electron->ecalEnergy)) return false;
    if(electron->isConv)                                            return false;
    if(electron->nMissingHits       > 0)                            return false;
    
  } else {
    if(fabs(electron->dEtaIn)       >= 0.005)                       return false;
    if(fabs(electron->dPhiIn)       >= 0.02)                        return false;
    if(electron->sieie              >= 0.03)                        return false;
    if(electron->hovere             >= 0.10)                        return false;
    if(fabs(electron->d0)           >= 0.02)                        return false;
    if(fabs(electron->dz)           >= 0.1)                         return false;
    if(fabs(1.0 - electron->eoverp) >= 0.05*(electron->ecalEnergy)) return false;
    if(electron->isConv)                                            return false;
    if(electron->nMissingHits       > 0)                            return false;
  }
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool passMuonSel(const baconhep::TMuon *muon)
{
  // Tight muon ID (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon)
  if(!(muon->typeBits & baconhep::kGlobal)) return false;
  if(!(muon->typeBits & baconhep::kPFMuon)) return false;
  if(muon->muNchi2    >= 10)  return false;
  if(muon->nValidHits <  1)   return false;
  if(muon->nMatchStn  <  2)   return false;
  if(fabs(muon->d0)   >= 0.2) return false;
  if(fabs(muon->dz)   >= 0.5) return false;
  if(muon->nPixHits   <  1)   return false;
  if(muon->nTkLayers  <  6)   return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso04 + TMath::Max(muon->neuHadIso04 + muon->gammaIso04 - 0.5*(muon->puIso04), double(0));
  if(iso >= 0.12*(muon->pt)) return false;
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool hasMuonOverlap(const baconhep::TElectron *electron, const TClonesArray *muonArr)
{
  const double ELE_MASS  = 0.000511;
  const double MUON_MASS = 0.105658369;
  
  bool hasOverlap = false;
  TLorentzVector vEle;
  vEle.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);
  for(int j=0; j<muonArr->GetEntriesFast(); j++) {
    const baconhep::TMuon *muon = (baconhep::TMuon*)muonArr->At(j);
    if(!(muon->typeBits & baconhep::kGlobal) && !(muon->typeBits & baconhep::kTracker)) continue;
    if(!(muon->typeBits & baconhep::kPFMuon)) continue;
    if(muon->pt <= 15) continue;
    
    TLorentzVector vMu;
    vMu.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);
    if(vMu.DeltaR(vEle) < 0.1) {
      hasOverlap = true;
      break;
    }
  }
  
  return hasOverlap;
}

//--------------------------------------------------------------------------------------------------
bool isV0Jets(const TClonesArray* genParArr)
{
  const int W_PDGID      = 24;  // W+
  const int Z_PDGID      = 23;  // Z
  const int TOP_PDGID    = 6;   // t  
  const int GLUON_PDGID  = 21;  // gluon
    
  assert(genParArr);
  
  const baconhep::TGenParticle *vpar=0;
  for(int j=0; j<genParArr->GetEntriesFast(); j++) { 
    const baconhep::TGenParticle *genp = (baconhep::TGenParticle*)genParArr->At(j);
    
    // Assumes boson is listed first and subsequent partons listed after
    if(genp->status==3 && ( abs(genp->pdgId)==W_PDGID || abs(genp->pdgId)==Z_PDGID )) { 
      vpar = genp;
    }
    
    // partons for V + N-jets are listed to have the same parent as the W
    if(vpar) { 
      if(genp->status==3 
         && ( abs(genp->pdgId)<TOP_PDGID || genp->pdgId==GLUON_PDGID )
	 && (genp->parent == vpar->parent)
      ) {
        return false;
      }
    }
  }
  
  return (vpar!=0);
}

//--------------------------------------------------------------------------------------------------
float computeMt(const float pt, const float phi, const float met, const float metphi)
{
  TLorentzVector vLep; vLep.SetPtEtaPhiM(pt,0,phi,0);
  TLorentzVector vMet; vMet.SetPtEtaPhiM(met,0,metphi,0);
  double dphi = vLep.DeltaPhi(vMet);
  return sqrt( 2.0* pt * met *(1.0-cos(dphi)) );
}

//--------------------------------------------------------------------------------------------------
float computeMT2W(const TLorentzVector &lep,
                  const std::vector<TLorentzVector> &bjetsv,
                  const std::vector<TLorentzVector> &jetsv, 
                  const float met, const float metphi)
{
  double minval = 8000.;

  const double upper_bound = 8000.;
  const double error_value = upper_bound-1.0;
  const double scan_step   = 0.5;
  mt2w_bisect::mt2w mt2w_event(upper_bound, error_value, scan_step);
  double pmiss[3] = { 0., met*cos(metphi), met*sin(metphi) };
  double pl[4] = { lep.E(), lep.Px(), lep.Py(), lep.Pz() };

  if(bjetsv.size()>=2) {
    for(unsigned int b1=0; b1<bjetsv.size(); b1++) {
      double pb1[4] = { bjetsv[b1].E(), bjetsv[b1].Px(), bjetsv[b1].Py(), bjetsv[b1].Pz() };

      for(unsigned int b2=0; b2<bjetsv.size(); b2++) {
        if(b1==b2) continue;
        double pb2[4] = { bjetsv[b2].E(), bjetsv[b2].Px(), bjetsv[b2].Py(), bjetsv[b2].Pz() };
 
        mt2w_event.set_momenta(pl,pb1,pb2,pmiss);
        double val = mt2w_event.get_mt2w();
        if(val < minval)
          minval = val;
      }
    }

  } else if(bjetsv.size()==1) {
    double pb1[4] = { bjetsv[0].E(), bjetsv[0].Px(), bjetsv[0].Py(), bjetsv[0].Pz() };
    for(unsigned int j=0; j<jetsv.size(); j++) {
      double pb2[4] = { jetsv[j].E(), jetsv[j].Px(), jetsv[j].Py(), jetsv[j].Pz() };

      mt2w_event.set_momenta(pl,pb1,pb2,pmiss);
      double val = mt2w_event.get_mt2w();
      if(val < minval)
        minval = val;
    }
  }

  return minval;
}

//--------------------------------------------------------------------------------------------------
float computeTTbarCorr(const float pt)
{
  // Reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#MC_SFs_Reweighting
  return exp(0.156 - 0.00137*pt);
}

//--------------------------------------------------------------------------------------------------
double effArea(const double eta)
{
  // effective area for PU correction (type==kEleGammaAndNeutralHadronIso03 in EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h)

  if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.130; }
  else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.137; }
  else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.067; }
  else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.089; }
  else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.107; }
  else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.110; }
  else                                             { return 0.138; }
}
