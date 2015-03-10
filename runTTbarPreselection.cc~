//================================================================================================
//
// Perform preselection semi-leptonic ttbar events and produce bacon bits
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => output bacon bits file name
//   argv[3] => dataset type: "mcsig", "mcttbar", "mcbkg", "e", "mu"
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


//=== FUNCTION DECLARATIONS ======================================================================================

bool passJetID(const baconhep::TJet *jet);
bool passEleID(const baconhep::TElectron *electron, const double rho);
bool passMuonID(const baconhep::TMuon *muon);
bool hasMuonOverlap(const baconhep::TElectron *electron, const TClonesArray *muonArr);
// finds and stores W daughters and b's
// [0] - [2]: daughters from W- and bbar
// [3] - [5]: daughters from W+ and b
// assumes input vectors are empty
void findDaughters(const TClonesArray* genParArr, 
                   std::vector<TLorentzVector> &dvecs,
		   std::vector<const baconhep::TGenParticle*> &dpars);

float computeMt(const float pt, const float phi, const float met, const float metphi);
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
    kMCBkg,    // Other MC background
    kEle,      // single electron data
    kMu        // single muon data
  };

  unsigned int dstype=0;
  if     (dstypename.compare("mcsig")==0)   { dstype = kMCSig; } 
  else if(dstypename.compare("mcttbar")==0) { dstype = kMCTTbar; } 
  else if(dstypename.compare("mcbkg")==0)   { dstype = kMCSig; } 
  else if(dstypename.compare("e")==0)       { dstype = kEle; }
  else if(dstypename.compare("mu")==0)      { dstype = kMu;  }
  assert(dstype>0);
  
  // Trigger bits mapping file
  std::string trigfilename = getenv("CMSSW_BASE");
  trigfilename += "/src/BaconAna/DataFormats/data/HLTFile_v0";
  baconhep::TTrigger trigger(trigfilename);
  
  // Cuts
  const unsigned int NJETS_CUT  = 4;
  // const unsigned int NBJETS_CUT = 0;
  const double ELE_PT_CUT       = 30;
  const double ELE_ETA_CUT      = 2.5;
  const double MUON_PT_CUT      = 30;
  const double MUON_ETA_CUT     = 2.1;
  const double JET_PT_CUT       = 30;
  const double JET_ETA_CUT      = 4.7;
  const double BTAG_WP          = 0.679;  // CSV medium working point
//  const double ANTIBTAG_WP      = 0.244;  // b-veto using CSV loose working point
  const double MATCHING_MAX_DR  = 0.3;
  const double WMASSLOW         = 40;
  const double WMASSHIGH        = 130;

  //
  // constants
  //
  const double ELE_MASS  =  0.000511;
  const double MUON_MASS =  0.105658369;
  const double W_MASS    = 80.4;
  
  const int ELE_PDGID    = 11;  // e-
  const int MUON_PDGID   = 13;  // mu-
  const int BOTTOM_PDGID = 5;   // b
  const int TOP_PDGID    = 6;   // t
  
  const double ECAL_GAP_LOW  = 1.4442;
  const double ECAL_GAP_HIGH = 1.566;
  
  
  // Print summary of selection cuts
  std::cout << " ===== Cuts ===== " << std::endl;
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
  // std::cout << " -- Number of b-tags >= " << NBJETS_CUT << std::endl; 
  // std::cout << " -- b-tag def: CSV   >  " << BTAG_WP << std::endl;
  std::cout << " -- MC matching: dR  <  " << MATCHING_MAX_DR << std::endl;
  std::cout << " -- dijet mass window:  " << WMASSLOW << " < Mjj < " << WMASSHIGH << std::endl;
  std::cout << std::endl;

  
  //--------------------------------------------------------------------------------------------------------------
  // Set up output file
  //==============================================================================================================
  unsigned int runNum, lumiSec, evtNum;                                      // event ID
  unsigned int metfilter;                                                    // MET filter bits
  unsigned int npv, npu;                                                     // number of PV / PU
  unsigned int njets, nbjets;                                                // jet multiplicity
  float        rho;
  float        scale1fb;                                                     // cross section scale factor per 1/fb
  float        ttWeight;                                                     // ttbar event weight (needed for MadGraph ttbar)
  float        pfmet, pfmetphi, pfmt;                                        // PF MET
  float        mvamet, mvametphi, mvamt;                                     // MVA MET
  
  int             lepId;                                                     // lepton PDG ID
  TLorentzVector *lepton=0;                                                  // lepton 4-vector
  
  // Notes: bjet1 and bjet2 are the leading and sub-leading jets
  //        with the two highest b-tagger value. jet1 and jet2
  //        are pairs of jets which are not bjet1 and bjet2, with
  //        jet1 pT > jet2 pT.
  //
  int             jet1flavGen,  jet2flavGen,  jet3flavGen,  jet4flavGen;    // PDG ID of gen-level quarks matched to jets in ttbar events
  float           jet1qgid,     jet2qgid,     jet3qgid,     jet4qgid;      // jet q/g discriminant
  float           jet1csv,      jet2csv,      jet3csv,      jet4csv;       // jet b-tag value
  float           jet1pmass,    jet2pmass,    jet3pmass,    jet4pmass;     // jet pruned mass
  float           jet1q,        jet2q,        jet3q,        jet4q;         // jet charge
  float           jet1q03,      jet2q03,      jet3q03,      jet4q03;       // jet charge with kappa=0.3
  float           jet1q2,       jet2q2,       jet3q2,       jet4q2;        // jet charge with kappa=2
  float           jet1ptD,      jet2ptD,      jet3ptD,      jet4ptD;       // jet pT-D
  float           jet1flavAlgo, jet2flavAlgo, jet3flavAlgo, jet4flavAlgo;  // jet flavor, algorithmic definition (MC only)
  float           jet1flavPhys, jet2flavPhys, jet3flavPhys, jet4flavPhys;  // jet flavor, physics definition (MC only)
  unsigned int    jet1npar,     jet2npar,     jet3npar,     jet4npar;      // jet particle multiplicity
  TLorentzVector *jet1vec=0,   *jet2vec=0,   *jet3vec=0,   *jet4vec=0;     // jet 4-vector
  TVector2       *jet1pull=0,  *jet2pull=0,  *jet3pull=0,  *jet4pull=0;    // jet pull vector

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
  outTree->Branch("scale1fb",  &scale1fb,  "scale1fb/F");
  outTree->Branch("ttWeight",  &ttWeight,  "ttWeight/F");
  outTree->Branch("pfmet",     &pfmet,     "pfmet/F");
  outTree->Branch("pfmetphi",  &pfmetphi,  "pfmetphi/F");
  outTree->Branch("pfmt",      &pfmt,      "pfmt/F");
  outTree->Branch("mvamet",    &mvamet,    "mvamet/F");
  outTree->Branch("mvametphi", &mvametphi, "mvametphi/F");
  outTree->Branch("mvamt",     &mvamt,     "mvamt/F");
  
  // lepton variables
  outTree->Branch("lepId", &lepId, "lepId/I");
  outTree->Branch("lepton", "TLorentzVector", &lepton);
  
  // jet 1 variables
  outTree->Branch("jet1flavGen",  &jet1flavGen,  "jet1flavGen/I");
  outTree->Branch("jet1qgid",     &jet1qgid,     "jet1qgid/F");
  outTree->Branch("jet1csv",      &jet1csv,      "jet1csv/F");
  outTree->Branch("jet1pmass",    &jet1pmass,    "jet1pmass/F");
  outTree->Branch("jet1q",        &jet1q,        "jet1q/F");
  outTree->Branch("jet1q03",      &jet1q03,      "jet1q03/F");
  outTree->Branch("jet1q2",       &jet1q2,       "jet1q2/F");
  outTree->Branch("jet1ptD",      &jet1ptD,      "jet1ptD/F");
  outTree->Branch("jet1flavGen",  &jet1flavGen, "jet1flavGen/F");
  outTree->Branch("jet1flavAlgo", &jet1flavAlgo, "jet1flavAlgo/F");
  outTree->Branch("jet1flavPhys", &jet1flavPhys, "jet1flavPhys/F");  
  outTree->Branch("jet1npar",     &jet1npar,     "jet1npar/i");
  outTree->Branch("jet1",     "TLorentzVector", &jet1vec);
  outTree->Branch("jet1pull", "TVector2",       &jet1pull);

  // jet 2 variables
  outTree->Branch("jet2flavGen",  &jet2flavGen,  "jet2flavGen/I");
  outTree->Branch("jet2qgid",     &jet2qgid,     "jet2qgid/F");
  outTree->Branch("jet2csv",      &jet2csv,      "jet2csv/F");
  outTree->Branch("jet2pmass",    &jet2pmass,    "jet2pmass/F");
  outTree->Branch("jet2q",        &jet2q,        "jet2q/F");
  outTree->Branch("jet2q03",      &jet2q03,      "jet2q03/F");
  outTree->Branch("jet2q2",       &jet2q2,       "jet2q2/F");
  outTree->Branch("jet2ptD",      &jet2ptD,      "jet2ptD/F");
  outTree->Branch("jet2flavGen",  &jet2flavGen,  "jet2flavGen/F");
  outTree->Branch("jet2flavAlgo", &jet2flavAlgo, "jet2flavAlgo/F");
  outTree->Branch("jet2flavPhys", &jet2flavPhys, "jet2flavPhys/F");  
  outTree->Branch("jet2npar",     &jet2npar,     "jet2npar/i");
  outTree->Branch("jet2",     "TLorentzVector", &jet2vec);
  outTree->Branch("jet2pull", "TVector2",       &jet2pull);

  // jet 3 variables
  outTree->Branch("jet3flavGen",  &jet3flavGen,  "jet3flavGen/I");
  outTree->Branch("jet3qgid",     &jet3qgid,     "jet3qgid/F");
  outTree->Branch("jet3csv",      &jet3csv,      "jet3csv/F");
  outTree->Branch("jet3pmass",    &jet3pmass,    "jet3pmass/F");
  outTree->Branch("jet3q",        &jet3q,        "jet3q/F");
  outTree->Branch("jet3q03",      &jet3q03,      "jet3q03/F");
  outTree->Branch("jet3q2",       &jet3q2,       "jet3q2/F");
  outTree->Branch("jet3ptD",      &jet3ptD,      "jet3ptD/F");
  outTree->Branch("jet3flavGen",  &jet3flavGen,  "jet3flavGen/F");
  outTree->Branch("jet3flavAlgo", &jet3flavAlgo, "jet3flavAlgo/F");
  outTree->Branch("jet3flavPhys", &jet3flavPhys, "jet3flavPhys/F");  
  outTree->Branch("jet3npar",     &jet3npar,     "jet3npar/i");
  outTree->Branch("jet3",     "TLorentzVector", &jet3vec);
  outTree->Branch("jet3pull", "TVector2",       &jet3pull);

  // jet 4 variables
  outTree->Branch("jet4flavGen",  &jet4flavGen,  "jet4flavGen/I");
  outTree->Branch("jet4qgid",     &jet4qgid,     "jet4qgid/F");
  outTree->Branch("jet4csv",      &jet4csv,      "jet4csv/F");
  outTree->Branch("jet4pmass",    &jet4pmass,    "jet4pmass/F");
  outTree->Branch("jet4q",        &jet4q,        "jet4q/F");
  outTree->Branch("jet4q03",      &jet4q03,      "jet4q03/F");
  outTree->Branch("jet4q2",       &jet4q2,       "jet4q2/F");
  outTree->Branch("jet4ptD",      &jet4ptD,      "jet4ptD/F");
  outTree->Branch("jet4flavGen",  &jet4flavGen,  "jet4flavGen/F");
  outTree->Branch("jet4flavAlgo", &jet4flavAlgo, "jet4flavAlgo/F");
  outTree->Branch("jet4flavPhys", &jet4flavPhys, "jet4flavPhys/F");  
  outTree->Branch("jet4npar",     &jet4npar,     "jet4npar/i");
  outTree->Branch("jet4",     "TLorentzVector", &jet4vec);
  outTree->Branch("jet4pull", "TVector2",       &jet4pull);
 
  //--------------------------------------------------------------------------------------------------------------
  // Process input file
  //==============================================================================================================
  baconhep::TEventInfo *info = new baconhep::TEventInfo();
  TClonesArray *eleArr  = new TClonesArray("baconhep::TElectron");
  TClonesArray *muonArr = new TClonesArray("baconhep::TMuon");
  TClonesArray *jetArr  = new TClonesArray("baconhep::TJet");
  TClonesArray *pvArr   = new TClonesArray("baconhep::TVertex");
  TClonesArray *genParArr = 0;
  if(dstype==kMCTTbar) { genParArr = new TClonesArray("baconhep::TGenParticle"); }
    
  std::cout << "Processing " << infilename << "..." << std::endl;    
  TFile *infile    = TFile::Open(infilename.c_str()); assert(infile);
  TTree *eventTree = (TTree*)infile->Get("Events");   assert(eventTree);
  
  hTotalEvents.Add((TH1D*)infile->Get("TotalEvents"));

  eventTree->SetBranchAddress("Info",     &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &eleArr);  TBranch *eleBr  = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon",     &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Jet05",    &jetArr);  TBranch *jetBr  = eventTree->GetBranch("Jet05");
  eventTree->SetBranchAddress("PV",       &pvArr);   TBranch *pvBr   = eventTree->GetBranch("PV");
  TBranch *genParBr = 0;
  if(dstype==kMCTTbar) {
    eventTree->SetBranchAddress("GenParticle", &genParArr);
    genParBr = eventTree->GetBranch("GenParticle");
  }
  
  double weight = 1;  // event weight for cross section normalization
  
  // Set up object to handle good run-lumi filtering if necessary
  bool hasJSON = false;
  baconhep::RunLumiRangeMap rlrm;
  if(jsonfilename.compare("none")!=0) {
    rlrm.addJSONFile(jsonfilename);
    hasJSON = true;
  }
  
  //
  // loop over events
  //
  for(unsigned int ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    infoBr->GetEntry(ientry);
  
    if(dstype==kEle || dstype==kMu) {
      if(hasJSON) {
        // JSON filter
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(!rlrm.hasRunLumi(rl)) continue;
      }
    } else {
      weight = 1000.*xsec/hTotalEvents.GetEntries();
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
      if(!passMuonID(muon))	          continue;    

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
    
      if(electron->pt	     <= ELE_PT_CUT)  continue;
      if(fabs(electron->eta) >= ELE_ETA_CUT) continue;
      if(!passEleID(electron, info->rhoIso)) continue;
      
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
  
    std::vector<TLorentzVector> dvecs;
    std::vector<const baconhep::TGenParticle*> dpars;
    if(dstype==kMCTTbar) {     
      //
      // find hadronic top at generator level
      //
      genParArr->Clear();
      genParBr->GetEntry(ientry);      
      findDaughters(genParArr,dvecs,dpars);
      // temp fix! should not be needed with latest mother find fix in bacon...
      if(dvecs.size()!=6 || dpars.size()!=6) continue;
    }
    
    //
    // jet selection
    //
    std::vector<double> minDR;
    std::vector<const baconhep::TJet*> minDR_jet;
    for(unsigned int i=0; i<6; i++){
      minDR.push_back(999);
      minDR_jet.push_back(0);
    }
    std::vector<const baconhep::TJet*> goodJets;
    jetArr->Clear();
    jetBr->GetEntry(ientry);
    
    for(int i=0; i<jetArr->GetEntriesFast(); i++) {
      const baconhep::TJet *jet = (baconhep::TJet*)jetArr->At(i);
    
      // lepton cleaning
      TLorentzVector vJet; vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
      if(vJet.DeltaR(vLep) < 0.5) continue;

      if(jet->pt        <= JET_PT_CUT)  continue;
      if(fabs(jet->eta) >= JET_ETA_CUT) continue;
      if(!passJetID(jet))               continue;

      goodJets.push_back(jet);
      
      // match
      if(dstype==kMCTTbar) {
        for(int j=0; j<6; j++) {
          if(fabs(dpars[j]->pdgId) > BOTTOM_PDGID) continue;
	  
	  float tmpDR = dvecs[j].DeltaR(vJet);
          if(tmpDR < minDR[j] && tmpDR < MATCHING_MAX_DR) { 
            minDR[j]     = tmpDR;
            minDR_jet[j] = jet; 
          }
        }
      }
    }

    if(goodJets.size()  < NJETS_CUT)  continue;
  
  
    ////////// === Event passes preselection! HURRAY! === //////////    
         
    // Read in PV branch to get number of primary vertices
    pvArr->Clear();
    pvBr->GetEntry(ientry);
    
    //
    // compute ttbar MC pT correction
    // Note: are cap at pT(top) = 400 GeV and the factor of 1.001 are standard prescriptions,
    //       or just for B2G-14-004? 
    //
    double ttWgt = 1;
    if(dstype==kMCTTbar) {
      // genParArr->Clear();
      // genParBr->GetEntry(ientry);
      
      double w1=1, w2=1;
      for(int i=0; i<genParArr->GetEntriesFast(); i++) {
        const baconhep::TGenParticle *p = (baconhep::TGenParticle*)genParArr->At(i);
        if(p->pdgId ==  TOP_PDGID && p->status == 3) { w1 = computeTTbarCorr(TMath::Min((float)400.,p->pt)); }
        if(p->pdgId == -TOP_PDGID && p->status == 3) { w2 = computeTTbarCorr(TMath::Min((float)400.,p->pt)); }
      }
      ttWgt *= 1.001*sqrt(w1*w2);
    }
    
    
    //
    // Make dijet pairs from remaining jets
    //
    for(unsigned int i1=0; i1<goodJets.size(); i1++) {
      const baconhep::TJet *jet1 = goodJets[i1];
      
      int flavGen1 = 0;
      if(dstype==kMCTTbar) {
        for(int j=0; j<6; j++) {
          if(!minDR_jet[j]) continue;
	  if(jet1 == minDR_jet[j]) { flavGen1 = dpars[j]->pdgId; }
        }
      }

      for(unsigned int i2=i1+1; i2<goodJets.size(); i2++) {
        const baconhep::TJet *jet2 = goodJets[i2];
	
	int flavGen2=0;
	if(dstype==kMCTTbar) {
          for(int j=0; j<6; j++) {
            if(!minDR_jet[j]) continue;
	    if(jet2 == minDR_jet[j]) { flavGen2 = dpars[j]->pdgId; }
          }
        }

	for(unsigned int i3=i2+1; i3<goodJets.size(); i3++){
	  const baconhep::TJet *jet3 = goodJets[i3];

	  int flavGen3=0;
	  if(dstype==kMCTTbar) {
	    for(int j=0; j<6; j++) {
	      if(!minDR_jet[j]) continue;
	      if(jet3 == minDR_jet[j]) { flavGen3 = dpars[j]->pdgId; }
	    }
	  }

	  for(unsigned int i4=i3+1; i4<goodJets.size(); i4++){
	    const baconhep::TJet *jet4 = goodJets[i4];
           
	    int flavGen4=0;
	    if(dstype==kMCTTbar) {
	      for(int j=0; j<6; j++) {
		if(!minDR_jet[j]) continue;
		if(jet4 == minDR_jet[j]) { flavGen4 = dpars[j]->pdgId; }
	      }
	    }

	    TLorentzVector vJet [4];
	    
	    vJet[0].SetPtEtaPhiM(jet1->pt, jet1->eta, jet1->phi, jet1->mass);
	    vJet[1].SetPtEtaPhiM(jet2->pt, jet2->eta, jet2->phi, jet2->mass);
	    vJet[2].SetPtEtaPhiM(jet3->pt, jet3->eta, jet3->phi, jet3->mass);
	    vJet[3].SetPtEtaPhiM(jet4->pt, jet4->eta, jet4->phi, jet4->mass);
	    
	    TLorentzVector vDijet [6];
	    vDijet[0] = vJet[0]+vJet[1];
	    vDijet[1] = vJet[0]+vJet[2];
	    vDijet[2] = vJet[0]+vJet[3];
	    vDijet[3] = vJet[1]+vJet[2];
	    vDijet[4] = vJet[1]+vJet[3];
	    vDijet[5] = vJet[2]+vJet[3];
	    
	    for(unsigned int j=0; j<6; j++){
	      if(fabs(vDijet[j].M() - W_MASS) > 20) continue;
	    }
	    
	    TVector2 vPull1; vPull1.Set(jet1->pullY, jet1->pullPhi);
	    TVector2 vPull2; vPull2.Set(jet2->pullY, jet2->pullPhi);
	    TVector2 vPull3; vPull3.Set(jet3->pullY, jet3->pullPhi);
	    TVector2 vPull4; vPull4.Set(jet4->pullY, jet4->pullPhi);
	    
	    //
	    // Fill output tree
	    //
	    runNum    = info->runNum;
	    lumiSec   = info->lumiSec;
	    evtNum    = info->evtNum;
	    metfilter = info->metFilterFailBits;
	    npv       = pvArr->GetEntriesFast();
	    npu       = info->nPUmean;
	    njets     = goodJets.size();
	    rho       = info->rhoJet;
	    scale1fb  = weight;
	    ttWeight  = ttWgt;
	    pfmet     = info->pfMETC;
	    pfmetphi  = info->pfMETCphi;
	    pfmt      = computeMt(vLep.Pt(), vLep.Phi(), info->pfMET, info->pfMETphi);
	    mvamet    = info->mvaMET0;
	    mvametphi = info->mvaMET0phi;
	    mvamt     = computeMt(vLep.Pt(), vLep.Phi(), info->mvaMET0, info->mvaMET0phi);
	    
	    lepId  = lid;
	    lepton = &vLep;
	    
	    jet1vec      = &vJet[0];
	    jet1pull     = &vPull1;
	    jet1flavGen  = flavGen1;
	    jet1qgid     = jet1->qgid;
	    jet1csv      = jet1->csv;
	    jet1pmass    = jet1->prunedm;
	    jet1q        = jet1->q;
	    jet1q03      = jet1->q03;
	    jet1q2       = jet1->q2;
	    jet1ptD      = jet1->ptD;
	    jet1flavAlgo = jet1->mcFlavor;
	    jet1flavPhys = jet1->mcFlavorPhys;
	    jet1npar     = jet1->nParticles;
	    
	    jet2vec      = &vJet[1];
	    jet2pull     = &vPull2;
	    jet2flavGen  = flavGen2;
	    jet2qgid     = jet2->qgid;
	    jet2csv      = jet2->csv;
	    jet2pmass    = jet2->prunedm;
	    jet2q        = jet2->q;
	    jet2q03      = jet2->q03;
	    jet2q2       = jet2->q2;
	    jet2ptD      = jet2->ptD;
 	    jet2flavAlgo = jet2->mcFlavor;
	    jet2flavPhys = jet2->mcFlavorPhys;
	    jet2npar     = jet2->nParticles;
	    
	    jet3vec      = &vJet[2];
	    jet3pull     = &vPull3;
	    jet3flavGen  = flavGen3;
	    jet3qgid     = jet3->qgid;
	    jet3csv      = jet3->csv;
	    jet3pmass    = jet3->prunedm;
	    jet3q        = jet3->q;
	    jet3q03      = jet3->q03;
	    jet3q2       = jet3->q2;
	    jet3ptD      = jet3->ptD;
	    jet3flavAlgo = jet3->mcFlavor;
	    jet3flavPhys = jet3->mcFlavorPhys;
	    jet3npar     = jet3->nParticles;
	    
	    jet4vec      = &vJet[3];
	    jet4pull     = &vPull4;
	    jet4flavGen  = flavGen4;
	    jet4qgid     = jet4->qgid;
	    jet4csv      = jet4->csv;
	    jet4pmass    = jet4->prunedm;
	    jet4q        = jet4->q;
	    jet4q03      = jet4->q03;
	    jet4q2       = jet4->q2;
	    jet4ptD      = jet4->ptD;
	    jet4flavAlgo = jet4->mcFlavor;
	    jet4flavPhys = jet4->mcFlavorPhys;
	    jet4npar     = jet4->nParticles;

	    outTree->Fill();
      	  }
	}
      }
    }
   
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
bool passJetID(const baconhep::TJet *jet)
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
     
  // PU Jet ID
  if     (0    <= fabs(jet->eta) && fabs(jet->eta) < 2.5  && jet->mva < -0.63) return false;
  else if(2.5  <= fabs(jet->eta) && fabs(jet->eta) < 2.75 && jet->mva < -0.60) return false;
  else if(2.75 <= fabs(jet->eta) && fabs(jet->eta) < 3    && jet->mva < -0.55) return false;
  else if(3    <= fabs(jet->eta) && fabs(jet->eta) < 5    && jet->mva < -0.45) return false;
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool passEleID(const baconhep::TElectron *electron, const double rho)
{
  if(!(electron->typeBits & baconhep::kEcalDriven))    return false;
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
bool passMuonID(const baconhep::TMuon *muon)
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
void findDaughters(const TClonesArray* genParArr,
                   std::vector<TLorentzVector> &dvecs,
		   std::vector<const baconhep::TGenParticle*> &dpars)
{
  const int W_PDGID      = 24;  // W+
  const int BOTTOM_PDGID = 5;   // b
  const int TOP_PDGID    = 6;   // t  
  
  assert(genParArr);
  
  std::vector<TLorentzVector> vecs[2];
  std::vector<const baconhep::TGenParticle*> pars[2];

  int tops[2], Ws_3[2], bs_3[2];
  for(int i=0; i<2; i++) { tops[i] = Ws_3[i] = bs_3[i] = -1; }

  // loop through gen particles
  // Note: assumes that parent particles are listed before their daughters
  for(int j=0; j<genParArr->GetEntriesFast(); j++) { 
    const baconhep::TGenParticle *genp = (baconhep::TGenParticle*)genParArr->At(j);

    if(genp->status == 3) {
      // find top quarks
      if(genp->pdgId == -TOP_PDGID) { tops[0] = j; }
      if(genp->pdgId ==  TOP_PDGID) { tops[1] = j; }
      
      // find W daughters from tops
      if(genp->pdgId == -W_PDGID && genp->parent == tops[0]) { Ws_3[0] = j; }
      if(genp->pdgId ==  W_PDGID && genp->parent == tops[1]) { Ws_3[1] = j; }
      
      // find b-quark daughters from tops
      if(genp->pdgId == -BOTTOM_PDGID && genp->parent == tops[0]) { bs_3[0] = j; }
      if(genp->pdgId ==  BOTTOM_PDGID && genp->parent == tops[1]) { bs_3[1] = j; }
      
      // find daughters from Ws
      for(int k=0; k<2; k++) {
    	if(genp->parent == Ws_3[k]) {
    	  TLorentzVector v; v.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
	  vecs[k].push_back(v);
	  pars[k].push_back(genp);
    	}
      }
    }
    if(pars[0].size()==2 && pars[1].size()==2) break;
  }
  
  // Add b's
  const baconhep::TGenParticle *b0 = (baconhep::TGenParticle*)genParArr->At(bs_3[0]);
  TLorentzVector vb0; vb0.SetPtEtaPhiM(b0->pt, b0->eta, b0->phi, b0->mass);
  vecs[0].push_back(vb0);
  pars[0].push_back(b0);
  
  const baconhep::TGenParticle *b1 = (baconhep::TGenParticle*)genParArr->At(bs_3[1]);
  TLorentzVector vb1; vb1.SetPtEtaPhiM(b1->pt, b1->eta, b1->phi, b1->mass);
  vecs[1].push_back(vb1);
  pars[1].push_back(b1);

  if(vecs[0].size()!=3 || pars[0].size()!=3 || vecs[1].size()!=3 || pars[1].size()!=3) return;  // temp fix! should not be needed with latest mother find fix in bacon...

  for(int k=0; k<2; k++) {
    for(unsigned int i=0; i<pars[k].size(); i++) {
      dvecs.push_back(vecs[k][i]);
      dpars.push_back(pars[k][i]);
    }
  }
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
