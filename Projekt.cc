#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//My contribution from 15/12/2022
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>





using namespace std;

typedef struct{
  Int_t event;
  Int_t track;
  Int_t pid;
  Int_t superlayer;
  Int_t station;
  Int_t wheel;
  Int_t sector;
  Int_t layer;
  Double_t localX;
  Double_t localY;
  Double_t localZ;
  Double_t globalX;
  Double_t globalY;
  Double_t globalZ;
  Double_t pt;
  Double_t phi;
  Double_t r;
  Double_t eta;
  Double_t tof;
  Double_t beta;
  Double_t distanceVTX;
  Double_t distanceIP;
  Double_t nextLocalX;
  Double_t nextLocalY;
  Double_t nextLocalZ;
  Double_t nextGlobalX;
  Double_t nextGlobalY;
  Double_t nextGlobalZ;
  Double_t globalHitDistDiff;
  Double_t betaNextHit;
  Double_t diffTOFGeant;
  Double_t diffTOFHitB;
  Double_t diffTOFNextB;
  bool samewheel;
  bool samesl;
  bool samesect;
  bool samest;
  } SimHitData;

typedef struct{
  Int_t pdgId;
  Int_t event;
  Double_t phi;
  Double_t theta;
  Double_t eta;
  Double_t pT;
  Double_t mass;
  Double_t vx;
  Double_t vy;
  Double_t vz;
  Double_t beta;
  Double_t invbeta;
}RecCandidateData;  

typedef struct{
  Int_t pdgId;
  Int_t event;
  Double_t phi;
  Double_t theta;
  Double_t eta;
  Double_t pT;
  Double_t mass;
  Double_t beta;
  Double_t invbeta;
}GenCandidateData; 

typedef struct{
  Double_t pgen;
  Double_t preco;
  Double_t invbetagen;
  Double_t invbetamanual;
  Double_t invbetareco;
  Double_t invbetarecoerror;
  Double_t event;
  Double_t deltaR;
  Double_t genPt;
  Double_t innerrecoPt;
  Double_t outerrecoPt;
  Double_t recoPt;
  Double_t genEta;
  Double_t recoEta;
  Int_t stationMatches;
}InverseBetaData;

typedef struct{
  Double_t pt;
  Double_t eta;
  Double_t p;
  Double_t invbeta;
}TotalMCData;

typedef struct{
  Double_t pt;
  Double_t eta;
  Double_t p;
  Double_t invbeta;
  Double_t deltaR;
}GMTMuonData;

//object definition
class Projekt : public edm::EDAnalyzer {
public:

  //constructor, function is called when new object is created
  explicit Projekt(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Projekt();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  unsigned int hscpCount;
  unsigned int nlines;
  Double_t hscpMass;
  unsigned int genMCCount;
  unsigned int genMCPassCount;
  unsigned int recCount;
  unsigned int recPassCount;
  unsigned int counter;
  unsigned int COUNTGEN;
  //added variables
  TFile *myRootFile;
  TH1F *histoGenDeltaR;
  TH1F *histoGenDeltaPhi;
  TH1F *histoGenDeltaEta;
  TH1F *histoDeltaR;
  TH1F *histoMinDeltaR;

  TH2F *histoMuonStationsTotal;
  TH2F *histoMuonStationsSelected;

  TEfficiency *histoPtEfficiency;
  TEfficiency *histoPEfficiency;
  TEfficiency *histoEtaEfficiency;
  TEfficiency *histoInvbetaEfficiency;
  TEfficiency *histoPtEfficiencyL1T;
  TEfficiency *histoPEfficiencyL1T;
  TEfficiency *histoEtaEfficiencyL1T;
  TEfficiency *histoInvbetaEfficiencyL1T;

  TTree *hscpTree;
  TTree *recoTree;
  TTree *generatedTree;
  TTree *invbetaTree;
  TTree *totalMCTree;
  TTree *gmtMuonTree;

  SimHitData hitData;
  RecCandidateData recData;
  GenCandidateData generatedData;
  InverseBetaData invbetaData;
  TotalMCData totalMCData;
  GMTMuonData gmtMuonData;

  ////////////////////////////////
// Definitions of various inputs // 
  //////////////////////////////
  edm::EDGetTokenT<edm::SimTrackContainer> inputSim;
  edm::EDGetTokenT<edm::SimVertexContainer> inputVtx;
  edm::EDGetTokenT<TrackingParticleCollection> inputTP;
  edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;
  edm::EDGetTokenT<vector<reco::GenParticle> > inputGP;

//My simHits contribution to the code:
  //edm::EDGetTokenT<vector<PSimHit>> inputHitsDT;
  edm::EDGetTokenT<vector<reco::Muon>> inputRecMuons;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonTimeExtra>> inputRecMuonsExtra;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> theGeometryToken;
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> theDTGeomToken;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> theCSCGeomToken;
  const edm::ESGetToken<RPCGeometry, MuonGeometryRecord> theRPCGeomToken;

  edm::EDGetTokenT<BXVector<l1t::Muon>> theGmtToken;  

};

//Definitions of various functions

/*
bool match(const TrackingParticle & tp, const l1t::TrackerMuon & gmt) {
  return (   (fabs(tp.pt()-gmt.trkPtr()->momentum().perp()) )/tp.pt() < 0.1
            && fabs(tp.phi()-gmt.trkPtr()->momentum().phi()) < 0.1
            && fabs(tp.eta()-gmt.trkPtr()->momentum().eta()) < 0.1 );
}*/

string print(const TrackingParticle & tp) {
  stringstream ss;
  ss << tp.pdgId()
     <<" pt: "<<tp.pt()
     <<" eta: "<<tp.eta()
     <<" phi: "<<tp.phi()
     <<" vtx[r,z]:  ["<<tp.parentVertex()->position().Rho() <<", "<<tp.parentVertex()->position().z()<<"]"
     <<" time: "<<tp.parentVertex()->position().T() 
     ; 
  return ss.str();
}

string print(const reco::GenParticle & gp) {
  stringstream ss;
  ss << gp.pdgId()
     <<" pt: "<<gp.pt()
     <<" eta: "<<gp.eta()
     <<" phi: "<<gp.phi()
     <<" eta test: " << gp.momentum().eta();
     ; 
  return ss.str();
}

bool isTightMuon(const reco::Muon& muon){
   if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;
   bool muID = muon::isGoodMuon(muon,muon::GlobalMuonPromptTight) &&
                                  (muon.numberOfMatchedStations() > 1);
   bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
                        muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
   bool ip = true;
   return muID && hits && ip;
}

const TrackingParticle & ancestor(const TrackingParticle & particle) {

  const TrackingVertexRef&  tpv = particle.parentVertex(); 
  if (tpv->nSourceTracks() == 0) return particle;
  const TrackingParticle & parent =  **(tpv->sourceTracks_begin());
  return ancestor(parent);
}

Projekt::Projekt(const edm::ParameterSet& conf) 
  : theConfig(conf), theEventCount(0), hscpCount(0), nlines(0), hscpMass(0), genMCCount(0), genMCPassCount(0), recCount(0), recPassCount(0), counter(0), COUNTGEN(0), theGeometryToken(esConsumes()), theDTGeomToken(esConsumes()), theCSCGeomToken(esConsumes()), theRPCGeomToken(esConsumes())
{
  cout <<" CTORXX" << endl;
//  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(theConfig.getParameter<edm::InputTag>("inputOMTF") );
 // inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
 // inputVtx =  consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));
//  inputGMT =  consumes< vector<l1t::TrackerMuon> >(edm::InputTag("gmtMuons"));
  inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
  inputGP  =  consumes< vector<reco::GenParticle> >(edm::InputTag("genParticles"));

  //nputHitsDT = consumes<vector<PSimHit>>(edm::InputTag("g4SimHits","MuonDTHits"));
  inputRecMuons = consumes<vector<reco::Muon>>(edm::InputTag("muons",""));
  inputRecMuonsExtra = consumes<edm::ValueMap<reco::MuonTimeExtra>>(edm::InputTag("muons","combined"));
  
  theGmtToken = consumes<BXVector<l1t::Muon>>(edm::InputTag("gmtStage2Digis","Muon"));
 //Order of input tags needs to be the same as in edmDumpEventContent

}


Projekt::~Projekt() 
{ 
  cout <<" DTOR" << endl;
}

void Projekt::beginJob()
{

  /////////////////////////////////////////////////////////
// FILES FOR PSEUDORAPIDITY ANALYSIS - files use simTracks //
  /////////////////////////////////////////////////////////

//FEVTSIM_stau_M200_full.root - FEVTSIM.root analysis
//stau_M432_analysis.root - stau_432.root analysis (GEN AND RECO ONLY)

  myRootFile = new TFile("stau_M432_analysis.root","RECREATE"); //remember to change name when changing datafiles!
  histoDeltaR = new TH1F("histoDeltaR","histoDeltaR",300,0,6);
  histoMinDeltaR = new TH1F("histoMinDeltaR","histoMinDeltaR",300,0,6);
  histoGenDeltaR = new TH1F("histoGenDeltaR","histoGenDeltaR",300,0,6);
  histoGenDeltaPhi = new TH1F("histoGenDeltaPhi","histoGenDeltaPhi",100,0,7);
  histoGenDeltaEta = new TH1F("histoGenDeltaEta","histoGenDeltaEta",200,0,6);

  histoMuonStationsTotal = new TH2F("histoMuonStationsTotal","histoMuonStationsTotal", 30, 0, 6, 6 , 0 , 6);
  histoMuonStationsSelected = new TH2F("histoMuonStationsSelected","histoMuonStationsSelected", 30, 0, 6 , 6, 0 , 6);

//Gen refers to total gen candidates, reco refers to gen candidates matched with reco candidates

  const Int_t binsPT = 24;
  Double_t binsPTEdges[binsPT+1] = {0.,20.,40.,60.,80.,100.,150.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.,900.,1000.,1200.,1400.,1600.,1800.,2000.,3000.};

  const Int_t binsP = 26;
  Double_t binsPEdges[binsP+1] = {0.,20.,40.,60.,80.,100.,150.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.,900.,1000.,1200.,1400.,1600.,1800.,2000.,3000.,4000.,5000.};

  const Int_t binsINVB = 19;
  Double_t binsINVBEdges[binsINVB+1] = {1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,5.0,6.0};
  histoPtEfficiency = new TEfficiency("histoPtEfficiency",";p_{T} [GeV];#epsilon",binsPT,binsPTEdges);
  histoPEfficiency = new TEfficiency("histoPEfficiency", ";p [GeV];#epsilon", binsP, binsPEdges);
  histoEtaEfficiency = new TEfficiency("histoEtaEfficiency", ";#eta;#epsilon", 100, -2.8, 2.8);
  histoInvbetaEfficiency = new TEfficiency("histoInvbetaEfficiency", ";1/#beta;#epsilon", binsINVB, binsINVBEdges);

  histoPtEfficiencyL1T = new TEfficiency("histoPtEfficiencyL1T",";p_{T} [GeV];#epsilon", binsPT, binsPTEdges);
  histoPEfficiencyL1T = new TEfficiency("histoPEfficiencyL1T", ";p [GeV];#epsilon", binsP, binsPEdges);
  histoEtaEfficiencyL1T = new TEfficiency("histoEtaEfficiencyL1T", ";#eta;#epsilon", 100, -2.8, 2.8);
  histoInvbetaEfficiencyL1T = new TEfficiency("histoInvbetaEfficiencyL1T", ";1/#beta;#epsilon", binsINVB, binsINVBEdges);


  recoTree = new TTree("recoTree", "recoTree");
  recoTree->Branch("PDGID", &(recData.pdgId),32000, 99);
  recoTree->Branch("EVENT", &(recData.event),32000, 99);
  recoTree->Branch("PHI", &(recData.phi),32000, 99);
  recoTree->Branch("THETA", &(recData.theta),32000, 99);
  recoTree->Branch("PT", &(recData.pT),32000, 99);
  recoTree->Branch("ETA", &(recData.eta),32000, 99);
  recoTree->Branch("MASS", &(recData.mass),32000, 99);
  recoTree->Branch("VX", &(recData.vx),32000, 99);
  recoTree->Branch("VY", &(recData.vy),32000, 99);
  recoTree->Branch("VZ", &(recData.vz),32000, 99);
  recoTree->Branch("BETA",&(recData.beta),32000,99);
  recoTree->Branch("INVBETA",&(recData.invbeta),32000,99);

  generatedTree = new TTree("generatedTree", "generatedTree");
  generatedTree->Branch("PDGID", &(generatedData.pdgId),32000, 99);
  generatedTree->Branch("EVENT", &(generatedData.event),32000, 99);
  generatedTree->Branch("PHI", &(generatedData.phi),32000, 99);
  generatedTree->Branch("THETA", &(generatedData.theta),32000, 99);
  generatedTree->Branch("PT", &(generatedData.pT),32000, 99);
  generatedTree->Branch("ETA", &(generatedData.eta),32000, 99);
  generatedTree->Branch("MASS", &(generatedData.mass),32000, 99);
  generatedTree->Branch("BETA", &(generatedData.beta),32000, 99);
  generatedTree->Branch("INVBETA", &(generatedData.invbeta),32000, 99);

  hscpTree = new TTree("simTree", "simTree");
  hscpTree->Branch("EV", &(hitData.event), 32000, 99);
  hscpTree->Branch("TR", &(hitData.track), 32000, 99);
  hscpTree->Branch("PID", &(hitData.pid), 32000, 99);
  hscpTree->Branch("SL", &(hitData.superlayer), 32000, 99);
  hscpTree->Branch("ST", &(hitData.station), 32000, 99);
  hscpTree->Branch("WH", &(hitData.wheel), 32000, 99);
  hscpTree->Branch("SECT", &(hitData.sector), 32000, 99);
  hscpTree->Branch("LYR", &(hitData.layer), 32000, 99);
  hscpTree->Branch("X_L", &(hitData.localX), 32000, 99);
  hscpTree->Branch("Y_L", &(hitData.localY), 32000, 99);
  hscpTree->Branch("Z_L", &(hitData.localZ), 32000, 99);
  hscpTree->Branch("X_G", &(hitData.globalX), 32000, 99);
  hscpTree->Branch("Y_G", &(hitData.globalY), 32000, 99);
  hscpTree->Branch("Z_G", &(hitData.globalZ), 32000, 99);
  hscpTree->Branch("PT", &(hitData.pt), 32000, 99);
  hscpTree->Branch("PHI", &(hitData.phi), 32000, 99);
  hscpTree->Branch("R", &(hitData.r), 32000, 99);
  hscpTree->Branch("ETA", &(hitData.eta), 32000, 99);
  hscpTree->Branch("TOF", &(hitData.tof), 32000, 99);
  hscpTree->Branch("BETA", &(hitData.beta), 32000, 99);
  hscpTree->Branch("DISTVTX", &(hitData.distanceVTX), 32000, 99);
  hscpTree->Branch("DISTIP", &(hitData.distanceIP), 32000, 99);
  hscpTree->Branch("NX_L", &(hitData.nextLocalX), 32000, 99);
  hscpTree->Branch("NY_L", &(hitData.nextLocalY), 32000, 99);
  hscpTree->Branch("NZ_L", &(hitData.nextLocalZ), 32000, 99);
  hscpTree->Branch("NX_G", &(hitData.nextGlobalX), 32000, 99);
  hscpTree->Branch("NY_G", &(hitData.nextGlobalY), 32000, 99);
  hscpTree->Branch("NZ_G", &(hitData.nextGlobalZ), 32000, 99);
  hscpTree->Branch("DISTDIFF_G", &(hitData.globalHitDistDiff), 32000, 99);
  hscpTree->Branch("NBETA", &(hitData.betaNextHit), 32000, 99);
  hscpTree->Branch("GEANTTOFDIFF", &(hitData.diffTOFGeant), 32000, 99);
  hscpTree->Branch("HITPTOFDIFF", &(hitData.diffTOFHitB), 32000, 99);
  hscpTree->Branch("NEXTPTOFFDIFF", &(hitData.diffTOFNextB), 32000, 99);
  hscpTree->Branch("SAMEWHEEL", &(hitData.samewheel), 32000, 99);
  hscpTree->Branch("SAMESL", &(hitData.samesl), 32000, 99);
  hscpTree->Branch("SAMESECT", &(hitData.samesect), 32000, 99);
  hscpTree->Branch("SAMEST", &(hitData.samest), 32000, 99);

  invbetaTree = new TTree("invbetaTree", "invbetaTree");
  invbetaTree->Branch("PGEN",&(invbetaData.pgen), 32000, 99);
  invbetaTree->Branch("PRECO",&(invbetaData.preco), 32000, 99);
  invbetaTree->Branch("INVBETAGEN", &(invbetaData.invbetagen), 32000, 99);
  invbetaTree->Branch("INVBETAMANUAL", &(invbetaData.invbetamanual), 32000, 99);
  invbetaTree->Branch("INVBETARECO", &(invbetaData.invbetareco), 32000,99);
  invbetaTree->Branch("INVBETARECOERROR",&(invbetaData.invbetarecoerror), 32000,99);
  invbetaTree->Branch("EVENT", &(invbetaData.event), 32000, 99);
  invbetaTree->Branch("DELTAR", &(invbetaData.deltaR), 32000, 99);
  invbetaTree->Branch("GENPT",&(invbetaData.genPt), 32000, 99);
  invbetaTree->Branch("INNERRECOPT", &(invbetaData.innerrecoPt), 32000,99);
  invbetaTree->Branch("OUTERRECOPT", &(invbetaData.outerrecoPt), 32000,99);
  invbetaTree->Branch("RECOPT",&(invbetaData.recoPt), 32000,99);
  invbetaTree->Branch("RECOETA",&(invbetaData.recoEta), 32000,99);
  invbetaTree->Branch("GENETA",&(invbetaData.genEta), 32000,99);
  invbetaTree->Branch("STMATCH",&(invbetaData.stationMatches), 32000,99);


  totalMCTree = new TTree("totalMCTree", "totalMCTree");
  totalMCTree->Branch("PT", &(totalMCData.pt), 32000, 99);
  totalMCTree->Branch("P", &(totalMCData.p), 32000, 99);
  totalMCTree->Branch("ETA", &(totalMCData.eta), 32000, 99);
  totalMCTree->Branch("INVBETA", &(totalMCData.invbeta), 32000, 99);

  gmtMuonTree = new TTree("gmtMuonTree", "gmtMuonTree");
  gmtMuonTree->Branch("PT", &(gmtMuonData.pt), 32000, 99);
  gmtMuonTree->Branch("P", &(gmtMuonData.p), 32000, 99);
  gmtMuonTree->Branch("ETA", &(gmtMuonData.eta), 32000, 99);
  gmtMuonTree->Branch("INVBETA", &(gmtMuonData.invbeta), 32000, 99);
  gmtMuonTree->Branch("DELTAR",&(gmtMuonData.deltaR), 32000, 99);

  ////////////////////////////  //////////////////////////////
// PSEUDORAPIDITY HISTOGRAMS SPECIFYING MUON SYSTEM REGIONS //
  //////////////////////////////////////////////////////////
  cout << "HERE Projekt::beginJob()" << endl;
}

void Projekt::endJob()
{
  cout << "MC HSCP: " << genMCCount << "PASS: " << genMCPassCount << endl;
  cout << "REC MUON: " << recCount << "PASS: " << recPassCount << endl;
  cout << "TEfficiency entries: " << counter << endl;
  cout << "Stau gen in detector: " << COUNTGEN << endl;
  cout << "Total nr of events: " << theEventCount << endl;

  recoTree->Write();
  hscpTree->Write();
  generatedTree->Write();
  invbetaTree->Write();
  totalMCTree->Write();
  gmtMuonTree->Write();

  histoDeltaR->Write();
  histoMinDeltaR->Write();
  histoGenDeltaR->Write();
  histoGenDeltaPhi->Write();
  histoGenDeltaEta->Write();

  histoPtEfficiency->Write();
  histoPEfficiency->Write();
  histoEtaEfficiency->Write();
  histoInvbetaEfficiency->Write();

  histoPtEfficiencyL1T->Write();
  histoPEfficiencyL1T->Write();
  histoEtaEfficiencyL1T->Write();
  histoInvbetaEfficiencyL1T->Write();

  histoMuonStationsTotal->Write();
  histoMuonStationsSelected->Write();         

  myRootFile->Close();

  delete histoDeltaR;
  delete histoMinDeltaR;
  delete histoGenDeltaR;
  delete histoGenDeltaPhi;
  delete histoGenDeltaEta;

  delete histoPtEfficiency;
  delete histoPEfficiency;
  delete histoEtaEfficiency;
  delete histoInvbetaEfficiency;

  delete histoPtEfficiencyL1T;
  delete histoPEfficiencyL1T;
  delete histoEtaEfficiencyL1T;
  delete histoInvbetaEfficiencyL1T;

  delete histoMuonStationsTotal;
  delete histoMuonStationsSelected;

  delete myRootFile;
  cout << "HERE Cwiczenie::endJob()" << endl;
}

void Projekt::analyze(
    const edm::Event& ev, const edm::EventSetup& es){
  cout << " -------------------------------- HERE Cwiczenie::analyze "<< endl;
  cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" useful event count:"<<++theEventCount << endl;

//  edm::Handle<vector<l1t::TrackerMuon> > gmtColl;
//  ev.getByToken(inputGMT, gmtColl);
//  const vector<l1t::TrackerMuon> & gmtMuons = *gmtColl.product();
//  histo->Fill(gmtMuons.size());
    
  ////////////////////////////////////////
// Assigning global geometry to analysis //
  //////////////////////////////////////

  auto const & globalGeometry = es.getData(theGeometryToken);  

  ////////////////////////////////////////////////////////////////////////////////////////////////////
// GENERATED PARTICLE ANALYSIS (MC LEVEL)                                                             //
// Definition of genParticle vector, which is needed at this stage to extract generated particle mass //
// (Value of mass for all HSCPs should be the same)                                                   //
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const vector<reco::GenParticle> & genParticles = ev.get(inputGP);
  Int_t genStauCandidates = 0;

  for (const auto & gp : genParticles) {
  	if (abs(gp.pdgId()) == 13 && gp.status() == 1){ //generated particle is a muon (which must be stable but that's a given)
      double muon_p = gp.pt()*cosh(gp.eta());
			double muon_phi = gp.phi();
      if (muon_phi <= 0){
				muon_phi += 2*M_PI;
			}
			double muon_beta = muon_p/gp.energy();
    }
    else if (abs(gp.pdgId())>1000000) { //generated particle is BSM particle (e.g. stau->pdgID=1000015
			if(gp.status()==1){ // particle is stable
        genMCCount++;
        hscpMass = gp.mass();
				double hscp_p = gp.pt()*cosh(gp.eta());
				double hscp_phi = gp.phi();
        double hscp_theta = gp.theta();
				if(hscp_phi <= 0){
					hscp_phi += 2*M_PI;
				}
				double hscp_beta = hscp_p/gp.energy();
      
      if(gp.pt()>50 && abs(gp.eta())<2.1 ){
      genMCPassCount++;  
      genStauCandidates++;
      generatedData.pdgId = gp.pdgId();
      generatedData.event = theEventCount;
      generatedData.phi = hscp_phi;
      generatedData.beta = hscp_beta;
      generatedData.invbeta = 1/hscp_beta;
      generatedData.theta = hscp_theta;
      generatedData.pT = gp.pt();
      generatedData.eta = gp.eta();
      generatedData.mass = hscpMass;
      generatedTree->Fill();
      }
    }
  }
  if(gp.status()==1){
    totalMCData.pt = gp.pt();
    totalMCData.p = gp.p();
    totalMCData.eta = gp.eta();
    totalMCData.invbeta = gp.energy()/gp.p();
    totalMCTree->Fill();
  }
}
   

  ///////////////////////////////////
// SimTrack/SimVtx vector definitions //
  //////////////////////////////////

/*  edm::Handle<edm::SimTrackContainer> simTrk;
  ev.getByToken(inputSim, simTrk);
  const vector<SimTrack>  & mySimTracks = *(simTrk.product());
  cout <<" SIMULATED TRACKS: "<<mySimTracks.size()<<endl;

  edm::Handle<edm::SimVertexContainer> simVtx;
  ev.getByToken(inputVtx, simVtx);
  const vector<SimVertex> & mySimVerts= *(simVtx.product());
  cout <<" SIMULATED VERTICES: "<<mySimVerts.size()<<endl;

  //////////////////
// SIMHIT ANALYSIS //
  ////////////////

 const vector<PSimHit> & simDTHits = ev.get(inputHitsDT);
  cout <<"Number of simulated DT hits in event: "<<simDTHits.size() << endl;

  for(const auto &track:mySimTracks){
    if(track.trackId()>10) break;
    cout << "Track: " <<track<< endl;
    cout << track.trackId() << ", " << track.vertIndex() << ", " << endl;
    cout << "Position rho = " << mySimVerts[track.vertIndex()].position().Rho() << ", Position z = " << mySimVerts[track.vertIndex()].position().z() << endl;
    cout << "Parent index: " <<  mySimVerts[track.vertIndex()].parentIndex() << endl; 
  }
  //////////////
// DT CHAMBERS //
  ////////////

  Int_t hitCount = 0;

  for (vector<PSimHit>::const_iterator iter=simDTHits.begin();iter<simDTHits.end();iter++){ //Iterator is from 0
    
    const PSimHit & hit = *iter;

    if(abs(hit.particleType())>1000000){
      DTLayerId dtDetLayerId(hit.detUnitId()); 
      cout << "======================================HIT SUMMARY==========================================" << std::endl;
      
      cout << dtDetLayerId << endl;
      Int_t eventNr = theEventCount;
	    Int_t pid = hit.particleType();
      Int_t trackNr = hit.trackId(); 
      Local3DPoint localPosition = hit.localPosition();
      GlobalPoint globalPosition = globalGeometry.idToDet(dtDetLayerId)->toGlobal(hit.localPosition());
      Double_t r = sqrt((globalPosition.x()*globalPosition.x())+(globalPosition.y()*globalPosition.y()));
      LocalVector p = hit.momentumAtEntry();
      Double_t pt = sqrt((p.y()*p.y())+(p.z()*p.z()));
      Double_t phi = globalPosition.phi();
      if(phi<0){phi = phi + 2*M_PI;}
      Double_t theta = globalPosition.theta();
      Double_t eta = -log(tan(abs(theta)/2));
      Double_t tof = hit.timeOfFlight();
      Double_t distanceIP = globalPosition.mag();

     math::XYZTLorentzVectorD hscpVertexV = mySimVerts[mySimTracks[hit.trackId()-1].vertIndex()].position();
      cout << "HIT VERTEX (x0,y0,z0) = (" << hscpVertexV << ")" << endl;
      Double_t distanceVTX = abs(distanceIP) - abs(hscpVertexV.mag());
      cout << "simHit: " << hit << endl;
      cout << "Track ID: " << hit.trackId() << " | Det Unit ID: " << hit.detUnitId() << " | PID: " << hit.particleType() << " | p: "<< hit.momentumAtEntry() <<" | phi: " << hit.phiAtEntry() << " | theta: " << hit.thetaAtEntry() << " | TOF: " << hit.timeOfFlight() << endl;
             
  ////////////////////////////
// TOF analysis (preliminary) // 
  ////////////////////////////

	   Double_t betaHit = p.mag()/sqrt((p.mag()*p.mag())+(hscpMass*hscpMass)); //.pabs() calculates the length of the vector pointing to the hit from (0,0,0) 
    															//mass and momenta in GeV (no unit conversion required)

	   cout << "#Beta(hit) = " << betaHit << endl;

    Double_t tofCalculated = distanceIP*1e9*0.01/(betaHit*TMath::C()); //multiplied by 1e9 and 0.01 for tof to be in ns! 
    cout << "HIT STATION: " << dtDetLayerId.station() << " || HIT SUPERLAYER: " << dtDetLayerId.superlayer() << "|| HIT LAYER: " << dtDetLayerId.layer() << endl;
    cout << "TOF FROM HIT = " << tof << " ns" << " || TOF CALCULATED = " << tofCalculated << " ns" << " || BETA CALCULATED = " << betaHit <<endl;
    cout << "LOCAL POSITION (x,y,z)cm = (" << localPosition.x() << ", " << localPosition.y() << ", " << localPosition.z() << ")" << endl;
    cout << "GLOBAL POSITION (x,y,z)cm = (" << globalPosition.x() << ", " << globalPosition.y() << ", " << globalPosition.z() << ")" << endl;
    cout << "HIT DISTANCE FROM IP = " << distanceIP << " cm" << endl;
    //cout << "HIT DISTANCE FROM VTX = " << distanceVTX << " cm" << endl;
    cout << "HIT MOMENTUM = " << p.mag() << " GeV" << endl;   

  //////////////////
// TOF BETWEEN HITS // 
  //////////////////

      if(iter!=simDTHits.end()){
        const PSimHit & nextHit = *(iter+1);
        if(hit.trackId()==nextHit.trackId()){
        DTLayerId nextDTDetLayerId(nextHit.detUnitId());
        
        Local3DPoint nextHitLocal = nextHit.localPosition();
        GlobalPoint nextHitGlobal = globalGeometry.idToDet(nextDTDetLayerId)->toGlobal(nextHitLocal);  

        Vector3DBase lHitDistV = nextHitLocal - localPosition;
        Vector3DBase gHitDistV = nextHitGlobal - globalPosition;
        cout << "GLOBAL INTERHIT DISTANCE (VECTOR) = " << gHitDistV << endl;

        Double_t gHitDist = sqrt(gHitDistV.x()*gHitDistV.x() + gHitDistV.y()*gHitDistV.y() + gHitDistV.z()*gHitDistV.z());
        Double_t lHitDist = sqrt(lHitDistV.x()*lHitDistV.x() + lHitDistV.y()*lHitDistV.y() + lHitDistV.z()*lHitDistV.z());

        cout << "DISTANCE FROM NEXT HIT (GLOBAL): " << gHitDist << " cm" << endl; 
        cout << "DISTANCE FROM NEXT HIT (LOCAL): " << lHitDist << " cm" << endl;
      
        Double_t betaNextHit = nextHit.pabs()/sqrt((nextHit.pabs()*nextHit.pabs())+(hscpMass*hscpMass));

        Double_t tofGivenDifference = nextHit.timeOfFlight() - tof;
        Double_t tofBetweenHitsCurrent = gHitDist*1e9*0.01/((betaHit)*TMath::C());
        Double_t tofBetweenHitsNext = gHitDist*1e9*0.01/((betaNextHit)*TMath::C());

        cout << "#DeltaTOF FROM HIT TOFs = " << tofGivenDifference << " ns || #DeltaTOF CALCULATED (hit #beta) = " << tofBetweenHitsCurrent << " ns || #DeltaTOF CALCULATED (next hit #beta) = " << tofBetweenHitsNext << endl;
        cout << "PROCESS TYPE: " << hit.processType() << endl;
        cout << "BETA AT HIT: " << hit.pabs()/sqrt((hit.pabs()*hit.pabs())+(hscpMass*hscpMass)) << endl;

        nlines++;
      
      hitData.event = eventNr;
      hitData.track = trackNr;
      hitData.pid = pid;
      hitData.superlayer = dtDetLayerId.superlayer();
      hitData.station = dtDetLayerId.station();
      hitData.wheel = dtDetLayerId.wheel();
      hitData.sector = dtDetLayerId.sector();
      hitData.layer = dtDetLayerId.layer();
      hitData.localX = localPosition.x();
      hitData.localY = localPosition.y();
      hitData.localZ = localPosition.z();
      hitData.globalX = globalPosition.x();
      hitData.globalY = globalPosition.y();
      hitData.globalZ = globalPosition.z();
      hitData.pt = pt;
      hitData.phi = phi;
      hitData.r = r;
      hitData.eta = eta;
      hitData.tof = tof;
      hitData.beta = betaHit;
      //hitData.distanceVTX = distanceVTX;
      hitData.distanceIP = distanceIP;
      hitData.nextLocalX = nextHitLocal.x();
      hitData.nextLocalY = nextHitLocal.y();
      hitData.nextLocalZ = nextHitLocal.z();
      hitData.nextGlobalX = nextHitGlobal.x();
      hitData.nextGlobalY = nextHitGlobal.y();
      hitData.nextGlobalZ = nextHitGlobal.z();
      hitData.globalHitDistDiff = gHitDist;
      hitData.betaNextHit = betaNextHit;
      hitData.diffTOFGeant = tofGivenDifference;
      hitData.diffTOFHitB = tofBetweenHitsCurrent;
      hitData.diffTOFNextB = tofBetweenHitsNext;

      if(nextDTDetLayerId.wheel() == dtDetLayerId.wheel()){hitData.samewheel = true;}
      else{hitData.samewheel = false;}
      if(nextDTDetLayerId.superlayer() == dtDetLayerId.superlayer()){hitData.samesl = true;}
      else{hitData.samesl = false;}
      if(nextDTDetLayerId.sector() == dtDetLayerId.sector()){hitData.samesect = true;}
      else{hitData.samesect = false;}
      if(nextDTDetLayerId.station() == dtDetLayerId.station()){hitData.samest = true;}
      else{hitData.samest = false;}

      hscpTree->Fill();


        }
      hitCount++; //has to be here in order for hit count to match numbering in hit ntuples :)
}
    }
  }*/
  //////////////////////////
// Simulated track analysis //
  //////////////////////////

  /*for (std::vector<SimTrack>::const_iterator it=mySimTracks.begin(); it<mySimTracks.end(); it++) {
    const SimTrack & track = *it;
    if ( track.type() == -99) continue;
    if ( track.vertIndex() != 0) continue;

    double phi_sim = track.momentum().phi(); //momentum azimutal angle
    double pt_sim = track.momentum().pt(); //transverse momentum
    double eta_sim = track.momentum().eta(); //pseudorapidity

    if(track.type()>1000000) hscpCount++;

  }*/

  //////////////////
// RECHIT ANALYSIS //
  ////////////////

  edm::Handle<edm::ValueMap<reco::MuonTimeExtra>> muonExtraHandle;
  ev.getByToken(inputRecMuonsExtra,muonExtraHandle);

  edm::Handle<vector<reco::Muon>> muonHandle;
  ev.getByToken(inputRecMuons, muonHandle);
  const vector<reco::Muon> & recoMuons = *muonHandle.product();
  const unsigned int eventRecMuons = recoMuons.size();
  cout << eventRecMuons << " muons in event" << endl;

  edm::ValueMap<reco::MuonTimeExtra> muonValueMap = *muonExtraHandle;
  

  cout << "=============RECO SUMMARY=============" << endl;

  for(vector<reco::Muon>::const_iterator iter = recoMuons.begin();iter<recoMuons.end();iter++){
  const reco::Muon & recoCandidate = *iter;
  Double_t phiReco = recoCandidate.phi();
  			if(phiReco <= 0){
					phiReco += 2*M_PI;
				}
  Double_t thetaReco = recoCandidate.theta();
  Double_t pTReco = recoCandidate.pt();
  Double_t etaReco = recoCandidate.eta();
  Int_t pdgIdReco = recoCandidate.pdgId();
  Double_t betaReco = recoCandidate.p()/sqrt(recoCandidate.p()*recoCandidate.p()+hscpMass*hscpMass);

  //cout << "No. chambers: " << recoCandidate.numberOfChambers() << endl;    

  //cout << "CANDIDATE " << recoMuonCount << " | PDG ID:" << pdgIdReco << " pT: " << pTReco <<  " phi: " << phiReco << " #eta: " << etaReco << "theta: " << thetaReco << endl;
  //cout << "            Vertex position (x,y,z): " <<  recoCandidate.vertex() << endl;
  recCount++;
if(pTReco>50 && abs(etaReco)<2.1 && isTightMuon(recoCandidate)){
  recData.pdgId = pdgIdReco;
  recData.event = theEventCount;
  recData.phi = phiReco;
  recData.theta = thetaReco;
  recData.pT = pTReco;
  recData.eta = etaReco;
  recData.beta = betaReco;
  recData.invbeta = 1/betaReco;
  recData.mass = recoCandidate.mass();
  recData.vx = recoCandidate.vx();
  recData.vy = recoCandidate.vy();
  recData.vz = recoCandidate.vz();
  recoTree->Fill();
  recPassCount++;
  }
  }

  ////////////////////////////
//  DELTA R GEN-RECO MATCHING // 
  ////////////////////////// 
    
  Int_t countingMatches = 0;
  for (const auto & gp : genParticles) {
      long unsigned int gpMatchCount = 0;
      Double_t genInvbeta = 1/(gp.p()/sqrt(gp.p()*gp.p() + hscpMass*hscpMass));

    if(gp.pt()>50 && abs(gp.eta())<2.1 && abs(gp.pdgId())>1000000 && gp.status()==1){
    
      for (const auto & gp2 : genParticles){
        if(gp2.pt()>50 && abs(gp2.eta())<2.1 && abs(gp2.pdgId())>1000000 && gp2.status()==1){
        if(gpMatchCount==1){
          continue;
        }
        Double_t deltaR = reco::deltaR(gp,gp2);
        Double_t deltaPhi = abs(gp.phi() - gp2.phi());
        Double_t deltaEta = abs(gp.eta() - gp2.eta());
        cout << "debug2: " << deltaR << endl;
        
        if(deltaR>0.1){
            cout << "ETA 1: " << gp.eta() << " PHI: " << gp.phi() << " PT: " << gp.pt() << " DELTA R: " << deltaR << endl;
            cout << "ETA 2: " << gp2.eta() << " PHI: " << gp2.phi() << " PT: " << gp2.pt() << " DELTA R: " << deltaR << endl;
            cout << print(gp) << endl;
            cout << print(gp2) << endl << endl;
        
            histoGenDeltaR->Fill(deltaR);
            histoGenDeltaPhi->Fill(deltaEta);
            histoGenDeltaEta->Fill(deltaPhi);
            gpMatchCount++;
            
            }
          }
      }
      
      Double_t minDeltaR = 100;
      Int_t recoMuonCount = 0;
      for(vector<reco::Muon>::const_iterator iter = recoMuons.begin();
      iter<recoMuons.end();iter++){
          reco::MuonRef muref(muonHandle,recoMuonCount);
          reco::MuonTimeExtra muonExtraData = (muonValueMap)[muref];
        const reco::Muon & recoCandidate = *iter;
        Double_t pTReco = recoCandidate.pt();
        Double_t etaReco = recoCandidate.eta();

        if(pTReco>50 && abs(etaReco)<2.1 && isTightMuon(recoCandidate) && abs(recoCandidate.pdgId())==13){
          cout << "NEW REC CANDIDATE" << endl;
          Double_t deltaR = reco::deltaR(gp,recoCandidate);
          cout << "here" << endl;
          cout << "delta R: " << deltaR << endl; 
          histoDeltaR->Fill(deltaR);
          if(deltaR<minDeltaR){
            minDeltaR = deltaR;
          }
          
          Double_t innerTrackPt = recoCandidate.innerTrack()->pt();
          Double_t outerTrackPt = recoCandidate.outerTrack()->pt();
          Double_t recoPt = recoCandidate.pt();
          Double_t genPt = gp.pt();

          invbetaData.innerrecoPt = innerTrackPt;
          invbetaData.outerrecoPt = outerTrackPt;
          invbetaData.genPt = genPt;
          invbetaData.recoPt = recoPt;
          Double_t genEta = gp.eta();
          Double_t recoEta = recoCandidate.eta();

          Double_t genP = gp.p();
          Double_t recoP = recoCandidate.p();
          Double_t recInvBetaMan = 1/(recoP/sqrt(recoP*recoP + hscpMass*hscpMass));
          cout << "recoP: " << recoP << "genP: " << genP << endl;
          Double_t recInvBetaRec = muonExtraData.inverseBeta();
          Double_t recInvBetaRecError = muonExtraData.inverseBetaErr();
          cout << "Reco p: " << recoCandidate.p() << "genInvBeta: " << genInvbeta << "recInvBetaMan: " << recInvBetaMan << "recInvBetaRec: " << recInvBetaRec << endl;
          invbetaData.event = theEventCount;
          invbetaData.pgen = genP;
          invbetaData.preco = recoP;
          invbetaData.invbetagen = genInvbeta;
          invbetaData.invbetamanual = recInvBetaMan;
          invbetaData.invbetareco = recInvBetaRec;
          invbetaData.invbetarecoerror = recInvBetaRecError;
          invbetaData.deltaR = deltaR;
          invbetaData.recoEta = recoEta;
          invbetaData.stationMatches = recoCandidate.numberOfMatchedStations();
          Int_t matchedStations = recoCandidate.numberOfMatchedStations();
          if(deltaR<0.1){
            countingMatches++;
          }
          invbetaData.genEta = genEta;
          invbetaTree->Fill();
                    
          }

          recoMuonCount++;
      }
      histoMinDeltaR->Fill(minDeltaR);
      cout << countingMatches << endl;
    }
  }


  //////////////////
// EFFICIENCY PLOTS //
  //////////////////
  for (const auto & gp : genParticles) {
    Double_t invbeta = 1/(gp.p()/sqrt(gp.p()*gp.p() + hscpMass*hscpMass));
    if(abs(gp.pdgId())>1000000 && gp.status()==1 && abs(gp.eta())<2.1){
      cout << "NEW GEN" << endl;
      COUNTGEN++;
      histoMuonStationsTotal->Fill(invbeta,0);
      //histoMuonStationsTotal->Fill(invbeta,1);
      histoMuonStationsTotal->Fill(invbeta,2);
      histoMuonStationsTotal->Fill(invbeta,3);
      histoMuonStationsTotal->Fill(invbeta,4);
      histoMuonStationsTotal->Fill(invbeta,5);
      Int_t rec = 0;
      Int_t matchedStations = 0;
    for(vector<reco::Muon>::const_iterator iter = recoMuons.begin();
      iter<recoMuons.end();iter++){
          const reco::Muon & recoCandidate = *iter;  
          Double_t deltaR = reco::deltaR(gp,recoCandidate);

          if(recoCandidate.pt()>50 && abs(recoCandidate.eta())<2.1 && isTightMuon(recoCandidate) && abs(recoCandidate.pdgId())==13 && deltaR<0.1){
              rec++;
              matchedStations = recoCandidate.numberOfMatchedStations();
            }
      }
      cout << "rec: " << rec << endl;
      if(rec==1){
        histoPtEfficiency->Fill(1,gp.pt());
        histoPEfficiency->Fill(1,gp.p());
        histoEtaEfficiency->Fill(1,gp.eta());
        histoInvbetaEfficiency->Fill(1,invbeta);
        histoMuonStationsSelected->Fill(invbeta,matchedStations);
        counter++;
      }
      else{
        histoPtEfficiency->Fill(0,gp.pt());
        histoPEfficiency->Fill(0,gp.p());
        histoEtaEfficiency->Fill(0,gp.eta());
        histoInvbetaEfficiency->Fill(0,invbeta);
        histoMuonStationsSelected->Fill(invbeta,0);
        counter++;
      }
    }
  }

  //////////////////////////////
// TRIGGER EFFICIENCY ANALYSIS //
  ////////////////////////////

  const l1t::MuonBxCollection & gmts = ev.get(theGmtToken); 
  int bxNumber = 0;

  edm::Handle<BXVector<l1t::Muon>> l1tMuonHandle;
  ev.getByToken(theGmtToken, l1tMuonHandle);
  const BXVector<l1t::Muon> & l1tMuons = *l1tMuonHandle.product();
  const unsigned int eventL1tMuons = l1tMuons.size();
  cout << eventL1tMuons << " L1T muons in event" << endl;

  for (const auto & gp : genParticles) {
    Double_t invbeta = 1/(gp.p()/sqrt(gp.p()*gp.p() + hscpMass*hscpMass));
    if(abs(gp.pdgId())>1000000 && gp.status()==1){
      cout << "NEW GEN" << endl;
      int l1tPt = 0;
      int l1tP = 0;
      int l1tEta = 0;
      int l1tInvbeta = 0;
        for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {

          const l1t::Muon & l1tMuon = *it;  
              gmtMuonData.pt = l1tMuon.pt();
              gmtMuonData.p = l1tMuon.p();
              gmtMuonData.eta = l1tMuon.eta();
              Double_t l1tMuonInvbeta = 1/(l1tMuon.p()/sqrt(l1tMuon.p()*l1tMuon.p() + hscpMass*hscpMass));
              gmtMuonData.invbeta = l1tMuonInvbeta;
              Double_t deltaR = reco::deltaR(gp,l1tMuon);
              gmtMuonData.deltaR = deltaR;
              gmtMuonTree->Fill();

          cout << "pt: " << l1tMuon.pt() << ", eta: " << l1tMuon.eta() << ", deltaR: " << deltaR << endl;
          
          if(abs(gp.eta())<1.6 && deltaR<0.4){ //22 GeV pT cut corresponds to cut in SingleMu trigger (i.e. most basic L1T trigger)
/*pT*/      l1tPt++;
            cout << "pT pass" << endl;
            }
          if(abs(gp.eta())<1.6 && deltaR<0.4){
/*p*/       l1tP++;
            cout << "p pass" << endl;
          } 
          if(deltaR<0.4){
/*eta*/     l1tEta++;
            cout << "eta pass" << endl;
          } 
          if(abs(gp.eta())<1.6 && deltaR<0.4){
/*invbeta*/ l1tInvbeta++;
            cout << "invbeta pass" << endl;
          }
      }
      if(l1tPt>=1){
        histoPtEfficiencyL1T->Fill(1,gp.pt());
      }
      else if(l1tPt==0){
        histoPtEfficiencyL1T->Fill(0,gp.pt());
      }
      if(l1tP>=1){
        histoPEfficiencyL1T->Fill(1,gp.p());
      }
      else if(l1tP==0){
        histoPEfficiencyL1T->Fill(0,gp.p());
      }
      if(l1tEta>=1){
        histoEtaEfficiencyL1T->Fill(1,gp.eta());
      }
      else if(l1tEta==0){
        histoEtaEfficiencyL1T->Fill(0,gp.eta());
      }
      if(l1tInvbeta>=1){
        histoInvbetaEfficiencyL1T->Fill(1,invbeta);
      }
      else if(l1tInvbeta==0){
        histoInvbetaEfficiencyL1T->Fill(0,invbeta);
      }
      }
    }
  } 


DEFINE_FWK_MODULE(Projekt);

