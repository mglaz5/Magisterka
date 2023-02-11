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
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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
  } HitData;

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
  //added variables
  TFile *myRootFile;
  TTree *hscpTree;
  HitData hitData;


  //////////////////////////////
// HSCP PSEUDORAPIDITY ANALYSIS // Here I have histograms which are sometimes only relevant for R hadrons, so I need to remember to comment them out!
  /////////////////////////////

  TH1D *histo_pseudorapidity; //in case of R hadrons, this is for *charged* R hadrons
  TH1D *histo_pseudorapidity_neutral; //only relevant for R hadrons - left in for staus (will just be empty - saves time changing code between datafiles)


  ////////////////////////////////////////  
// KINEMATIC PARAMETERS FOR HSCP AND MUON //
  ////////////////////////////////////////

  TH1D *histo_pdgCount;
  TH1D *histo_pdgCount_HSCP;
  
  TH1D *histo_stau_pt;
  TH1D *histo_stau_eta;
  TH1D *histo_stau_pl;
  TH1D *histo_stau_p;
  TH1D *histo_stau_phi;
  TH1D *histo_stau_beta;
  TH1D *histo_stau_invbeta;

  TH1D *histo_muon_pt;
  TH1D *histo_muon_eta;
  TH1D *histo_muon_pl;
  TH1D *histo_muon_p;
  TH1D *histo_muon_phi;
  TH1D *histo_muon_beta;
  TH1D *histo_muon_invbeta;

  ////////////////////////////////
// Definitions of various inputs // 
  //////////////////////////////
  edm::EDGetTokenT<edm::SimTrackContainer> inputSim;
  edm::EDGetTokenT<edm::SimVertexContainer> inputVtx;
  edm::EDGetTokenT<TrackingParticleCollection> inputTP;
  edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;
  edm::EDGetTokenT<vector<reco::GenParticle> > inputGP;

//My simHits contribution to the code:
  edm::EDGetTokenT<vector<PSimHit>> inputHitsDT;
  edm::EDGetTokenT<vector<PSimHit>> inputHitsCSC;
  edm::EDGetTokenT<vector<PSimHit>> inputHitsRPC;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> theGeometryToken;
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> theDTGeomToken;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> theCSCGeomToken;
  const edm::ESGetToken<RPCGeometry, MuonGeometryRecord> theRPCGeomToken;

};
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

const TrackingParticle & ancestor(const TrackingParticle & particle) {

  const TrackingVertexRef&  tpv = particle.parentVertex(); 
  if (tpv->nSourceTracks() == 0) return particle;
  const TrackingParticle & parent =  **(tpv->sourceTracks_begin());
  return ancestor(parent);
}

Projekt::Projekt(const edm::ParameterSet& conf) 
  : theConfig(conf), theEventCount(0), hscpCount(0), nlines(0), hscpMass(0), theGeometryToken(esConsumes()), theDTGeomToken(esConsumes()), theCSCGeomToken(esConsumes()), theRPCGeomToken(esConsumes())
{
  cout <<" CTORXX" << endl;
//  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(theConfig.getParameter<edm::InputTag>("inputOMTF") );
  inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits")); //czy jest jakas roznica pomiedzy tym a cColl??
  inputVtx =  consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));
//  inputGMT =  consumes< vector<l1t::TrackerMuon> >(edm::InputTag("gmtMuons"));
  inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
  inputGP  =  consumes< vector<reco::GenParticle> >(edm::InputTag("genParticles"));
//Inlusion of simHits, firstly just from DT chambers, now for all chamber types (minus present GEMS)
  inputHitsDT = consumes<vector<PSimHit>>(edm::InputTag("g4SimHits","MuonDTHits"));
  inputHitsCSC = consumes<vector<PSimHit>>(edm::InputTag("g4SimHits","MuonCSCHits"));
  inputHitsRPC = consumes<vector<PSimHit>>(edm::InputTag("g4SimHits","MuonRPCHits"));
  
  
 //Note to self: order of input tags needs to be the same as in edmDumpEventContent

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
//RECOSIM_rhadron_muonsystem.root - combination of both RECOSIM files
//stau432_stau_muonsystem.root - not possible bc file has no simulated tracks
//FEVTSIM_1_stau_M432_full.root - for new FEVTSIM file received on 13/12/2022 (mass 432 MeV, stau) - overwritten by accident by FEVTSIM.root

  myRootFile=new TFile("FEVTSIM_stau_M200_full.root","RECREATE"); //remember to change name when changing datafiles!
  
  hscpTree = new TTree("simTree", "simTree");
  hscpTree->Branch("EV", &(hitData.event), 32000, 99);
  hscpTree->Branch("TR", &(hitData.track), 32000, 99);
  hscpTree->Branch("PID", &(hitData.pid), 32000, 99);
  hscpTree->Branch("SL", &(hitData.superlayer), 32000, 99);
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
  hscpTree->Branch("HITDISTDIFF_G", &(hitData.globalHitDistDiff), 32000, 99);
  hscpTree->Branch("NBETA", &(hitData.betaNextHit), 32000, 99);
  hscpTree->Branch("GEANTTOFDIFF", &(hitData.diffTOFGeant), 32000, 99);
  hscpTree->Branch("CURRENTBTOFDIFF", &(hitData.diffTOFHitB), 32000, 99);
  hscpTree->Branch("NEXTBTOFFDIFF", &(hitData.diffTOFNextB), 32000, 99);

  ////////////////////////////  //////////////////////////////
// PSEUDORAPIDITY HISTOGRAMS SPECIFYING MUON SYSTEM REGIONS //
  //////////////////////////////////////////////////////////

  const Int_t nBins = 20; //this method remains from a previous version of the code, could've been defined within histogram definition...
  Double_t edges[nBins+1] = {-2.1,-1.89,-1.68,-1.47,-1.26,-1.05,-0.84,-0.63,-0.42,-0.21,0,0.21,0.42,0.63,0.84,1.05,1.26,1.47,1.68,1.89,2.1};
  histo_pseudorapidity = new TH1D("histo_pseudorapidity","Charged HSCP count in MTF ranges;Simulated #eta;#Events",nBins,edges);
  histo_pseudorapidity_neutral = new TH1D("histo_pseudorapidity_neutral","Neutral HSCP count in MTF ranges;Simulated #eta;#Events",nBins,edges);

  ////////////////////////////////////////////
// KINEMATIC HISTOGRAMS for  HSCPs AND MUONS //
  ///////////////////////////////////////////

  histo_pdgCount = new TH1D("histo_pdgCount","PID Count;PID;#Events",12,-2000000,2000000);
  histo_pdgCount_HSCP = new TH1D("histo_pdgCount_HSCP","PID Count for HSCP only;PID;#Events",2,-2000000,2000000);

  histo_stau_pt = new TH1D("histo_stau_pt","Generated HSCP p_{T}; Generated p_{T} [GeV]; #Events",30,0.,3000.);
  histo_stau_eta = new TH1D("histo_stau_eta","Generated HSCP #eta; Generated #eta; #Events",100,-5.,5.);
  histo_stau_pl = new TH1D("histo_stau_pl","Generated HSCP p_{L};Genereated p_{L} [GeV]; #Events",60,-3000.,3000.);
  histo_stau_p = new TH1D("histo_stau_p","Generated HSCP p;Generated p [GeV];#Events",50,0.,5000.);
  histo_stau_phi = new TH1D("histo_stau_phi","Generated HSCP #phi;Generated #phi [rad];#Events",37,0.,2*M_PI);
  histo_stau_beta = new TH1D("histo_stau_beta","Generated HSCP #beta;Generated #beta;#Events",20,0.,1.);
  histo_stau_invbeta = new TH1D("histo_stau_invbeta","Generated HSCP 1/#beta;Generated 1/#beta;#Events",40,1.,10.);

  histo_muon_pt = new TH1D("histo_muon_pt","Genereated muon p_{T}; Generated p_{T} [GeV]; #Events",100,0.,100.);
  histo_muon_eta = new TH1D("histo_muon_eta","Generated muon #eta; Generated #eta; #Events",100,-5.,5.);
  histo_muon_pl = new TH1D("histo_muon_pl","Generated muon p_{L};Generated p_{L} [GeV]; #Events",60,-3000.,3000.);
  histo_muon_p = new TH1D("histo_muon_p","Generated muon p;Generated p [GeV];#Events",50,0.,5000.);
  histo_muon_phi = new TH1D("histo_muon_phi","Generated muon #phi;Generated #phi [rad];#Events",37,0.,2*M_PI);
  histo_muon_beta = new TH1D("histo_muon_beta","Generated muon #beta;Generated #beta;#Events",20,0.,1.);
  histo_muon_invbeta = new TH1D("histo_muon_invbeta","Generated muon 1/#beta;Generated 1/#beta;#Events",40,1.,10.);

  

  cout << "HERE Projekt::beginJob()" << endl;
}

void Projekt::endJob()
{
  //write histogram data
  cout << "HSCP COUNT: " << hscpCount << endl; 
  cout << "HSCP COUNT INVERSE: " << 1./hscpCount << endl;
  TH1D *histo_pseudorapidity_ratio = (TH1D*) histo_pseudorapidity->Clone(); //in case of R hadrons, this is for *charged* R hadrons
  histo_pseudorapidity_ratio -> SetTitle("Ratio of charged HSCPs in MTF ranges:total HSCP count;Simulated #eta;HSCP_{MTF}/HSCP_{TOTAL}");
  double scaling = 1./hscpCount; //PROBLEM - gives zero, SOLUTION - 1. not 1 (otherwise any fraction resulting from int/int gives zero)
  histo_pseudorapidity_ratio -> Scale(scaling);
  
  TH1D *histo_pseudorapidity_neutral_ratio = (TH1D*) histo_pseudorapidity_neutral->Clone(); //only relevant for R hadrons
  histo_pseudorapidity_neutral_ratio -> SetTitle("Ratio of neutral HSCPs in MTF ranges:total HSCP count;Simulated #eta;HSCP_{MTF}/HSCP_{TOTAL}");
  histo_pseudorapidity_neutral_ratio -> Scale(scaling);

  histo_pseudorapidity->Write(); 
  histo_pseudorapidity_neutral->Write();
  histo_pseudorapidity_ratio->Write();
  histo_pseudorapidity_neutral_ratio->Write();


  cout << "Event count: " << theEventCount << endl;  
  cout << "Length of nTuple: " << nlines << endl;
  cout << "HSCP mass: " << hscpMass << endl;

  histo_pdgCount->Write();
  histo_pdgCount_HSCP->Write();

  histo_stau_pt->Write();
  histo_stau_eta->Write();
  histo_stau_pl->Write();
  histo_stau_p->Write();
  histo_stau_phi->Write();
  histo_stau_beta->Write();
  histo_stau_invbeta->Write();

  histo_muon_pt->Write();
  histo_muon_eta->Write();
  histo_muon_pl->Write();
  histo_muon_p->Write();
  histo_muon_phi->Write();
  histo_muon_beta->Write();
  histo_muon_invbeta->Write();

  hscpTree->Write();
  myRootFile->Close();
  
  delete histo_pdgCount;
  delete histo_pdgCount_HSCP;
  
  delete histo_pseudorapidity;
  delete histo_pseudorapidity_neutral;
  delete histo_pseudorapidity_ratio;
  delete histo_pseudorapidity_neutral_ratio;

  delete histo_stau_pt;
  delete histo_stau_eta;
  delete histo_stau_pl;
  delete histo_stau_p;
  delete histo_stau_phi;
  delete histo_stau_beta;
  delete histo_stau_invbeta;

  delete histo_muon_pt;
  delete histo_muon_eta;
  delete histo_muon_pl;
  delete histo_muon_p;
  delete histo_muon_phi;
  delete histo_muon_beta;
  delete histo_muon_invbeta;
  
  delete myRootFile;
  cout << "HERE Cwiczenie::endJob()" << endl;
}

void Projekt::analyze(
    const edm::Event& ev, const edm::EventSetup& es){
  cout << " -------------------------------- HERE Cwiczenie::analyze "<< endl;
  cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" useful event count:"<<++theEventCount << endl;
 // bool debug = true;

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

    //cout <<"Number of Gen Particles: "<<genParticles.size() << endl;
  for (const auto & gp : genParticles) {
		if(gp.status()==1){	
			histo_pdgCount->Fill(gp.pdgId());
			//cout<<"type: "<<gp.pdgId()<<" pt_gen: "<<gp.pt()<<" eta_gen: "<<gp.eta()<<" phi_gen: "<<gp.phi()<<endl;
		}
  	if (abs(gp.pdgId()) == 13 && gp.status() == 1){ //generated particle is a muon (which must be stable but that's a given)
			histo_muon_pt->Fill(gp.pt());
      double muon_pl = gp.pt()*sinh(gp.eta());
      histo_muon_pl->Fill(muon_pl);
      double muon_p = gp.pt()*cosh(gp.eta());
			histo_muon_p->Fill(muon_p);
			if(muon_p>2){
			  histo_muon_eta->Fill(gp.eta());
			}
			double muon_phi = gp.phi();
      if (muon_phi <= 0){
				muon_phi += 2*M_PI;
			}
			histo_muon_phi->Fill(muon_phi);
			double muon_beta = muon_p/gp.energy();
			histo_muon_beta->Fill(muon_beta);
      histo_muon_invbeta->Fill(1/muon_beta);
    }
    else if (abs(gp.pdgId())>1000000) { //generated particle is BSM particle (e.g. stau->pdgID=1000015
      //cout << "Particle with PDG ID: " << gp.pdgId() << " and status: " << gp.status() << "\nValid HSCP candidate generated!" << endl;
			if(gp.status()==1){
        hscpMass = gp.mass();
        histo_pdgCount_HSCP->Fill(gp.pdgId());
				histo_stau_pt->Fill(gp.pt());
				double hscp_pl = gp.pt()*sinh(gp.eta());
				histo_stau_pl->Fill(hscp_pl);
				double hscp_p = gp.pt()*cosh(gp.eta());
				histo_stau_p->Fill(hscp_p);
                if(hscp_p>0.5){
				  histo_stau_eta->Fill(gp.eta());
                }
				double hscp_phi = gp.phi();
				if(hscp_phi <= 0){
					hscp_phi += 2*M_PI;
				}
				histo_stau_phi->Fill(hscp_phi);
				double hscp_beta = hscp_p/gp.energy();
        histo_stau_beta->Fill(hscp_beta);
				histo_stau_invbeta->Fill(1/hscp_beta);
     	}
			
      //cout << "\tAssociated Vertex: " << "[" << gp.vx() << "," << gp.vy() << "," << gp.vz() << "]" << endl;
      //cout << "\tNumber of Mothers: " << gp.numberOfMothers() << endl;
      
      //for(long unsigned int i=0; i<gp.numberOfMothers();i++){
       // long unsigned int no = i+1;
        //cout << "\t\tMother #: " << no << ": " << gp.mother(i)->pdgId() << endl;
      //} 
      
      //cout << "\tNumber of Daughters: " << gp.numberOfDaughters() << endl;
      //if(gp.numberOfDaughters()==0){
        //continue;
      //}
      //for(long unsigned int j=0; j<gp.numberOfDaughters(); j++){ 
       // long unsigned int no2 = j+1;
       // cout << "\t\tDaughter #" << no2 << ": " << gp.daughter(j)->pdgId() << endl;
     
    //  cout << endl;
//else{
 //     continue;
    }
	}
   
    /*else{
      cout << "Particle with PDG ID: " << gp.pdgId() << "\nGenerated particle does not qualify as HSCP candidate." << endl;
    }*/

  ///////////////////////////////////
// SimTrack/SimVtx vector definitions //
  //////////////////////////////////

  edm::Handle<edm::SimTrackContainer> simTrk;
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

  const vector<PSimHit> & simCSCHits = ev.get(inputHitsCSC);
  cout <<"Number of simulated CSC hits in event: "<<simCSCHits.size() << endl;

  /*for(const auto &vertex:mySimVerts){
    cout << vertex.parentIndex() << endl;
  }

  for(const auto &track:mySimTracks){
    if(track.trackId()>10) break;
    cout << "Track: " <<track<< endl;
    cout << track.trackId() << ", " << track.vertIndex() << ", " << endl;
    cout << "Position rho = " << mySimVerts[track.vertIndex()].position().Rho() << ", Position z = " << mySimVerts[track.vertIndex()].position().z() << endl;
    cout << "Parent index: " <<  mySimVerts[track.vertIndex()].parentIndex() << endl; 
  }*/
  //////////////
// DT CHAMBERS //
  ////////////

  Int_t hitCount = 0;

  for (vector<PSimHit>::const_iterator iter=simDTHits.begin();iter<simDTHits.end();iter++){ //Iterator is from 0
    
    const PSimHit & hit = *iter;

    if(abs(hit.particleType())>1000000){

//Definition of globsal position of detector in which hit recorded
      //GlobalPoint detPosition = globalGeometry.idToDet(dtDetChamberId)->position();

      DTLayerId dtDetLayerId(hit.detUnitId()); 
      Int_t eventNr = theEventCount;
	    Int_t pid = hit.particleType();
      Int_t trackNr = hit.trackId(); 
      Int_t subDetectorId = dtDetLayerId.subdetId();
      Int_t station = dtDetLayerId.station();   
      Int_t superlayer = dtDetLayerId.superLayer();  
      cout << "Check superlayer: " << superlayer << endl;
      Local3DPoint localPosition = hit.localPosition();
      GlobalPoint globalPosition = globalGeometry.idToDet(dtDetLayerId)->toGlobal(hit.localPosition());
      Double_t r = sqrt((globalPosition.x()*globalPosition.x())+(globalPosition.y()*globalPosition.y()));
      LocalVector p = hit.momentumAtEntry();
      Double_t pt = sqrt((p.y()*p.y())+(p.z()*p.z()));
      Double_t phi = globalPosition.phi();
      if(phi<0){phi = phi + 2*M_PI;}
      Double_t theta = globalPosition.theta();
      Double_t eta = -log(tan(abs(phi)/2));
      Double_t tof = hit.timeOfFlight();
      Double_t distanceIP = globalPosition.mag();

      math::XYZTLorentzVectorD hscpVertexV = mySimVerts[mySimTracks[hit.trackId()-1].vertIndex()].position();
      cout << "HIT VERTEX (x0,y0,z0) = (" << hscpVertexV << ")" << endl;
      Double_t distanceVTX = abs(distanceIP) - abs(hscpVertexV.mag());

      cout << "=========================================================================================" << std::endl;
      cout << "HIT SUMMARY" << endl;
      cout << "simHit: " << hit << endl;
      cout << "Track ID: " << hit.trackId() << " | Det Unit ID: " << hit.detUnitId() << " | PID: " << hit.particleType() << " | p: "<< hit.momentumAtEntry() <<" | phi: " << hit.phiAtEntry() << " | theta: " << hit.thetaAtEntry() << " | TOF: " << hit.timeOfFlight() << endl;
             
  ////////////////////////////
// TOF analysis (preliminary) // 
  ////////////////////////////

	   Double_t betaHit = p.mag()/sqrt((p.mag()*p.mag())+(hscpMass*hscpMass)); //.pabs() calculates the length of the vector pointing to the hit from (0,0,0) 
    															//mass and momenta in GeV (no unit conversion required)

	   cout << "#Beta(hit) = " << betaHit << endl;

    Double_t tofCalculated = distanceIP*1e9*0.01/(betaHit*TMath::C()); //multiplied by 1e9 and 0.01 for tof to be in ns! 
    cout << "HIT STATION: " << station << " || HIT SUPERLAYER: " << superlayer << "|| HIT LAYER: " << dtDetLayerId.layer() << endl;
    cout << "TOF FROM HIT = " << tof << " ns" << " || TOF CALCULATED = " << tofCalculated << " ns" << " || BETA CALCULATED = " << betaHit <<endl;
    cout << "LOCAL POSITION (x,y,z)cm = (" << localPosition.x() << ", " << localPosition.y() << ", " << localPosition.z() << ")" << endl;
    cout << "GLOBAL POSITION (x,y,z)cm = (" << globalPosition.x() << ", " << globalPosition.y() << ", " << globalPosition.z() << ")" << endl;
    cout << "HIT DISTANCE FROM IP = " << distanceIP << " cm" << endl;
    cout << "HIT DISTANCE FROM VTX = " << distanceVTX << " cm" << endl;
    cout << "HIT MOMENTUM = " << p.mag() << " GeV" << endl;   

  //////////////////
// TOF BETWEEN HITS // 
  //////////////////

      if(iter!=simDTHits.end()){
        const PSimHit & nextHit = *(iter+1);
        if(hit.trackId()==nextHit.trackId()){
         DTLayerId nextDTDetLayerId(nextHit.detUnitId());
        Int_t nextSuperLayer = nextDTDetLayerId.superLayer();

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
      hitData.superlayer = superlayer;
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
      hitData.distanceVTX = distanceVTX;
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
      hscpTree->Fill();


        }
      hitCount++; //has to be here in order for hit count to match numbering in hit ntuples :)
}
    }
  }
  //////////////////////////
// Simulated track analysis //
  //////////////////////////

  for (std::vector<SimTrack>::const_iterator it=mySimTracks.begin(); it<mySimTracks.end(); it++) {
    const SimTrack & track = *it;
    if ( track.type() == -99) continue;
    if ( track.vertIndex() != 0) continue;

    double phi_sim = track.momentum().phi(); //momentum azimutal angle
    double pt_sim = track.momentum().pt(); //transverse momentum
    double eta_sim = track.momentum().eta(); //pseudorapidity

    if(track.type()>1000000) hscpCount++;

  }
}
/*    bool muon = false;
    bool matched = false;
    if ( abs(track.type()) == 13 && pt_sim > 1.0) muon = true;
    for (unsigned i=0; i<  gmtMuons.size(); i++) {
        if (   fabs(pt_sim-gmtMuons[i].trkPtr()->momentum().perp()) < 0.5
            && fabs(phi_sim-gmtMuons[i].trkPtr()->momentum().phi()) < 0.1
            && fabs(eta_sim-gmtMuons[i].trkPtr()->momentum().eta()) < 0.1) matched = true;
    }
    if (debug && (muon || matched) ) {
     if (debug) {
      if (muon)    std::cout <<"MUON  "; else std::cout<<"      ";
      if (matched) std::cout <<"MATCH "; else std::cout<<"      "; */
     /* std::cout <<" trackId: " <<track.trackId() 
          <<" type: "<<track.type() // 13 or -13 is a muon
          << " pt_sim: " << pt_sim <<" eta_sim: "<<eta_sim<<" phi_sim: "<<phi_sim
          <<" vtx: "<<track.vertIndex();
      if (track.vertIndex() < static_cast<int>(mySimVerts.size()) ) {
         double vr = mySimVerts[track.vertIndex()].position().Rho();
         double vz = mySimVerts[track.vertIndex()].position().z();
         double z0 = vz-track.momentum().z()/pt_sim*vr;  
         std::cout <<"vert[r,z]: ["<< vr <<", "<< vz <<"], z0: "<<z0
                   <<", parent: "<< mySimVerts[track.vertIndex()].parentIndex();
      }
      std::cout << std::endl; 
    }*/
  
  //if(simMuonCount!=1) { cout<<"    Simulated muon count != 1"<<endl; return; } 

/* std::cout <<" L1 MUONS: "<<std::endl;
  edm::Handle<l1t::RegionalMuonCandBxCollection> l1Omtf;
  ev.getByToken(inputOMTF, l1Omtf);
  int bxNumber = -1;
  for (l1t::RegionalMuonCandBxCollection::const_iterator it = l1Omtf.product()->begin(bxNumber);
       it != l1Omtf.product()->end(bxNumber); ++it) {
    omtfMuonCount++;
    pt_omtf  =  (it->hwPt()-1.)/2.;
    phi_omtf = ( (15.+it->processor()*60.)/360. + it->hwPhi()/576. ) *2*M_PI; 
    if (phi_omtf > 1*M_PI) phi_omtf -=  2*M_PI;
    eta_omtf = it->hwEta()/240.*2.26;
    if (0) std::cout<<" Processor : "<<it->processor()  <<" pT: "<<it->hwPt()<<" eta: "<<it->hwEta()<<" phi: "<<it->hwPhi()<< std::endl;
    std::cout<<" pT: "<< pt_omtf <<" phi: "<<phi_omtf<<" eta: "<<eta_omtf<<std::endl;
  }
  if(omtfMuonCount != 1) cout<<"    OMTF muon count != 1"<<endl;
  if(omtfMuonCount==1) histo->Fill(phi_omtf);
*/
 
/*
  edm::Handle<TrackingVertexCollection> tvColl;
  ev.getByToken(inputTV, tvColl);
  const TrackingVertexCollection & myTV = *(tvColl.product());
  std::cout<<" tracking Vertices: " << myTV.size() << std::endl;

  edm::Handle<TrackingVertexCollection> tv0Coll;
  ev.getByToken(inputTV0, tv0Coll);
  const TrackingVertexCollection & myTV0 = *(tv0Coll.product());
  std::cout<<" initial Vertices: " << myTV0.size() << std::endl;
  for (const auto & tv0 : myTV0 ) {
      std::cout <<" vertex   : "<< tv0.position()<<" nSimVtx: "<< tv0.nG4Vertices() <<" sTk: "<< tv0.nSourceTracks() <<" dTk: "<<tv0.nDaughterTracks() << std::endl; 
  }
*/

/*
  edm::Handle<TrackingParticleCollection> tpColl;
  ev.getByToken(inputTP, tpColl);
  const TrackingParticleCollection & myTP = *(tpColl.product());
  std::cout<<" TRACKING PARTICLES: " << myTP.size() << std::endl;
  for (const auto & tp : myTP) {
    bool gmtTrk = false;
    //for (const auto & l1m : gmtMuons) {
     // if (edm::match(tp, l1m)) gmtTrk=true;
    //}
    if (gmtTrk && tp.pt()>10.) std::cout <<"**MATCH "<<print(tp)<<std::endl;
    if ( abs( tp.pdgId())!=13  || tp.pt() < 1. || tp.parentVertex()->position().Rho()>200. ||  fabs(tp.parentVertex()->position().T()) > 24.) continue;

//    const TrackingVertexRef&  tpv = tp.parentVertex();
//    const reco::GenParticleRefVector& genParticles = tp.genParticles();
//    std::cout << "Muon : " << tp.pt() <<" tracks: "<< tp.g4Tracks().size()<<" genSize: "<<genParticles.size()<< tpv->position()<<std::endl;
//    std::cout <<" TPV   : "<< tpv->position().T()<<" nSimVtx: "<< tpv->nG4Vertices()<<" nGenVtx: "<< tpv->nGenVertices() 
//                           <<" sTk: "<< tpv->nSourceTracks() <<" dTk: "<<tpv->nDaughterTracks(); 
//    if (tpv->g4Vertices().size()>0) std::cout<<" parentIndex: "<< tpv->g4Vertices()[0].parentIndex();
//    std::cout << std::endl; 
      
//    for (TrackingVertex::tp_iterator it=tpv->sourceTracks_begin(); it != tpv->sourceTracks_end();++it) {
//      const TrackingParticle & t = **it;
//      std::cout<< "PARENT: "<< t.pdgId()<<" pt: "<<t.pt()<<" vtx: "<<t.parentVertex()->position()<<std::endl;
//      const TrackingVertexRef & tpv = t.parentVertex();
//    std::cout <<" TPV   : "<< tpv->position().T()<<" nSimVtx: "<< tpv->nG4Vertices()<<" nGenVtx: "<< tpv->nGenVertices() 
//                           <<" sTk: "<< tpv->nSourceTracks() <<" dTk: "<<tpv->nDaughterTracks(); 
//    if (tpv->g4Vertices().size()>0) std::cout<<" parentIndex: "<< tpv->g4Vertices()[0].parentIndex();
//    std::cout << std::endl; 
//        if (t.parentVertex()->g4Vertices().size() >0) {
//        const SimVertex & vtx=t.parentVertex()->g4Vertices()[0];
//   double vr = mySimVerts[track.vertIndex()].position().Rho();
//   double vz = mySimVerts[track.vertIndex()].position().z();
//   double z0 = vz-track.momentum().z()/pt_sim*vr;
//   std::cout <<"vert[r,z]: ["<< vr <<", "<< vz <<"], z0: "<<z0
//             <<", parent: "<< mySimVerts[track.vertIndex()].parentIndex();
//        std::cout <<"Parent SimVtx "<< vtx.position()<<" parentTrack: "<<vtx.parentIndex()<<std::endl;
//        } 
//    }
*/
      
/*    if (debug)std::cout <<"   MUON "<<print(tp)<<std::endl; 
    const TrackingParticle & muAn = ancestor(tp);
    bool matched = false; 
    for (const auto & l1m : gmtMuons) {
      if (match(muAn, l1m)) matched =true;
    }
    if (debug) {
      if (matched) std::cout <<"******* "; else std::cout <<"        ";
      std::cout <<print(muAn) << std::endl;
    }

  }

  //std::cout <<"size of GMT muons: "<< gmtMuons.size() << std::endl;
  //for (const auto & l1m : gmtMuons) {
    const edm::Ptr< TTTrack<Ref_Phase2TrackerDigi_> >&  pTTTrack = l1m.trkPtr();
    if (debug) std::cout <<" HERE GMTz0: "<<l1m.hwZ0()<<" GMTpT: "<<l1m.hwPt()<<" trkZ0: "<<pTTTrack->z0()<< std::endl;
    if (pTTTrack->momentum().perp() > 10.)
    std::cout <<" ###### momentum: "<< pTTTrack->momentum() <<" pT: "<<pTTTrack->momentum().perp()<<" eta: "<< pTTTrack->momentum().eta()<<" phi: "<< pTTTrack->momentum().phi()<<" pca: "<<pTTTrack->POCA()<< std::endl;
  }
  */
  
  
  //write std io
    //std::cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" useful event count:"<<++theEventCount << std::endl;
    

DEFINE_FWK_MODULE(Projekt);

