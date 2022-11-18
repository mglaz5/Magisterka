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

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TFile.h"
#include <sstream>




using namespace std;


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
  //added variables
  TFile *myRootFile;
  TH1D *histo_pdgCount;
  
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

  edm::EDGetTokenT<edm::SimTrackContainer> inputSim;
  edm::EDGetTokenT<edm::SimVertexContainer> inputVtx;
  edm::EDGetTokenT<TrackingParticleCollection> inputTP;
  edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;
  edm::EDGetTokenT<vector<reco::GenParticle> > inputGP;


};
/*
bool match(const TrackingParticle & tp, const l1t::TrackerMuon & gmt) {
  return (   (fabs(tp.pt()-gmt.trkPtr()->momentum().perp()) )/tp.pt() < 0.1
            && fabs(tp.phi()-gmt.trkPtr()->momentum().phi()) < 0.1
            && fabs(tp.eta()-gmt.trkPtr()->momentum().eta()) < 0.1 );
}*/

std::string print(const TrackingParticle & tp) {
  std::stringstream ss;
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
  : theConfig(conf), theEventCount(0) 
{
  cout <<" CTORXX" << endl;
//  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(theConfig.getParameter<edm::InputTag>("inputOMTF") );
  inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
  inputVtx =  consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));
//  inputGMT =  consumes< vector<l1t::TrackerMuon> >(edm::InputTag("gmtMuons"));
  inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
  inputGP  =  consumes< vector<reco::GenParticle> >(edm::InputTag("genParticles"));
}


Projekt::~Projekt() 
{ 
  cout <<" DTOR" << endl;
}

void Projekt::beginJob()
{
  //make a new Root file
  myRootFile=new TFile("rhadron_mion_analiza.root","RECREATE");
  //create relevant histograms
  histo_pdgCount = new TH1D("histo_pdgCount","PID Count;PID;#Events",12,-3000000,3000000);
  
  histo_stau_pt = new TH1D("histo_stau_pt","Generated R hadron p_{T}; Generated p_{T} [GeV]; #Events",30,0.,3000.);
  histo_stau_eta = new TH1D("histo_stau_eta","Generated R hadron #eta; Generated #eta; #Events",100,-5.,5.);
  histo_stau_pl = new TH1D("histo_stau_pl","Generated R hadron p_{L};Genereated p_{L} [GeV]; #Events",60,-3000.,3000.);
  histo_stau_p = new TH1D("histo_stau_p","Generated R hadron p;Generated p [GeV];#Events",50,0.,5000.);
  histo_stau_phi = new TH1D("histo_stau_phi","Generated R hadron #phi;Generated #phi [rad];#Events",37,0.,2*M_PI);
  histo_stau_beta = new TH1D("histo_stau_beta","Generated R hadron #beta;Generated #beta;#Events",20,0.,1.);
  histo_stau_invbeta = new TH1D("histo_stau_invbeta","Generated R hadron 1/#beta;Generated 1/#beta;#Events",40,1.,10.);

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
  histo_pdgCount->Write();  

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

  myRootFile->Close();
  delete histo_pdgCount;
  
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
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  bool debug = true;

//  edm::Handle<vector<l1t::TrackerMuon> > gmtColl;
//  ev.getByToken(inputGMT, gmtColl);
//  const vector<l1t::TrackerMuon> & gmtMuons = *gmtColl.product();
//  histo->Fill(gmtMuons.size());

/*  edm::Handle<edm::SimTrackContainer> simTrk;
  ev.getByToken(inputSim, simTrk);
  const std::vector<SimTrack>  & mySimTracks = *(simTrk.product());
  std::cout <<" SIMULATED TRACKS: "<<mySimTracks.size()<<std::endl;

  edm::Handle<edm::SimVertexContainer> simVtx;
  ev.getByToken(inputVtx, simVtx);
  const std::vector<SimVertex> & mySimVerts= *(simVtx.product());
  std::cout <<" SIMULATED VERTICES: "<<mySimVerts.size()<<std::endl;



  for (std::vector<SimTrack>::const_iterator it=mySimTracks.begin(); it<mySimTracks.end(); it++) {
    const SimTrack & track = *it;
    if ( track.type() == -99) continue;
    if ( track.vertIndex() != 0) continue;
    sucessful muon, add to countG
    simMuonCount++;

    double phi_sim = track.momentum().phi(); //momentum azimutal angle
    double pt_sim = track.momentum().pt(); //transverse momentum
    double eta_sim = track.momentum().eta(); //pseudorapidity

    bool muon = false;
    bool matched = false;
    if ( abs(track.type()) == 13 && pt_sim > 1.0) muon = true;
    for (unsigned i=0; i<  gmtMuons.size(); i++) {
        if (   fabs(pt_sim-gmtMuons[i].trkPtr()->momentum().perp()) < 0.5
            && fabs(phi_sim-gmtMuons[i].trkPtr()->momentum().phi()) < 0.1
            && fabs(eta_sim-gmtMuons[i].trkPtr()->momentum().eta()) < 0.1) matched = true;
    }
//    if (debug && (muon || matched) ) {
    if (debug) {
      if (muon)    std::cout <<"MUON  "; else std::cout<<"      ";
      if (matched) std::cout <<"MATCH "; else std::cout<<"      "; 
      std::cout <<" trackId: " <<track.trackId() 
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
    } 
  } 
//  if(simMuonCount!=1) { cout<<"    Simulated muon count != 1"<<endl; return; }
*/
  const std::vector<reco::GenParticle> & genParticles = ev.get(inputGP); 
  std::cout <<"Number of Gen Particles: "<<genParticles.size() << std::endl;
  for (const auto & gp : genParticles) {
		if(gp.status()==1){	
			histo_pdgCount->Fill(gp.pdgId());
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
      std::cout << "Particle with PDG ID: " << gp.pdgId() << " and status: " << gp.status() << "\nValid HSCP candidate generated!" << std::endl;
			if(gp.status()==1){
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
			
      std::cout << "\tAssociated Vertex: " << "[" << gp.vx() << "," << gp.vy() << "," << gp.vz() << "]" << std::endl;
      std::cout << "\tNumber of Mothers: " << gp.numberOfMothers() << std::endl;
      
      for(long unsigned int i=0; i<gp.numberOfMothers();i++){
        long unsigned int no = i+1;
        std::cout << "\t\tMother #: " << no << ": " << gp.mother(i)->pdgId() << std::endl;
      } 
      
      std::cout << "\tNumber of Daughters: " << gp.numberOfDaughters() << std::endl;
      if(gp.numberOfDaughters()==0){
        continue;
      }
      for(long unsigned int j=0; j<gp.numberOfDaughters(); j++){ 
        long unsigned int no2 = j+1;
        std::cout << "\t\tDaughter #" << no2 << ": " << gp.daughter(j)->pdgId() << std::endl;
     
      std::cout << std::endl;
	}
}
    
    /*else{
      std::cout << "Particle with PDG ID: " << gp.pdgId() << "\nGenerated particle does not qualify as HSCP candidate." << std::endl;
    }*/
    else{
      continue;
    }
  }

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
//  cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" useful event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Projekt);
