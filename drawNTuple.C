#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH2F.h"  
#include "TAttMarker.h"
#include "TString.h"

void drawNTuple(){

    //std::unique_ptr<TFile> myROOTFile(TFile::Open("FEVTSIM_stau_M200_full.root", "READ"));
    //std::unique_ptr<TNtuple> myNTuple(myROOTFile->Get<TNtuple>("hscpNTuple"));

    //Float_t x,y;
    //myNTuple->SetBranchAddress("x", &x);
    //myNTuple->SetBranchAddress("y", &y);

    //Int_t minEvents = myNTuple->GetMinimum("Event");
    Int_t maxEvents = hscpNTuple->GetMaximum("Event"); //gives maximum value of event, given that first event is 0
    /*std::cout << "Total Events: " << maxEvents << std::endl;
    std::cout<< "First Event: " <<  minEvents << std::endl;

    Int_t nEntries = Int_t(myNTuple->GetEntries());
    std::cout << "NTuple entries: " << nEntries << std::endl;*/

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->Divide(2,1);
    TCanvas *c2 = new TCanvas("c2","c2");
    c2->Divide(2,1);

    c1->cd(1);
    TH2F *histo = new TH2F("histo","histo",16,-800,800,16,-800,800);
    histo->SetMarkerColor(kRed);
    hscpNTuple->Draw("y:x>>histo","Event==1");

    c1->cd(2);
    TH2F *histoL = new TH2F("histoL","histoL",22,-1100,1100,16,-800,800);
    histoL->SetMarkerColor(kRed);
    hscpNTuple->Draw("y:z>>histoL","Event==1");


    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo2","Event==2");
    TH2F *histo2 = (TH2F*) gROOT->FindObject("histo2");
    c1->cd(1);
    histo2->SetMarkerColor(kBlue);
    histo2->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo2L","Event==2");
    TH2F *histo2L = (TH2F*) gROOT->FindObject("histo2L");
    c1->cd(2);
    histo2L->SetMarkerColor(kBlue);
    histo2L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo3","Event==3");
    TH2F *histo3 = (TH2F*) gROOT->FindObject("histo3");
    c1->cd(1);
    histo3->SetMarkerColor(kGreen);
    histo3->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo3L","Event==3");
    TH2F *histo3L = (TH2F*) gROOT->FindObject("histo3L");
    c1->cd(2);
    histo3L->SetMarkerColor(kGreen);
    histo3L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo4","Event==4");
    TH2F *histo4 = (TH2F*) gROOT->FindObject("histo4");
    c1->cd(1);
    histo4->SetMarkerColor(kMagenta);
    histo4->Draw("same");
  
    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo4L","Event==4");
    TH2F *histo4L = (TH2F*) gROOT->FindObject("histo4L");
    c1->cd(2);
    histo4L->SetMarkerColor(kMagenta);
    histo4L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("x:y>>histo5","Event==5");
    TH2F *histo5 = (TH2F*) gROOT->FindObject("histo5");
    c1->cd(1);
    histo5->SetMarkerColor(kOrange);
    histo5->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("x:z>>histo5L","Event==5");
    TH2F *histo5L = (TH2F*) gROOT->FindObject("histo5L");
    c1->cd(2);
    histo5L->SetMarkerColor(kOrange);
    histo5L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo6","Event==6");
    TH2F *histo6 = (TH2F*) gROOT->FindObject("histo6");
    c1->cd(1);
    histo6->SetMarkerColor(30);
    histo6->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo6L","Event==6");
    TH2F *histo6L = (TH2F*) gROOT->FindObject("histo6L");
    c1->cd(2);
    histo6L->SetMarkerColor(30);
    histo6L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo7","Event==7");
    TH2F *histo7 = (TH2F*) gROOT->FindObject("histo7");
    c1->cd(1);
    histo7->SetMarkerColor(kSpring);
    histo7->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo7L","Event==7");
    TH2F *histo7L = (TH2F*) gROOT->FindObject("histo7L");
    c1->cd(2);
    histo7L->SetMarkerColor(kSpring);
    histo7L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo8","Event==8");
    TH2F *histo8 = (TH2F*) gROOT->FindObject("histo8");
    c1->cd(1);
    histo8->SetMarkerColor(kAzure);
    histo8->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo8L","Event==8L");
    TH2F *histo8L = (TH2F*) gROOT->FindObject("histo8L");
    c1->cd(2);
    histo8L->SetMarkerColor(kAzure);
    histo8L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo9","Event==9");
    TH2F *histo9 = (TH2F*) gROOT->FindObject("histo9");
    c1->cd(1);
    histo9->SetMarkerColor(kViolet);
    histo9->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo9L","Event==9L");
    TH2F *histo9L = (TH2F*) gROOT->FindObject("histo9L");
    c1->cd(2);
    histo9L->SetMarkerColor(kViolet);
    histo9L->Draw("same");

    c2->cd(1);
    hscpNTuple->Draw("y:x>>histo10","Event==10");
    TH2F *histo10 = (TH2F*) gROOT->FindObject("histo10");
    c1->cd(1);
    histo10->SetMarkerColor(46);
    histo10->Draw("same");

    c2->cd(2);
    hscpNTuple->Draw("y:z>>histo10L","Event==10");
    TH2F *histo10L = (TH2F*) gROOT->FindObject("histo10L");
    c1->cd(2);
    histo10L->SetMarkerColor(46);
    histo10L->Draw("same");

    /*TH2D *histos[maxEvents];

    for(Int_t i=2;i<maxEvents+1;i++){
      c2->cd(0);
      hscpNTuple->Draw((TString::Format("x:y>>histo%i",i)),(TString::Format("Event==%i",i)));
      histos[i] = (TH2D*) gROOT->FindObject((TString::Format("x:y>>histo%i",i)));
      c1->cd(0);
      histos[i]->SetMarkerColor(10+i);
      histos[i]->Draw("same");
    }*/
  }
    


    

