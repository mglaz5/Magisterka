#include "TCanvas.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TEfficiency.h"

void genrecofigures(){
  ////////////////
// GEN V RECO ETA //
  ////////////////

    /*TCanvas *c1 = new TCanvas("c1","c1",10,64,700,500);
    c1->Range(-3.5,-17.325,3.5,155.925);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    TH1F *genEta = new TH1F("genEta", "genEta", 100, -2.8, 2.8);
    generatedTree->Draw("ETA>>genEta","","goff");
    genEta->SetLineColor(2);
    genEta->SetFillColor(2);
    genEta->GetXaxis()->SetTitle("#eta");
    genEta->GetYaxis()->SetTitle("#Candidates");
    genEta->GetXaxis()->SetRange(1,100);
    genEta->GetXaxis()->SetLabelFont(42);
    genEta->GetXaxis()->SetTitleOffset(1);
    genEta->GetXaxis()->SetTitleFont(42);
    genEta->GetYaxis()->SetLabelFont(42);
    genEta->GetYaxis()->SetTitleFont(42);
    genEta->SetStats(false);
    genEta->SetTitle("");
    genEta->Draw("");

    TH1F *recoEta = new TH1F("recoEta", "recoEta", 100, -2.8, 2.8);
    recoTree->Draw("ETA>>recoEta","","goff");
    recoEta->SetLineColor(8);
    recoEta->SetLineWidth(2);
    recoEta->GetXaxis()->SetTitle("#eta");
    recoEta->GetYaxis()->SetTitle("#Candidates");
    recoEta->GetXaxis()->SetRange(1,100);
    recoEta->GetXaxis()->SetLabelFont(42);
    recoEta->GetXaxis()->SetTitleOffset(1);
    recoEta->GetXaxis()->SetTitleFont(42);
    recoEta->GetYaxis()->SetLabelFont(42);
    recoEta->GetYaxis()->SetTitleFont(42);
    recoEta->Draw("same");

    auto legend1 = new TLegend(0.31,0.21,0.69,0.41);
    legend1->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4, #DeltaR < 0.1","C");
    legend1->AddEntry(genEta,"Generated candidates","l");
    legend1->AddEntry(recoEta,"Reconstructed (tight) candidates","l");
    legend1->Draw();

    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);
    c1->SaveAs("figures_2.4_50/stau_M432_eta_genrecofull.png");
    
    
  ////////////////
// GEN V RECO PT //
  //////////////

    TCanvas *c2 = new TCanvas("c2","c2",10,64,700,500);
    c2->Range(-3.5,-17.325,3.5,155.925);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetBorderSize(2);
    c2->SetFrameBorderMode(0);
    c2->SetFrameBorderMode(0);

    TH1F *genPt = new TH1F("genPt", "genPt", 100,0,2250);
    generatedTree->Draw("PT>>genPt","","goff");
    genPt->SetLineColor(2);
    genPt->SetFillColor(2);
    genPt->GetXaxis()->SetTitle("p_{T} [GeV]");
    genPt->GetYaxis()->SetTitle("#Candidates");
    genPt->GetXaxis()->SetRange(1,100);
    genPt->GetXaxis()->SetLabelFont(42);
    genPt->GetXaxis()->SetTitleOffset(1);
    genPt->GetXaxis()->SetTitleFont(42);
    genPt->GetYaxis()->SetLabelFont(42);
    genPt->GetYaxis()->SetTitleFont(42);
    genPt->SetStats(false);
    genPt->SetTitle("");
    genPt->Draw("");

    TH1F *recoPt = new TH1F("recoPt", "recoPt", 100,0,2250);
    recoTree->Draw("PT>>recoPt","","goff");
    recoPt->SetLineColor(8);
    recoPt->SetLineWidth(2);
    recoPt->GetXaxis()->SetTitle("p_{T} [GeV]");
    recoPt->GetYaxis()->SetTitle("#Candidates");
    recoPt->GetXaxis()->SetRange(1,100);
    recoPt->GetXaxis()->SetLabelFont(42);
    recoPt->GetXaxis()->SetTitleOffset(1);
    recoPt->GetXaxis()->SetTitleFont(42);
    recoPt->GetYaxis()->SetLabelFont(42);
    recoPt->GetYaxis()->SetTitleFont(42);
    recoPt->Draw("same"); #to have stats of both histograms, use sames instead of same

    auto legend2 = new TLegend(0.41,0.52,0.80,0.72);
    legend2->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4, #DeltaR < 0.1","C");
    legend2->AddEntry(genPt,"Generated candidates","l");
    legend2->AddEntry(recoPt,"Reconstructed (tight) candidates","l");
    legend2->Draw();

    c2->Modified();
    c2->cd();
    c2->SetSelected(c2);
    c2->SaveAs("figures_2.4_50/stau_M432_pt_genrecofull.png");

  ////////////////
// GEN V RECO PHI //
  ////////////////

    TCanvas *c3 = new TCanvas("c3","c3",10,64,700,500);
    c3->Range(-3.5,-17.325,3.5,155.925);
    c3->SetFillColor(0);
    c3->SetBorderMode(0);
    c3->SetBorderSize(2);
    c3->SetFrameBorderMode(0);
    c3->SetFrameBorderMode(0);

    TH1F *genPhi = new TH1F("genPhi", "genPhi", 100,0,6.9);
    generatedTree->Draw("PHI>>genPhi","DELTAR<0.1","goff");
    genPhi->SetLineColor(2);
    genPhi->SetFillColor(2);
    genPhi->GetXaxis()->SetTitle("#phi [rad]");
    genPhi->GetYaxis()->SetTitle("#Candidates");
    genPhi->GetXaxis()->SetRange(1,100);
    genPhi->GetXaxis()->SetLabelFont(42);
    genPhi->GetXaxis()->SetTitleOffset(1);
    genPhi->GetXaxis()->SetTitleFont(42);
    genPhi->GetYaxis()->SetLabelFont(42);
    genPhi->GetYaxis()->SetTitleFont(42);
    genPhi->SetTitle("");
    genPhi->SetStats(false);
    genPhi->Draw("");

    TH1F *recoPhi = new TH1F("recoPhi", "recoPhi", 100,0,6.9);
    recoTree->Draw("PHI>>recoPhi","DELTAR<0.1","goff");
    recoPhi->SetLineColor(8);
    recoPhi->SetLineWidth(2);
    recoPhi->GetXaxis()->SetTitle("#phi [rad]");
    recoPhi->GetYaxis()->SetTitle("#Candidates");
    recoPhi->GetXaxis()->SetRange(1,100);
    recoPhi->GetXaxis()->SetLabelFont(42);
    recoPhi->GetXaxis()->SetTitleOffset(1);
    recoPhi->GetXaxis()->SetTitleFont(42);
    recoPhi->GetYaxis()->SetLabelFont(42);
    recoPhi->GetYaxis()->SetTitleFont(42);
    recoPhi->Draw("same");

    auto legend3 = new TLegend(0.24,0.26,0.62,0.46);
    legend3->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4, #DeltaR < 0.1","C");
    legend3->AddEntry(genPhi,"Generated candidates","l");
    legend3->AddEntry(recoPhi,"Reconstructed (tight) candidates","l");
    legend3->Draw();

    c3->Modified();
    c3->cd();
    c3->SetSelected(c3);
    c3->SaveAs("figures_2.4_50/stau_M432_phi_genreco.png");

  //////////////////
// GEN V RECO(VM) P //
  //////////////////  

    TCanvas *c4 = new TCanvas("c4","c4",10,64,700,500);
    c4->Range(-3.5,-17.325,3.5,155.925);
    c4->SetFillColor(0);
    c4->SetBorderMode(0);
    c4->SetBorderSize(2);
    c4->SetFrameBorderMode(0);
    c4->SetFrameBorderMode(0);

    TH1F *genP = new TH1F("genP", "genP", 100,0,4300);
    invbetaTree->Draw("PGEN>>genP","DELTAR<0.1","goff");
    genP->SetLineColor(2);
    genP->SetFillColor(2);
    genP->GetXaxis()->SetTitle("p [GeV]");
    genP->GetYaxis()->SetTitle("#Entries");
    genP->GetXaxis()->SetRange(1,100);
    genP->GetXaxis()->SetLabelFont(42);
    genP->GetXaxis()->SetTitleOffset(1);
    genP->GetXaxis()->SetTitleFont(42);
    genP->GetYaxis()->SetLabelFont(42);
    genP->GetYaxis()->SetTitleFont(42);
    genP->SetTitle("");
    genP->SetStats(false);
    genP->Draw("");

    TH1F *recoP = new TH1F("recoP", "recoP", 100,0,4300);
    invbetaTree->Draw("PRECO>>recoP","DELTAR<0.1","goff");
    recoP->SetLineColor(8);
    recoP->SetLineWidth(2);
    recoP->GetXaxis()->SetTitle("p [GeV]");
    recoP->GetYaxis()->SetTitle("#Entries");
    recoP->GetXaxis()->SetRange(1,100);
    recoP->GetXaxis()->SetLabelFont(42);
    recoP->GetXaxis()->SetTitleOffset(1);
    recoP->GetXaxis()->SetTitleFont(42);
    recoP->GetYaxis()->SetLabelFont(42);
    recoP->GetYaxis()->SetTitleFont(42);
    recoP->Draw("same");

    auto legend4 = new TLegend(0.39,0.54,0.77,0.74);
    legend4->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4, #DeltaR < 0.1","C");
    legend4->AddEntry(genP,"Generated candidates","l");
    legend4->AddEntry(recoP,"Reconstructed (tight) candidates (reconstructed #beta)","l");
    legend4->Draw();

    c4->Modified();
    c4->cd();
    c4->SetSelected(c4);
    c4->SaveAs("figures_2.4_50/stau_M432_p_genreco.png");  

  ////////////////////
// DELTA R PLOT FULL //      
  //////////////////

    TCanvas *c5 = new TCanvas("c5","c5",10,64,700,500);
    c5->Range(-3.5,-17.325,3.5,155.925);
    c5->SetFillColor(0);
    c5->SetBorderMode(0);
    c5->SetBorderSize(2);
    c5->SetFrameBorderMode(0);
    c5->SetFrameBorderMode(0);

    TH1F *fullDeltaR = new TH1F("fullDeltaR", "fullDeltaR", 100,0,6.5);
    invbetaTree->Draw("DELTAR>>fullDeltaR","","goff");
    fullDeltaR->SetLineWidth(2);
    fullDeltaR->GetXaxis()->SetTitle("#DeltaR");
    fullDeltaR->GetYaxis()->SetTitle("#Entries");
    fullDeltaR->GetXaxis()->SetRange(1,100);
    fullDeltaR->GetXaxis()->SetLabelFont(42);
    fullDeltaR->GetXaxis()->SetTitleOffset(1);
    fullDeltaR->GetXaxis()->SetTitleFont(42);
    fullDeltaR->GetYaxis()->SetLabelFont(42);
    fullDeltaR->GetYaxis()->SetTitleFont(42);
    fullDeltaR->Draw("");

    histoDeltaR->SetStats(11);
    histoDeltaR->SetLineColor(2);
    histoDeltaR->SetLineWidth(2);
    histoDeltaR->GetXaxis()->SetTitle("#DeltaR");
    histoDeltaR->GetYaxis()->SetTitle("#Entries");
    histoDeltaR->GetXaxis()->SetRange(1,205);
    histoDeltaR->GetXaxis()->SetLabelFont(42);
    histoDeltaR->GetXaxis()->SetTitleOffset(1);
    histoDeltaR->GetXaxis()->SetTitleFont(42);
    histoDeltaR->GetYaxis()->SetLabelFont(42);
    histoDeltaR->GetYaxis()->SetTitleFont(42);
    histoDeltaR->SetTitle("");
    histoDeltaR->Draw("");

    c5->Modified();
    c5->cd();
    c5->SetSelected(c5);
    c5->SaveAs("figures_2.5_30/stau_M432_deltaR.png");    

  /////////////////////////////////
// GEN V RECO V RECO(VM)  INVBETA //
  ////////////////////////////////  

    TCanvas *c7 = new TCanvas("c7","c7",10,64,700,500);
    c7->Range(-3.5,-17.325,3.5,155.925);
    c7->SetFillColor(0);
    c7->SetBorderMode(0);
    c7->SetBorderSize(2);
    c7->SetFrameBorderMode(0);
    c7->SetFrameBorderMode(0);

    TH1F *genBeta = new TH1F("genBeta", "genBeta", 100,0,6.5);
    invbetaTree->Draw("INVBETAGEN>>genBeta","DELTAR<0.1","goff");
    genBeta->SetLineColor(2);
    genBeta->SetFillColor(2);
    genBeta->GetXaxis()->SetTitle("1/#beta");
    genBeta->GetYaxis()->SetTitle("#Entries");
    genBeta->GetXaxis()->SetRange(1,100);
    genBeta->GetXaxis()->SetLabelFont(42);
    genBeta->GetXaxis()->SetTitleOffset(1);
    genBeta->GetXaxis()->SetTitleFont(42);
    genBeta->GetYaxis()->SetLabelFont(42);
    genBeta->GetYaxis()->SetTitleFont(42);
    genBeta->SetTitle("");
    genBeta->SetStats(0);
    genBeta->Draw("");

    TH1F *recoBeta = new TH1F("recoBeta", "recoBeta", 100,0,6.5);
    invbetaTree->Draw("INVBETARECO>>recoBeta","DELTAR<0.1","goff");
    recoBeta->SetLineColor(4);
    recoBeta->SetLineWidth(2);
    recoBeta->GetXaxis()->SetTitle("1/#beta");
    recoBeta->GetYaxis()->SetTitle("#Entries");
    recoBeta->GetXaxis()->SetRange(1,100);
    recoBeta->GetXaxis()->SetLabelFont(42);
    recoBeta->GetXaxis()->SetTitleOffset(1);
    recoBeta->GetXaxis()->SetTitleFont(42);
    recoBeta->GetYaxis()->SetLabelFont(42);
    recoBeta->GetYaxis()->SetTitleFont(42);
    recoBeta->SetStats(0);
    recoBeta->Draw("same");

    TH1F *recoBetaVM = new TH1F("recoBetaVM", "recoBetaVM", 100,0,6.5);
    invbetaTree->Draw("INVBETAVM>>recoBetaVM","DELTAR<0.1","goff");
    recoBetaVM->SetLineColor(8);
    recoBetaVM->SetLineWidth(2);
    recoBetaVM->GetXaxis()->SetTitle("1/#beta");
    recoBetaVM->GetYaxis()->SetTitle("#Entries");
    recoBetaVM->GetXaxis()->SetRange(1,100);
    recoBetaVM->GetXaxis()->SetLabelFont(42);
    recoBetaVM->GetXaxis()->SetTitleOffset(1);
    recoBetaVM->GetXaxis()->SetTitleFont(42);
    recoBetaVM->GetYaxis()->SetLabelFont(42);
    recoBetaVM->GetYaxis()->SetTitleFont(42);
    recoBetaVM->SetStats(0); 
    recoBetaVM->Draw("same");

    auto legend7 = new TLegend(0.37,0.54,0.75,0.74);
    legend7->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4, #DeltaR < 0.1","C");
    legend7->AddEntry(genBeta,"Generated candidates","l");
    legend7->AddEntry(recoBeta,"Reconstructed (tight) candidates (manual #beta)","l");
    legend7->AddEntry(recoBetaVM,"Reconstructed (tight) candidates (reconstructed #beta)","l");
    legend7->Draw();

    //c7->SetLogx();
    c7->Modified();
    c7->cd();
    c7->SetSelected(c7);
    c7->SaveAs("figures_2.4_50/stau_M432_invbeta_genreco.png");  

  ////////////////////////
// 1/beta gen/reco manual //
  ////////////////////////

  ////////////////////
// GEN V RECO(VM) PT //
  ////////////////// 

    TCanvas *c8 = new TCanvas("c8","c8",10,64,700,500);
    c8->Range(-3.5,-17.325,3.5,155.925);
    c8->SetFillColor(0);
    c8->SetBorderMode(0);
    c8->SetBorderSize(2);
    c8->SetFrameBorderMode(0);
    c8->SetFrameBorderMode(0);

    TH1F *genPT = new TH1F("genPT", "genPT", 100,0,2500);
    invbetaTree->Draw("RECOPT>>genPT","DELTAR<0.1","goff");
    genPT->SetLineColor(2);
    genPT->SetFillColor(2);
    genPT->GetXaxis()->SetTitle("p [GeV]");
    genPT->GetYaxis()->SetTitle("#Entries");
    genPT->GetXaxis()->SetLabelFont(42);
    genPT->GetXaxis()->SetTitleOffset(1);
    genPT->GetXaxis()->SetTitleFont(42);
    genPT->GetYaxis()->SetLabelFont(42);
    genPT->GetYaxis()->SetTitleFont(42);
    genPT->SetTitle("");
    genPT->SetStats(false);
    genPT->Draw("");

    TH1F *recoPT = new TH1F("recoPT", "recoPT", 100,0,2500);
    invbetaTree->Draw("GENPT>>recoPT","DELTAR<0.1","goff");
    recoPT->SetLineColor(8);
    recoPT->SetLineWidth(2);
    recoPT->GetXaxis()->SetTitle("p [GeV]");
    recoPT->GetYaxis()->SetTitle("#Entries");
    recoPT->GetXaxis()->SetLabelFont(42);
    recoPT->GetXaxis()->SetTitleOffset(1);
    recoPT->GetXaxis()->SetTitleFont(42);
    recoPT->GetYaxis()->SetLabelFont(42);
    recoPT->GetYaxis()->SetTitleFont(42);
    recoPT->Draw("same");

    auto legend8 = new TLegend(0.39,0.54,0.77,0.74);
    legend8->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4, #DeltaR < 0.1","C");
    legend8->AddEntry(genPt,"Generated candidates","l");
    legend8->AddEntry(recoPt,"Reconstructed (tight) candidates (reconstructed #beta)","l");
    legend8->Draw();

    c8->Modified();
    c8->cd();
    c8->SetSelected(c8);
    c8->SaveAs("figures_2.4_50/stau_M432_pT_genreco.png");
  ////////////////////////////////////////////////////
// EFFICIENCY PLOTS (pT, p, 1/beta, eta) MATCHED ONLY //     
  ////////////////////////////////////////////////////
    
    TCanvas *c9 = new TCanvas("c9","c9",10,64,700,500);
    c9->Range(-3.5,-17.325,3.5,155.925);
    c9->SetFillColor(0);
    c9->SetBorderMode(0);
    c9->SetBorderSize(2);
    c9->SetFrameBorderMode(0);
    c9->SetFrameBorderMode(0);


    pEff->Write();
    pTEfficiency->Draw();

    c9->Modified();
    c9->cd();
    c9->SetSelected(c9);
    c9->SaveAs("figures_2.4_50/stau_M432_pT_EFFICIENCY.png");
//-----------------------------------------------------------------
    TCanvas *c10 = new TCanvas("c10","c10",10,64,700,500);
    c10->Range(-3.5,-17.325,3.5,155.925);
    c10->SetFillColor(0);
    c10->SetBorderMode(0);
    c10->SetBorderSize(2);
    c10->SetFrameBorderMode(0);
    c10->SetFrameBorderMode(0);

    TH1F *efficiency_p = (TH1F*)recoP->Clone();
    efficiency_p->SetStats(false);
    efficiency_p->Divide(genP);
    efficiency_p->SetLineColor(2);
    efficiency_p->SetLineWidth(2);
    efficiency_p->GetXaxis()->SetTitle("p [GeV]");
    efficiency_p->GetYaxis()->SetTitle("Efficiency");
    efficiency_p->GetXaxis()->SetLabelFont(42);
    efficiency_p->GetXaxis()->SetTitleOffset(1);
    efficiency_p->GetXaxis()->SetTitleFont(42);
    efficiency_p->GetYaxis()->SetLabelFont(42);
    efficiency_p->GetYaxis()->SetTitleFont(42);
    efficiency_p->Draw();
    efficiency_p->SetTitle(" ");

    c10->Modified();
    c10->cd();
    c10->SetSelected(c10);
    c10->SaveAs("figures_2.4_50/stau_M432_p_EFFICIENCY.png");
//-----------------------------------------------------------------
    TCanvas *c11 = new TCanvas("c11","c11",10,64,700,500);
    c11->Range(-3.5,-17.325,3.5,155.925);
    c11->SetFillColor(0);
    c11->SetBorderMode(0);
    c11->SetBorderSize(2);
    c11->SetFrameBorderMode(0);
    c11->SetFrameBorderMode(0);


    c11->Modified();
    c11->cd();
    c11->SetSelected(c11);
    c11->SaveAs("figures_2.4_50/stau_M432_invbeta_EFFICIENCY.png");
//-----------------------------------------------------------------
    TCanvas *c12 = new TCanvas("c12","c12",10,64,700,500);
    c12->Range(-3.5,-17.325,3.5,155.925);
    c12->SetFillColor(0);
    c12->SetBorderMode(0);
    c12->SetBorderSize(2);
    c12->SetFrameBorderMode(0);
    c12->SetFrameBorderMode(0);

    TH1F *genEta2 = new TH1F("genEta2","genEta2",100, -2.8, 2.8);
    invbetaTree->Draw("GENETA>>genEta2","DELTAR<0.1","goff");
    genEta2->SetLineColor(2);
    genEta2->SetFillColor(2);
    genEta2->GetXaxis()->SetTitle("#eta");
    genEta2->GetYaxis()->SetTitle("#Entries");
    genEta2->GetXaxis()->SetLabelFont(42);
    genEta2->GetXaxis()->SetTitleOffset(1);
    genEta2->GetXaxis()->SetTitleFont(42);
    genEta2->GetYaxis()->SetLabelFont(42);
    genEta2->GetYaxis()->SetTitleFont(42);
    genEta2->SetStats(0);
    genEta2->SetTitle(" ");
    genEta2->Draw("");

    TH1F *recoEta2 = new TH1F("recoEta2","recoEta2",100, -2.8, 2.8);
    invbetaTree->Draw("RECOETA>>recoEta2","DELTAR<0.1","goff");
    recoEta2->SetLineColor(8);
    recoEta2->SetLineWidth(2);
    recoEta2->GetXaxis()->SetTitle("#eta");
    recoEta2->GetYaxis()->SetTitle("#Entries");
    recoEta2->GetXaxis()->SetLabelFont(42);
    recoEta2->GetXaxis()->SetTitleOffset(1);
    recoEta2->GetXaxis()->SetTitleFont(42);
    recoEta2->GetYaxis()->SetLabelFont(42);
    recoEta2->GetYaxis()->SetTitleFont(42);
    recoEta2->Draw("same");

    auto legend12 = new TLegend(0.31,0.21,0.69,0.41);
    legend12->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4, #DeltaR < 0.1","C");
    legend12->AddEntry(genEta2,"Generated candidates","l");
    legend12->AddEntry(recoEta2,"Reconstructed (tight) candidates (reconstructed #beta)","l");
    legend12->Draw();
    
    c12->Modified();
    c12->cd();
    c12->SetSelected(c12);
    c12->SaveAs("figures_2.4_50/stau_M432_eta_genreco.png");
//-----------------------------------------------------------------------  
    TCanvas *c13 = new TCanvas("c13","c13",10,64,700,500);
    c13->Range(-3.5,-17.325,3.5,155.925);
    c13->SetFillColor(0);
    c13->SetBorderMode(0);
    c13->SetBorderSize(2);
    c13->SetFrameBorderMode(0);
    c13->SetFrameBorderMode(0);   
    
    c13->Modified();
    c13->cd();
    c13->SetSelected(c13);
    c13->SaveAs("figures_2.4_50/stau_M432_eta_EFFICIENCY.png");

  ///////////////////
// MATCHED STATIONS //   
  //////////////////

    TCanvas *c14 = new TCanvas("c14","c14",10,64,700,500);
    c14->Range(-3.5,-17.325,3.5,155.925);
    c14->SetFillColor(0);
    c14->SetBorderMode(0);
    c14->SetBorderSize(2);
    c14->SetFrameBorderMode(0);
    c14->SetFrameBorderMode(0);
    
    TH2F *histoMuonStationRatio = (TH2F*)histoMuonStationsSelected->Clone("histoMuonStationRatio");
    histoMuonStationRatio->Divide(histoMuonStationsSelected,histoMuonStationsTotal);
    gStyle->SetPaintTextFormat("4.3f");
    histoMuonStationRatio->GetXaxis()->SetTitle("1/#beta_{GEN}");
    histoMuonStationRatio->GetYaxis()->SetTitle("#Stations/Total candidates");
    histoMuonStationRatio->Draw("BOX TEXT90");

    c14->Modified();
    c14->cd();
    c14->SetSelected(c14);
    c14->SaveAs("stau_MatchedStationRatio.png");

  ////////////
// DELTAR L1 //
  //////////

  TCanvas *c15 = new TCanvas("c15", "c15",10,64,700,500);
  c15->Range(-3.5,-17.325,3.5,155.925);
  c15->SetFillColor(0);
  c15->SetBorderMode(0);
  c15->SetBorderSize(2);
  c15->SetFrameBorderMode(0);
  c15->SetFrameBorderMode(0);

  TH1F *deltaRGenL1 = new TH1F("deltaRGenL1","deltaRGenL1",160, 0, 8);
  gmtMuonTree->Draw("DELTAR>>deltaRGenL1","","goff");
  deltaRGenL1->SetLineColor(2);
  deltaRGenL1->SetLineWidth(2);
  deltaRGenL1->SetFillColor(2);
  deltaRGenL1->SetFillStyle(3004);
  deltaRGenL1->GetXaxis()->SetTitle("#DeltaR");
  deltaRGenL1->GetYaxis()->SetTitle("#Entries");
  deltaRGenL1->GetXaxis()->SetLabelFont(42);
  deltaRGenL1->GetXaxis()->SetTitleOffset(1);
  deltaRGenL1->GetXaxis()->SetTitleFont(42);
  deltaRGenL1->GetYaxis()->SetLabelFont(42);
  deltaRGenL1->GetYaxis()->SetTitleFont(42);
  deltaRGenL1->Draw("");
  deltaRGenL1->SetStats(0);

  c15->Modified();
  c15->cd();
  c15->SetSelected(c15);
  c15->SaveAs("deltaR_gentol1.png");

  // LOG! //
  TCanvas *c15_2 = new TCanvas("c15_2", "c15_2",10,64,700,500);
  c15_2->Range(-3.5,-17.325,3.5,155.925);
  c15_2->SetFillColor(0);
  c15_2->SetBorderMode(0);
  c15_2->SetBorderSize(2);
  c15_2->SetFrameBorderMode(0);
  c15_2->SetFrameBorderMode(0);

  deltaRGenL1->Draw("");
  deltaRGenL1->SetStats(0);
  TLine *line = new TLine(0.4,0,0.4,240000);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  c15_2->Modified();
  c15_2->SetLogy();
  c15_2->cd();
  c15_2->SetSelected(c15_2);
  c15_2->SaveAs("deltaR_gentoreco.png");*/

  ////////////
// DELTAR L1 //
  //////////
  TCanvas *c16 = new TCanvas("c16", "c16",10,64,700,500);
  c16->Range(-3.5,-17.325,3.5,155.925);
  c16->SetFillColor(0);
  c16->SetBorderMode(0);
  c16->SetBorderSize(2);
  c16->SetFrameBorderMode(0);
  c16->SetFrameBorderMode(0);

  TH1F *deltaRGenL1 = new TH1F("deltaRGenL1","deltaRGenL1",124, 0, 6.2);
  gmtMuonTree->Draw("DELTAR>>deltaRGenLReco","","goff");
  deltaRGenL1->SetLineColor(2);
  deltaRGenL1->SetLineWidth(2);
  deltaRGenL1->SetFillColor(2);
  deltaRGenL1->SetFillStyle(3004);
  deltaRGenL1->GetXaxis()->SetTitle("#DeltaR_GEN^L1");
  deltaRGenL1->GetYaxis()->SetTitle("#Entries");
  deltaRGenL1->GetXaxis()->SetLabelFont(42);
  deltaRGenL1->GetXaxis()->SetTitleOffset(1);
  deltaRGenL1->GetXaxis()->SetTitleFont(42);
  deltaRGenL1->GetYaxis()->SetLabelFont(42);
  deltaRGenL1->GetYaxis()->SetTitleFont(42);
  deltaRGenL1->Draw("");
  deltaRGenL1->SetStats(0);

  c16->Modified();
  c16->cd();
  c16->SetSelected(c16);
  c16->SaveAs("figures/deltaR_gentoL1.png");

  //////////////////////
// ETA v INVBETA v EFF //
  /////////////////////

    TCanvas *c21 = new TCanvas("c21","c21",10,64,900,500);
    c21->Range(-3.5,-17.325,3.5,155.925);
    c21->SetFillColor(0);
    c21->SetBorderMode(0);
    c21->SetBorderSize(2);
    c21->SetRightMargin(0.13);
    c21->SetFrameBorderMode(0);
    c21->SetFrameBorderMode(0);

    TH1F *InvbetaEtaSelected = (TH1F*)EfficiencyInvbetaEtaSelected->Clone();
    InvbetaEtaSelected->SetStats(false);
    InvbetaEtaSelected->Divide(EfficiencyInvbetaEtaAll);
    InvbetaEtaSelected->SetLineColor(2);
    InvbetaEtaSelected->SetLineWidth(2);
    InvbetaEtaSelected->GetXaxis()->SetTitle("#eta^{GEN}");
    InvbetaEtaSelected->GetYaxis()->SetTitle("1/#beta^{GEN}");
    InvbetaEtaSelected->GetZaxis()->SetTitle("#epsilon");
    InvbetaEtaSelected->GetXaxis()->SetLabelFont(42);
    InvbetaEtaSelected->GetXaxis()->SetTitleOffset(1);
    InvbetaEtaSelected->GetXaxis()->SetTitleFont(42);
    InvbetaEtaSelected->GetYaxis()->SetLabelFont(42);
    InvbetaEtaSelected->GetYaxis()->SetTitleFont(42);
    InvbetaEtaSelected->GetZaxis()->SetLabelFont(42);
    InvbetaEtaSelected->GetZaxis()->SetTitleFont(42);
    InvbetaEtaSelected->GetZaxis()->SetLimits(0.,1.);
    InvbetaEtaSelected->Draw("COLZ");
    InvbetaEtaSelected->SetTitle(" ");

    TLine *lineEtaLowestLower = new TLine(-2.8,1,2.8,1);
    lineEtaLowestLower->SetLineColor(4);
    lineEtaLowestLower->SetLineWidth(2);
    lineEtaLowestLower->SetLineStyle(2);
    lineEtaLowestLower->Draw();

    TLine *lineEtaLowestHigher = new TLine(-2.8,1.1,2.8,1.1);
    lineEtaLowestHigher->SetLineColor(4);
    lineEtaLowestHigher->SetLineWidth(2);
    lineEtaLowestHigher->SetLineStyle(2);
    lineEtaLowestHigher->Draw();

    TLine *lineEtaHighestLower = new TLine(-2.8,1.2,2.8,1.2);
    lineEtaHighestLower->SetLineColor(2);
    lineEtaHighestLower->SetLineWidth(2);
    lineEtaHighestLower->SetLineStyle(2);
    lineEtaHighestLower->Draw();

    TLine *lineEtaHighestHigher = new TLine(-2.8,1.3,2.8,1.3);
    lineEtaHighestHigher->SetLineColor(2);
    lineEtaHighestHigher->SetLineWidth(2);
    lineEtaHighestHigher->SetLineStyle(2);
    lineEtaHighestHigher->Draw();

    auto legend21 = new TLegend(0.45,0.55,0.85,0.8);
    legend21->SetHeader("Full detector #eta range (excluding forward region)","C");
    legend21->AddEntry(lineEtaLowestLower,"Problematic efficiency bin (1/#beta->[1.0,1.1])","l");
    legend21->AddEntry(lineEtaHighestLower,"Highest efficiency bin (1/#beta->[1.2,1.3])","l");
    legend21->Draw();
    c21->Modified();
    c21->cd();
    c21->SetSelected(c21);
    c21->SaveAs("efficiencyInvbetaEtaFULL.png");


  //////////////////////
// PHI v INVBETA v EFF //
  /////////////////////

    TCanvas *c22Barrel = new TCanvas("c22Barrel","c22Barrel",10,64,900,500);
    c22Barrel->Range(-3.5,-17.325,3.5,155.925);
    c22Barrel->SetFillColor(0);
    c22Barrel->SetBorderMode(0);
    c22Barrel->SetBorderSize(2);
    c22Barrel->SetRightMargin(0.13);
    c22Barrel->SetFrameBorderMode(0);
    c22Barrel->SetFrameBorderMode(0);

    TH1F *InvbetaPhiSelectedBarrel = (TH1F*)EfficiencyInvbetaPhiSelectedBarrel->Clone();
    InvbetaPhiSelectedBarrel->SetStats(false);
    InvbetaPhiSelectedBarrel->Divide(EfficiencyInvbetaPhiAllBarrel);
    InvbetaPhiSelectedBarrel->SetLineColor(2);
    InvbetaPhiSelectedBarrel->SetLineWidth(2);
    InvbetaPhiSelectedBarrel->GetXaxis()->SetTitle("#phi^{GEN}");
    InvbetaPhiSelectedBarrel->GetYaxis()->SetTitle("1/#beta");
    InvbetaPhiSelectedBarrel->GetZaxis()->SetTitle("#epsilon");
    InvbetaPhiSelectedBarrel->GetXaxis()->SetLabelFont(42);
    InvbetaPhiSelectedBarrel->GetXaxis()->SetTitleOffset(1);
    InvbetaPhiSelectedBarrel->GetXaxis()->SetTitleFont(42);
    InvbetaPhiSelectedBarrel->GetYaxis()->SetLabelFont(42);
    InvbetaPhiSelectedBarrel->GetYaxis()->SetTitleFont(42);
    InvbetaPhiSelectedBarrel->GetZaxis()->SetLabelFont(42);
    InvbetaPhiSelectedBarrel->GetZaxis()->SetTitleFont(42);
    InvbetaPhiSelectedBarrel->GetZaxis()->SetLimits(0.,1.);
    InvbetaPhiSelectedBarrel->Draw("COLZ");
    InvbetaPhiSelectedBarrel->SetTitle(" ");

    TLine *linePhiLowestLowerBarrel = new TLine(0,1,2*M_PI,1);
    linePhiLowestLowerBarrel->SetLineColor(4);
    linePhiLowestLowerBarrel->SetLineWidth(2);
    linePhiLowestLowerBarrel->SetLineStyle(2);
    linePhiLowestLowerBarrel->Draw();

    TLine *linePhiLowestHigherBarrel = new TLine(0,1.1,2*M_PI,1.1);
    linePhiLowestHigherBarrel->SetLineColor(4);
    linePhiLowestHigherBarrel->SetLineWidth(2);
    linePhiLowestHigherBarrel->SetLineStyle(2);
    linePhiLowestHigherBarrel->Draw();

    TLine *linePhiHighestLowerBarrel= new TLine(0,1.2,2*M_PI,1.2);
    linePhiHighestLowerBarrel->SetLineColor(2);
    linePhiHighestLowerBarrel->SetLineWidth(2);
    linePhiHighestLowerBarrel->SetLineStyle(2);
    linePhiHighestLowerBarrel->Draw();

    TLine *linePhiHighestHigherBarrel = new TLine(0,1.3,2*M_PI,1.3);
    linePhiHighestHigherBarrel->SetLineColor(2);
    linePhiHighestHigherBarrel->SetLineWidth(2);
    linePhiHighestHigherBarrel->SetLineStyle(2);
    linePhiHighestHigherBarrel->Draw();

    auto legend22Barrel = new TLegend(0.4,0.45,0.8,0.7);
    legend22Barrel->SetHeader("Barrel (|#eta|<0.83), #phi binwidth = 30 rad","C");
    legend22Barrel->AddEntry(linePhiLowestLowerBarrel,"Problematic efficiency bin (1/#beta->[1.0,1.1])","l");
    legend22Barrel->AddEntry(linePhiHighestLowerBarrel,"Highest efficiency bin (1/#beta->[1.2,1.3])","l");
    legend22Barrel->Draw();

    c22Barrel->Modified();
    c22Barrel->cd();
    c22Barrel->SetSelected(c22Barrel);
    c22Barrel->SaveAs("efficiencyInvbetaPhiBarrel.png");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    
    TCanvas *c22Endcap = new TCanvas("c22Endcap","c22Endcap",10,64,900,500);
    c22Endcap->Range(-3.5,-17.325,3.5,155.925);
    c22Endcap->SetFillColor(0);
    c22Endcap->SetBorderMode(0);
    c22Endcap->SetBorderSize(2);
    c22Endcap->SetRightMargin(0.13);
    c22Endcap->SetFrameBorderMode(0);
    c22Endcap->SetFrameBorderMode(0);

    TH1F *InvbetaPhiSelectedEndcap = (TH1F*)EfficiencyInvbetaPhiSelectedEndcap->Clone();
    InvbetaPhiSelectedEndcap->SetStats(false);
    InvbetaPhiSelectedEndcap->Divide(EfficiencyInvbetaPhiAllEndcap);
    InvbetaPhiSelectedEndcap->SetLineColor(2);
    InvbetaPhiSelectedEndcap->SetLineWidth(2);
    InvbetaPhiSelectedEndcap->GetXaxis()->SetTitle("#phi^{GEN}");
    InvbetaPhiSelectedEndcap->GetYaxis()->SetTitle("1/#beta^{GEN}");
    InvbetaPhiSelectedEndcap->GetZaxis()->SetTitle("#epsilon");
    InvbetaPhiSelectedEndcap->GetXaxis()->SetLabelFont(42);
    InvbetaPhiSelectedEndcap->GetXaxis()->SetTitleOffset(1);
    InvbetaPhiSelectedEndcap->GetXaxis()->SetTitleFont(42);
    InvbetaPhiSelectedEndcap->GetYaxis()->SetLabelFont(42);
    InvbetaPhiSelectedEndcap->GetYaxis()->SetTitleFont(42);
    InvbetaPhiSelectedEndcap->GetZaxis()->SetLabelFont(42);
    InvbetaPhiSelectedEndcap->GetZaxis()->SetTitleFont(42);
    InvbetaPhiSelectedEndcap->GetZaxis()->SetLimits(0.,1.);
    InvbetaPhiSelectedEndcap->Draw("COLZ");
    InvbetaPhiSelectedEndcap->SetTitle(" ");

    TLine *linePhiLowestLowerEndcap = new TLine(0,1,2*M_PI,1);
    linePhiLowestLowerEndcap->SetLineColor(4);
    linePhiLowestLowerEndcap->SetLineWidth(2);
    linePhiLowestLowerEndcap->SetLineStyle(2);
    linePhiLowestLowerEndcap->Draw();

    TLine *linePhiLowestHigherEndcap = new TLine(0,1.1,2*M_PI,1.1);
    linePhiLowestHigherEndcap->SetLineColor(4);
    linePhiLowestHigherEndcap->SetLineWidth(2);
    linePhiLowestHigherEndcap->SetLineStyle(2);
    linePhiLowestHigherEndcap->Draw();

    TLine *linePhiHighestLowerEndcap= new TLine(0,1.2,2*M_PI,1.2);
    linePhiHighestLowerEndcap->SetLineColor(2);
    linePhiHighestLowerEndcap->SetLineWidth(2);
    linePhiHighestLowerEndcap->SetLineStyle(2);
    linePhiHighestLowerEndcap->Draw();

    TLine *linePhiHighestHigherEndcap = new TLine(0,1.3,2*M_PI,1.3);
    linePhiHighestHigherEndcap->SetLineColor(2);
    linePhiHighestHigherEndcap->SetLineWidth(2);
    linePhiHighestHigherEndcap->SetLineStyle(2);
    linePhiHighestHigherEndcap->Draw();

    auto legend22Endcap = new TLegend(0.4,0.45,0.8,0.7);
    legend22Endcap->SetHeader("Endcap (|#eta|>1.24), #phi bindwith = 20 rad","C");
    legend22Endcap->AddEntry(linePhiLowestLowerEndcap,"Problematic efficiency bin (1/#beta->[1.0,1.1])","l");
    legend22Endcap->AddEntry(linePhiHighestLowerEndcap,"Highest efficiency bin (1/#beta->[1.2,1.3])","l");
    legend22Endcap->Draw();

    c22Endcap->Modified();
    c22Endcap->cd();
    c22Endcap->SetSelected(c22Endcap);
    c22Endcap->SaveAs("efficiencyInvbetaPhiEndcap.png");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    TCanvas *c22Overlap = new TCanvas("c22Overlap","c22Overlap",10,64,900,500);
    c22Overlap->Range(-3.5,-17.325,3.5,155.925);
    c22Overlap->SetFillColor(0);
    c22Overlap->SetBorderMode(0);
    c22Overlap->SetBorderSize(2);
    c22Overlap->SetRightMargin(0.13);
    c22Overlap->SetFrameBorderMode(0);
    c22Overlap->SetFrameBorderMode(0);

    TH1F *InvbetaPhiSelectedOverlap = (TH1F*)EfficiencyInvbetaPhiSelectedOverlap->Clone();
    InvbetaPhiSelectedOverlap->SetStats(false);
    InvbetaPhiSelectedOverlap->Divide(EfficiencyInvbetaPhiAllOverlap);
    InvbetaPhiSelectedOverlap->SetLineColor(2);
    InvbetaPhiSelectedOverlap->SetLineWidth(2);
    InvbetaPhiSelectedOverlap->GetXaxis()->SetTitle("#phi^{GEN}");
    InvbetaPhiSelectedOverlap->GetZaxis()->SetTitle("#epsilon");
    InvbetaPhiSelectedOverlap->GetXaxis()->SetLabelFont(42);
    InvbetaPhiSelectedOverlap->GetXaxis()->SetTitleOffset(1);
    InvbetaPhiSelectedOverlap->GetXaxis()->SetTitleFont(42);
    InvbetaPhiSelectedOverlap->GetYaxis()->SetLabelFont(42);
    InvbetaPhiSelectedOverlap->GetYaxis()->SetTitleFont(42);
    InvbetaPhiSelectedOverlap->GetZaxis()->SetLabelFont(42);
    InvbetaPhiSelectedOverlap->GetZaxis()->SetTitleFont(42);
    InvbetaPhiSelectedOverlap->GetZaxis()->SetLimits(0.,1.);
    InvbetaPhiSelectedOverlap->Draw("COLZ");
    InvbetaPhiSelectedOverlap->SetTitle(" ");

    TLine *linePhiLowestLowerOverlap = new TLine(0,1,2*M_PI,1);
    linePhiLowestLowerOverlap->SetLineColor(4);
    linePhiLowestLowerOverlap->SetLineWidth(2);
    linePhiLowestLowerOverlap->SetLineStyle(2);
    linePhiLowestLowerOverlap->Draw();

    TLine *linePhiLowestHigherOverlap = new TLine(0,1.1,2*M_PI,1.1);
    linePhiLowestHigherOverlap->SetLineColor(4);
    linePhiLowestHigherOverlap->SetLineWidth(2);
    linePhiLowestHigherOverlap->SetLineStyle(2);
    linePhiLowestHigherOverlap->Draw();

    TLine *linePhiHighestLowerOverlap= new TLine(0,1.2,2*M_PI,1.2);
    linePhiHighestLowerOverlap->SetLineColor(2);
    linePhiHighestLowerOverlap->SetLineWidth(2);
    linePhiHighestLowerOverlap->SetLineStyle(2);
    linePhiHighestLowerOverlap->Draw();

    TLine *linePhiHighestHigherOverlap = new TLine(0,1.3,2*M_PI,1.3);
    linePhiHighestHigherOverlap->SetLineColor(2);
    linePhiHighestHigherOverlap->SetLineWidth(2);
    linePhiHighestHigherOverlap->SetLineStyle(2);
    linePhiHighestHigherOverlap->Draw();

    auto legend22Overlap = new TLegend(0.4,0.45,0.8,0.7);
    legend22Overlap->SetHeader("Overlap (0.83<|#eta|<1.24), #phi binwidth = 60 rad","C");
    legend22Overlap->AddEntry(linePhiLowestLowerOverlap,"Problematic efficiency bin (1/#beta->[1.0,1.1])","l");
    legend22Overlap->AddEntry(linePhiHighestLowerOverlap,"Highest efficiency bin (1/#beta->[1.2,1.3])","l");
    legend22Overlap->Draw();

    c22Overlap->Modified();
    c22Overlap->cd();
    c22Overlap->SetSelected(c22Overlap);
    c22Overlap->SaveAs("efficiencyInvbetaPhiOverlap.png");

/////////////////////////////////////////////////////////////////////////////////////////////////////   

    TCanvas *c23 = new TCanvas("c23","c23",10,64,900,500);
    c23->Range(-3.5,-17.325,3.5,155.925);
    c23->SetFillColor(0);
    c23->SetBorderMode(0);
    c23->SetBorderSize(2);
    c23->SetRightMargin(0.13);
    c23->SetFrameBorderMode(0);
    c23->SetFrameBorderMode(0);

    TH2F *plot = (TH2F*)EfficiencyInvbetaEtaAll->Clone();
    plot->SetStats(false);
    plot->SetLineColor(2);
    plot->SetLineWidth(2);
    plot->GetXaxis()->SetTitle("#eta");
    plot->GetXaxis()->SetRangeUser(-4.,4.);
    plot->GetYaxis()->SetRangeUser(0.,8.);
    plot->GetYaxis()->SetTitle("1/#beta");
    plot->GetZaxis()->SetTitle("#Candidates");
    plot->GetZaxis()->SetLimits(0.,1.);
    plot->GetXaxis()->SetLabelFont(42);
    plot->GetXaxis()->SetTitleOffset(1);
    plot->GetXaxis()->SetTitleFont(42);
    plot->GetYaxis()->SetLabelFont(42);
    plot->GetYaxis()->SetTitleFont(42);
    plot->GetZaxis()->SetLabelFont(42);
    plot->GetZaxis()->SetTitleFont(42);
    plot->Draw("COLZ");
    plot->SetTitle(" ");

    c23->Modified();
    c23->cd();
    c23->SetSelected(c23);
    c23->SaveAs("EfficiencyInvbetaEtaAll.png");

/////////////////////////////////////////////////////////////////////////////////////////////////////

}

