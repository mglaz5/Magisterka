#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"

void genrecofigures(){
  ////////////////
// GEN V RECO ETA //
  ////////////////

    TCanvas *c1 = new TCanvas("c1","c1",10,64,700,500);
    c1->Range(-3.5,-17.325,3.5,155.925);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    TH1F *genEta = new TH1F("genEta", "genEta", 100, -2.8, 2.8);
    generatedTree->Draw("ETA>>genEta","","goff");
    genEta->SetLineColor(2);
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
    recoEta->SetLineColor(4);
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
    legend1->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4","C");
    legend1->AddEntry(genEta,"Generated candidates","l");
    legend1->AddEntry(recoEta,"Reconstructed (tight) candidates","l");
    legend1->Draw();

    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);
    c1->SaveAs("figures/stau_M432_eta_genreco.png");
    
    
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
    recoPt->SetLineColor(4);
    recoPt->GetXaxis()->SetTitle("p_{T} [GeV]");
    recoPt->GetYaxis()->SetTitle("#Candidates");
    recoPt->GetXaxis()->SetRange(1,100);
    recoPt->GetXaxis()->SetLabelFont(42);
    recoPt->GetXaxis()->SetTitleOffset(1);
    recoPt->GetXaxis()->SetTitleFont(42);
    recoPt->GetYaxis()->SetLabelFont(42);
    recoPt->GetYaxis()->SetTitleFont(42);
    recoPt->Draw("same");

    auto legend2 = new TLegend(0.41,0.52,0.80,0.72);
    legend2->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4","C");
    legend2->AddEntry(genPt,"Generated candidates","l");
    legend2->AddEntry(recoPt,"Reconstructed (tight) candidates","l");
    legend2->Draw();

    c2->Modified();
    c2->cd();
    c2->SetSelected(c2);
    c2->SaveAs("figures/stau_M432_pt_genreco.png");

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
    generatedTree->Draw("PHI>>genPhi","","goff");
    genPhi->SetLineColor(2);
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
    recoTree->Draw("PHI>>recoPhi","","goff");
    recoPhi->SetLineColor(4);
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
    legend3->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4","C");
    legend3->AddEntry(genPhi,"Generated candidates","l");
    legend3->AddEntry(recoPhi,"Reconstructed (tight) candidates","l");
    legend3->Draw();

    c3->Modified();
    c3->cd();
    c3->SetSelected(c3);
    c3->SaveAs("figures/stau_M432_phi_genreco.png");

  ////////////////////////////////////
// GEN V RECO V RECO(VM) INVERSE BETA //
  ////////////////////////////////////  

   /* TCanvas *c4 = new TCanvas("c4","c4",10,64,700,500);
    c4->Range(-3.5,-17.325,3.5,155.925);
    c4->SetFillColor(0);
    c4->SetBorderMode(0);
    c4->SetBorderSize(2);
    c4->SetFrameBorderMode(0);
    c4->SetFrameBorderMode(0);

    TH1F *genInvBeta = new TH1F("genInvBeta", "genInvBeta", 100,0,6.5);
    generatedTree->Draw("INVBETA>>genInvBeta","","goff");
    genInvBeta->SetLineColor(2);
    genInvBeta->SetFillColor(2);
    genInvBeta->GetXaxis()->SetTitle("1/#beta");
    genInvBeta->GetYaxis()->SetTitle("#Entries");
    genInvBeta->GetXaxis()->SetRange(1,100);
    genInvBeta->GetXaxis()->SetLabelFont(42);
    genInvBeta->GetXaxis()->SetTitleOffset(1);
    genInvBeta->GetXaxis()->SetTitleFont(42);
    genInvBeta->GetYaxis()->SetLabelFont(42);
    genInvBeta->GetYaxis()->SetTitleFont(42);
    genInvBeta->SetTitle("");
    genInvBeta->SetStats(false);
    genInvBeta->Draw("");

    TH1F *recoInvBeta = new TH1F("recoInvBeta", "recoInvBeta", 100,0,6.5);
    recoTree->Draw("INVBETA>>recoInvBeta","","goff");
    recoInvBeta->SetLineColor(4);
    recoInvBeta->SetLineWidth(2);
    recoInvBeta->GetXaxis()->SetTitle("1/#beta");
    recoInvBeta->GetYaxis()->SetTitle("#Entries");
    recoInvBeta->GetXaxis()->SetRange(1,100);
    recoInvBeta->GetXaxis()->SetLabelFont(42);
    recoInvBeta->GetXaxis()->SetTitleOffset(1);
    recoInvBeta->GetXaxis()->SetTitleFont(42);
    recoInvBeta->GetYaxis()->SetLabelFont(42);
    recoInvBeta->GetYaxis()->SetTitleFont(42);
    recoInvBeta->Draw("same");

    TH1F *recoInvBetaVM = new TH1F("recoInvBetaVM", "recoInvBetaVM", 100,0,6.5);
    recoTree->Draw("INVBETAEXTRA>>recoInvBetaVM","","goff");
    recoInvBetaVM->SetLineColor(8);
    recoInvBetaVM->SetLineWidth(2);
    recoInvBetaVM->GetXaxis()->SetTitle("1/#beta");
    recoInvBetaVM->GetYaxis()->SetTitle("#Entries");
    recoInvBetaVM->GetXaxis()->SetRange(1,100);
    recoInvBetaVM->GetXaxis()->SetLabelFont(42);
    recoInvBetaVM->GetXaxis()->SetTitleOffset(1);
    recoInvBetaVM->GetXaxis()->SetTitleFont(42);
    recoInvBetaVM->GetYaxis()->SetLabelFont(42);
    recoInvBetaVM->GetYaxis()->SetTitleFont(42);
    recoInvBetaVM->Draw("same");

    auto legend4 = new TLegend(0.39,0.54,0.77,0.74);
    legend4->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4","C");
    legend4->AddEntry(genInvBeta,"Generated candidates","l");
    legend4->AddEntry(recoInvBeta,"Reconstructed (tight) candidates (manual #beta)","l");
    legend4->AddEntry(recoInvBetaVM,"Reconstructed (tight) candidates (reconstructed #beta)","l");
    legend4->Draw();

    c4->Modified();
    c4->cd();
    c4->SetSelected(c4);
    c4->SaveAs("figures/stau_M432_invbeta_genrecoVM.png");    */

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
    c5->SaveAs("figures/stau_M432_deltaR.png");    

  ////////////////////
// MIN DELTA R PLOT  //      
  //////////////////

    TCanvas *c6 = new TCanvas("c6","c6",10,64,700,500);
    c6->Range(-3.5,-17.325,3.5,155.925);
    c6->SetFillColor(0);
    c6->SetBorderMode(0);
    c6->SetBorderSize(2);
    c6->SetFrameBorderMode(0);
    c6->SetFrameBorderMode(0);

    histoMinDeltaR->Draw();
    histoMinDeltaR->SetStats(11);
    histoMinDeltaR->SetLineColor(2);
    histoMinDeltaR->SetLineWidth(2);
    histoMinDeltaR->GetXaxis()->SetTitle("#DeltaR");
    histoMinDeltaR->GetYaxis()->SetTitle("#Entries");
    histoMinDeltaR->GetXaxis()->SetRange(1,205);
    histoMinDeltaR->GetXaxis()->SetLabelFont(42);
    histoMinDeltaR->GetXaxis()->SetTitleOffset(1);
    histoMinDeltaR->GetXaxis()->SetTitleFont(42);
    histoMinDeltaR->GetYaxis()->SetLabelFont(42);
    histoMinDeltaR->GetYaxis()->SetTitleFont(42);
    histoMinDeltaR->SetTitle("");
    histoMinDeltaR->Draw("");   

    c6->Modified();
    c6->cd();
    c6->SetSelected(c6);
    c6->SaveAs("figures/stau_M432_minDeltaR.png");  

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
    invbetaTree->Draw("INVBETAGEN>>genBeta","","goff");
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
    invbetaTree->Draw("INVBETARECO>>recoBeta","","goff");
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
    invbetaTree->Draw("INVBETAVM>>recoBetaVM","","goff");
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
    legend7->SetHeader("Conditions: p_{T} > 50 GeV, |#eta| < 2.4","C");
    legend7->AddEntry(genBeta,"Generated candidates","l");
    legend7->AddEntry(recoBeta,"Reconstructed (tight) candidates (manual #beta)","l");
    legend7->AddEntry(recoBetaVM,"Reconstructed (tight) candidates (reconstructed #beta)","l");
    legend7->Draw();

    c7->Modified();
    c7->cd();
    c7->SetSelected(c7);
    c7->SaveAs("figures/stau_M432_invbeta_genrecoVM.png");  

  ////////////////////////
// 1/beta gen/reco manual //
  ////////////////////////

    /*TCanvas *c8 = new TCanvas("c8","c8",10,64,700,500);
    c8->Range(-3.5,-17.325,3.5,155.925);
    c8->SetFillColor(0);
    c8->SetBorderMode(0);
    c8->SetBorderSize(2);
    c8->SetFrameBorderMode(0);
    c8->SetFrameBorderMode(0);

    TH2F *genrecomanual = new TH2F("genrecomanual","genrecomanual",)


    c8->Modified();
    c8->cd();
    c8->SetSelected(c8);
    c8->SaveAs("figures/stau_M432_invbetadiff.png");   */   
  
  ////////////////////
// 1/beta gen/reco VM // 
  ////////////////////

}
