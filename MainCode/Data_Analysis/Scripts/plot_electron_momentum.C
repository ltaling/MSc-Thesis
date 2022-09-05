// This is a file to run through the example ROOT file included in Afro Papadopoulou's tutorial on e4nu. After
// following Afro's tutorial and cloning the github repository, we import the file that is on page four (4) of
// Afro's slides: /pnfs/genie/persistent/users/apapadop/e4v_Workshop_C12_1GeV_TestSample.root.
//
// Author: Jacob A. Smith
// Date of Creation: 8/26/2021 11:27am US Eastern


#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void plot_electron_momentum(){
    TFile *f = new TFile("/pnfs/genie/persistent/users/apapadop/e4v_Workshop_C12_1GeV_TestSample.root");
    TTree *t1 = (TTree*)f->Get("gst");
    Double_t pl;
    bool qel, mec, res, dis; // quasi-elastic, mesonic current, resonance, and deep inelastic scattering, resp.
    TCanvas *c = new TCanvas();

    t1->SetBranchAddress("pl", &pl);
    t1->SetBranchAddress("qel", &qel);
    t1->SetBranchAddress("mec", &mec);
    t1->SetBranchAddress("res", &res);
    t1->SetBranchAddress("dis", &dis);

    // create the histograms
    TH1F *hpl = new TH1F("hpl", "Lepton Momentum Distribution", 50, 0, 1.2);
    TH1F *hpl_qel = new TH1F("hpl_qel", "QEL", 50, 0, 1.2);
    TH1F *hpl_mec = new TH1F("hpl_mec", "MEC", 50, 0, 1.2);
    TH1F *hpl_res = new TH1F("hpl_res", "RES", 50, 0, 1.2);
    TH1F *hpl_dis = new TH1F("hpl_dis", "DIS", 50, 0, 1.2);

    // color code histograms
    hpl_qel->SetLineColor( kBlack);
    hpl_mec->SetLineColor( kGreen);
    hpl_res->SetLineColor( kRed);
    hpl_dis->SetLineColor( kMagenta);

    Long64_t n_entries = t1->GetEntries();

    // for all the entries in the ROOT file, calculate and histogram the lepton (electron) momentum
    for( Long64_t i = 0; i < n_entries; i++){
        t1->GetEntry( i);

        if( qel) hpl_qel->Fill( pl);
        else if( mec) hpl_mec->Fill( pl);
        else if( res) hpl_res->Fill( pl);
        else hpl_dis->Fill( pl); // if the event is not qel, mec, or res, it must be dis 

        // fill the general Lepton Momentum Distribution
        hpl->Fill( pl);
    }

    // draw the histogram
    hpl->Draw("hist"); // this might not display the histogram because of permissions errors on the GENIE GPVM
    hpl_qel->Draw("hist same");
    hpl_mec->Draw("hist same");
    hpl_res->Draw("hist same");
    hpl_dis->Draw("hist same");    

    // create a legend for the histogram
    TLegend *legend = new TLegend( 0.775, 0.45, 0.9, 0.75);
    legend->AddEntry( hpl_qel, "QEL", "l");
    legend->AddEntry( hpl_mec, "MEC", "l");
    legend->AddEntry( hpl_res, "RES", "l");
    legend->AddEntry( hpl_dis, "DIS", "l");
    legend->Draw();

    c->SaveAs("/genie/app/users/smithja/output/plot_electron_momentum/hpl_separated.pdf");
}
