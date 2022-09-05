// Author: Jacob Smith (smithja)
// Creation Date: 10/1/2021

// These include statements were in the files that I received Graham.
// Delete/Modify them at your own peril.
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

#include "parse_file.C"

void h2_el_prot_theta( std::string file_name){
    
    std::string info ("");
    std::string cuts ("");

    // import the ROOT data file from the output of doing ./genie_analysis in the e4nu directory
    TFile *input1 = TFile::Open( TString( file_name));

    // Grab the TList in the ROOT File for info on the run
    TList *Run_Info = (TList*)input1->Get("Run_Info");

    std::vector<std::string> parse = parse_file(Run_Info);

    info = parse[0];
    cuts = parse[1];


    // histogram initialization
    TH2F* h2_el_prot_theta = (TH2F*)input1->Get(TString::Format("h2_el_prot_theta"));

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // format the histogram
    h2_el_prot_theta->SetStats( 0);
    h2_el_prot_theta->GetZaxis()->SetRangeUser( 0, h2_el_prot_theta->GetMaximum()); // do not assign a color to bins with <=0 events in them
    h2_el_prot_theta->GetXaxis()->SetRangeUser( 20, 60);
    h2_el_prot_theta->GetYaxis()->SetRangeUser( 0, 180);
    h2_el_prot_theta->GetXaxis()->SetTitle("#theta_{electron}");
    h2_el_prot_theta->GetYaxis()->SetTitle("#theta_{proton}");
    h2_el_prot_theta->GetXaxis()->CenterTitle( true);
    h2_el_prot_theta->GetYaxis()->CenterTitle( true);

    h2_el_prot_theta->SetTitle(TString(info)+TString(cuts)); 
    

    // draw the histogram and save it as a PDF file
    h2_el_prot_theta->Draw("colz");
    c->SetLeftMargin( 0.14);
    c->SetRightMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();

    std::string save_file_path ("invalid");
    save_file_path = "../output/h2_el_prot_theta/"+TString(info)+TString(cuts)+"_MottXSecEq1.pdf"; 
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));
}

