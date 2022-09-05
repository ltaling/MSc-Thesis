// File to plot the phi and theta projections of the h2_prot_theta_phi histogram 
// as well as the h2_prot_theta_phi histogram itself.
//
// Author: Jacob Smith (smithja)
// Creation Date: 09/28/2021

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
void prot_phiAndTheta_plot( std::string file_name){
    // declare variables that will be used in the plotting/formatting of histograms and file names
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
    TH2F* h2_prot_theta_phi = (TH2F*)input1->Get(TString::Format("h2_prot_theta_phi"));
    TH1D* h1_prot_phi = h2_prot_theta_phi->ProjectionX();
    TH1D* h1_prot_theta = h2_prot_theta_phi->ProjectionY();

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // format the histograms
    h2_prot_theta_phi->SetStats( 0);
    h2_prot_theta_phi->GetZaxis()->SetRangeUser( 0, h2_prot_theta_phi->GetMaximum()); // do not assign a color to bins with <=0 events in them
    h2_prot_theta_phi->GetXaxis()->SetRangeUser( 0, 360); // phi should be between 0 and 360
    h2_prot_theta_phi->GetYaxis()->SetRangeUser( 0, 180); // theta should be between 0 and 360
    h2_prot_theta_phi->GetXaxis()->SetTitle("#phi_{proton} [deg]");
    h2_prot_theta_phi->GetYaxis()->SetTitle("#theta_{proton} [deg]");
    h2_prot_theta_phi->GetXaxis()->CenterTitle( true);
    h2_prot_theta_phi->GetYaxis()->CenterTitle( true);
    h2_prot_theta_phi->SetTitle(TString(info)+TString(cuts));

    h1_prot_phi->SetStats( 0);
    h1_prot_phi->GetXaxis()->SetRangeUser( 0, 360);
    h1_prot_phi->GetXaxis()->SetTitle("#phi_{proton} [deg]");
    h1_prot_phi->GetYaxis()->SetTitle("Events");
    h1_prot_phi->GetXaxis()->CenterTitle( true);
    h1_prot_phi->GetYaxis()->CenterTitle( true);
    h1_prot_phi->SetTitle(TString(info)+TString(cuts));

    h1_prot_theta->SetStats( 0);h1_prot_theta->GetXaxis()->SetRangeUser( 0, 180); 
    h1_prot_theta->GetXaxis()->SetTitle("#theta_{proton} [deg]");
    h1_prot_theta->GetYaxis()->SetTitle("Events");
    h1_prot_theta->GetXaxis()->CenterTitle( true);
    h1_prot_theta->GetYaxis()->CenterTitle( true);
    h1_prot_theta->SetTitle(TString(info)+TString(cuts));

   
    // create the file paths for each of the plots according to the cut_case
    std::string h2_prot_theta_phi_save_path ("invalid");
    std::string h1_prot_theta_save_path ("invalid");
    std::string h1_prot_phi_save_path ("invalid");

    h2_prot_theta_phi_save_path = "../output/h2_prot_theta_phi/"+TString(info)+TString(cuts)+"MottXSecEq1.pdf"; 
    h1_prot_theta_save_path = "../output/h1_prot_theta/"+TString(info)+TString(cuts)+"MottXSecEq1.pdf"; 
    h1_prot_phi_save_path = "../output/h1_prot_phi/"+TString(info)+TString(cuts)+"MottXSecEq1.pdf"; 

  
    // draw the theta-phi histogram and save it as a PDF file
    h2_prot_theta_phi->Draw("colz");
    c->SetLeftMargin( 0.14);
    c->SetRightMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();
    c->SaveAs( TString::Format("%s", h2_prot_theta_phi_save_path.c_str()));

    // draw the phi histogram and save it as a PDF file
    h1_prot_phi->Draw("hist");
    c->SetLeftMargin( 0.14);
    c->SetRightMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();
    c->SaveAs( TString::Format("%s", h1_prot_phi_save_path.c_str()));

    // draw the theta histogram and save it as a PDF file
    h1_prot_theta->Draw("hist");
    c->SetLeftMargin( 0.14);
    c->SetRightMargin( 0.14);
    c->SetBottomMargin(	0.14);
    c->Update();
    c->SaveAs( TString::Format("%s", h1_prot_theta_save_path.c_str()));

    // smithja: this section of code integrates the distributions and prints
    //          it to the screen
    int num_events = h2_prot_theta_phi->Integral();
    cout << "// --- PROGRAM OUTPUT -----------------------------------------------------------------------" << endl;
    cout << "There are " << num_events << " many total events in the h2_prot_theta_phi histogram." << endl;
    cout << "// ------------------------------------------------------------------------------------------" << endl << endl; // separate program output
}
