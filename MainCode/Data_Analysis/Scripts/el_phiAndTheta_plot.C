// File to plot the phi and theta projections of the h2_Electron_Theta_Phi histogram 
// as well as the h2_Electron_Theta_Phi histogram itself.
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
#include <TList.h>
#include <TParameter.h>
#include <TNamed.h>

#include <iostream>
#include <vector>
#include <regex>

#include "parse_file.C"


// MAIN FUNCTION ---------------------------------------------------------------------------------------------------------
void el_phiAndTheta_plot( std::string file_name){

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

    // Electron theta cuts
    double theta_el_lb = ((TParameter<double>*)Run_Info->FindObject("El theta lb"))->GetVal();
    double theta_el_ub = ((TParameter<double>*)Run_Info->FindObject("El theta ub"))->GetVal();


    // histogram initialization
    TH2F* h2_Electron_Theta_Phi = (TH2F*)input1->Get(TString::Format("h2_Electron_Theta_Phi"));
    TH1D* h1_el_phi = h2_Electron_Theta_Phi->ProjectionX();
    TH1D* h1_el_theta = h2_Electron_Theta_Phi->ProjectionY();

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // format the histograms
    h2_Electron_Theta_Phi->SetStats( 0);
    h2_Electron_Theta_Phi->GetZaxis()->SetRangeUser( 0, h2_Electron_Theta_Phi->GetMaximum()); // do not assign a color to bins with <=0 events in them
    h2_Electron_Theta_Phi->GetXaxis()->SetRangeUser( 0, 360); // phi should be between 0 and 360
    h2_Electron_Theta_Phi->GetYaxis()->SetRangeUser( 20, 60); // only plot theta values that are simulated
    h2_Electron_Theta_Phi->SetAxisRange(20,60,"Y");
    //h2_Electron_Theta_Phi->SetAxisRange(20,60,"X"); 
   h2_Electron_Theta_Phi->GetXaxis()->SetTitle("#phi_{electron} [deg]");
    h2_Electron_Theta_Phi->GetYaxis()->SetTitle("#theta_{electron} [deg]");
    h2_Electron_Theta_Phi->GetXaxis()->CenterTitle( true);
    h2_Electron_Theta_Phi->GetYaxis()->CenterTitle( true);
    h2_Electron_Theta_Phi->SetTitle(TString(info)+TString(cuts));

    h1_el_phi->SetStats( 0);
    h1_el_phi->GetXaxis()->SetRangeUser( 0, 360);
    h1_el_phi->GetXaxis()->SetTitle("#phi_{electron} [deg]");
    h1_el_phi->GetYaxis()->SetTitle("Events");
    h1_el_phi->GetXaxis()->CenterTitle( true);
    h1_el_phi->GetYaxis()->CenterTitle( true);
    h1_el_phi->SetTitle(TString(info)+TString(cuts));

    h1_el_theta->SetStats( 0);
    h1_el_theta->GetXaxis()->SetRangeUser( theta_el_lb, theta_el_ub); 
    h1_el_theta->GetXaxis()->SetTitle("#theta_{electron} [deg]");
    h1_el_theta->GetYaxis()->SetTitle("Events");
    h1_el_theta->GetXaxis()->CenterTitle( true);
    h1_el_theta->GetYaxis()->CenterTitle( true);
    h1_el_theta->SetTitle(TString(info)+TString(cuts));


    // create the file paths for each of the plots according to the cut_case
    std::string h2_Electron_Theta_Phi_save_path ("invalid");
    std::string h1_el_theta_save_path ("invalid");
    std::string h1_el_phi_save_path ("invalid");

    h2_Electron_Theta_Phi_save_path = "../output/h2_Electron_Theta_Phi/"+TString(info)+TString(cuts)+"MottXSecEq1.pdf"; 
    h1_el_theta_save_path = "../output/h1_el_theta/"+TString(info)+TString(cuts)+"MottXSecEq1_noCuts.pdf"; 
    h1_el_phi_save_path = "../output/h1_el_phi/"+TString(info)+TString(cuts)+"MottXSecEq1_noCuts.pdf"; 


    // draw the theta-phi histogram and save it as a PDF file
    h2_Electron_Theta_Phi->Draw("colz");
    c->SetLeftMargin( 0.14);
    //if (cut_case.compare("noCuts") != 0) { c->SetTopMargin( 0.16); }
    c->SetRightMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();
    c->SaveAs( TString::Format("%s", h2_Electron_Theta_Phi_save_path.c_str()));

    // draw the phi histogram and save it as a PDF file
    h1_el_phi->Draw("hist");
    c->SetLeftMargin( 0.14);
    //if (cut_case.compare("noCuts") != 0) { c->SetTopMargin( 0.16); }
    c->SetRightMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();
    c->SaveAs( TString::Format("%s", h1_el_phi_save_path.c_str()));

    // draw the theta histogram and save it as a PDF file
    h1_el_theta->Draw("hist");
    c->SetLeftMargin( 0.14);
    //if (cut_case.compare("noCuts") != 0) { c->SetTopMargin( 0.16); }
    c->SetRightMargin( 0.14);
    c->SetBottomMargin(	0.14);
    c->Update();
    c->SaveAs( TString::Format("%s", h1_el_theta_save_path.c_str()));

    // smithja: this section of code integrates the distributions and prints
    //          it to the screen
    int num_events = h2_Electron_Theta_Phi->Integral();
    cout << "// --- PROGRAM OUTPUT ---------------------------------------------------------------------------" << endl;
    cout << "There are " << num_events << " many total events in the h2_Electron_Theta_Phi histogram." << endl;
    cout << "// ----------------------------------------------------------------------------------------------" << endl << endl; // separate program output
}
