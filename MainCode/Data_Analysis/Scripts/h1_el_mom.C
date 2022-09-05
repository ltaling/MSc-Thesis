// File to plot histograms of the electron momentum for a simulation. 
//
// Author: Jacob Smith (smithja)
// Creation Date: 11/02/2021
//
// Notes:
// - This file was adapted from Graham Chambers' (gchamber) file called el_mom_plot.C located at
//   /genie/app/users/gchamber/e4nu/tutorial_code/. To better indicate which quantities are being
//   plotted and if the quantities are separated into the interaction components, I have modified
//   this file and created several others to plot similar distributions.

// These include statements were here when I received Graham's code. Delete/Modify them at
// your own peril.
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
#include <string.h>
#include <stdio.h>

#include "parse_file.C"

void h1_el_mom( std::string file_name){
    
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
    TH1F* h1_el_momentum[6];

    h1_el_momentum[0] = (TH1F*)input1->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__5"));
    h1_el_momentum[1] = (TH1F*)input1->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__4"));
    h1_el_momentum[2] = (TH1F*)input1->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__3"));
    h1_el_momentum[3] = (TH1F*)input1->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__2"));
    h1_el_momentum[4] = (TH1F*)input1->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__1"));
    h1_el_momentum[5] = (TH1F*)input1->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__0"));

    // compile all the sectors into one histogram
    for( int i = 1; i < 6; i++) { h1_el_momentum[5]->Add(h1_el_momentum[5-i]); }

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // format the histogram
    h1_el_momentum[5]->Sumw2();
    h1_el_momentum[5]->Rebin( 25);
//    h1_el_momentum[5]->Scale( 1./h1_el_momentum[5]->Integral(),"width"); // normalize the histogram to 1
    h1_el_momentum[5]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
    h1_el_momentum[5]->GetXaxis()->SetTitle("Electron Momentum [GeV/c]");
    h1_el_momentum[5]->GetYaxis()->SetTitle("Events");
    h1_el_momentum[5]->GetXaxis()->SetRangeUser( 0, 3);
    h1_el_momentum[5]->GetXaxis()->CenterTitle( true);
    h1_el_momentum[5]->GetYaxis()->CenterTitle( true);

    h1_el_momentum[5]->SetTitle(TString(info)+TString(cuts));
  
    // draw the histogram
    h1_el_momentum[5]->Draw("HIST");
    h1_el_momentum[5]->Draw("E SAME");

    // crop the margins of the canvas
    c->SetLeftMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();

    // save this histogram as a PDF file
    std::string save_file_path ("invalid");


    save_file_path = "../output/h1_el_mom/"+TString(info)+TString(cuts)+"MottXSecEq1.pdf"; 
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));

    // smithja: this section of code integrates the distribution and prints it to the screen
    int num_events = h1_el_momentum[5]->Integral();

    cout << "// --- PROGRAM OUTPUT ------------------------------------------------------" << endl;
    cout << "There are " << num_events << " many events in the h1_el_mom histogram." << endl;
    cout << "// -------------------------------------------------------------------------" << endl << endl; // separate program output
}
