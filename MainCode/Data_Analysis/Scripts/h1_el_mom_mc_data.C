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

void h1_el_mom_mc_data( int range){
    
    // declare variables that will be used in the plotting/formatting of histograms and file names
    std::string info1 ("");
    std::string cuts1 ("");
    std::string info2 ("");
    std::string cuts2 ("");

    std::string file_name1;
    std::string file_name2;

    if (range == 1) {
        file_name1 = "/genie/app/users/nsteinbe/e4nu/CLAS/GENIE/incl_range1_Genie_1_C12_2.261000.root";
        file_name2 = "/genie/app/users/nsteinbe/e4nu/CLAS/DATA/incl_range1_Data__C12_2.261000.root";
    }
    else if (range == 2) {
        file_name1 = "/genie/app/users/nsteinbe/e4nu/CLAS/GENIE/incl_range2_Genie_1_C12_2.261000.root";
        file_name2 = "/genie/app/users/nsteinbe/e4nu/CLAS/DATA/incl_range2_Data__C12_2.261000.root";
    }
    else {
        file_name1 = "/genie/app/users/nsteinbe/e4nu/CLAS/GENIE/incl_range3_Genie_1_C12_2.261000.root";
        file_name2 = "/genie/app/users/nsteinbe/e4nu/CLAS/DATA/incl_range3_Data__C12_2.261000.root";
    }

    TFile *input1 = TFile::Open( TString(file_name1));
    TFile *input2 = TFile::Open( TString(file_name2));

    // Grab the TList in the ROOT File for info on the run
    TList *Run_Info1 = (TList*)input1->Get("Run_Info");
    TList *Run_Info2 = (TList*)input2->Get("Run_Info");

    std::vector<std::string> parse1 = parse_file(Run_Info1);
    std::vector<std::string> parse2 = parse_file(Run_Info2);

    info1 = parse1[0];
    cuts1 = parse1[1];
    info2 = parse2[0];
    cuts2 = parse2[1];

    std::string run_type1 = std::string(((TNamed*)Run_Info1->FindObject("Run type"))->GetTitle());
    std::string run_type2 = std::string(((TNamed*)Run_Info2->FindObject("Run type"))->GetTitle());

    // histogram initialization
    TH1F* h1_el_momentum_mc[6][4];
    TH1F* h1_el_momentum_data[6];

    for(int i = 0; i < 6; i++) {
        for(int j = 0; j <= 4; j++) {
            h1_el_momentum_mc[i][j] = (TH1F*)input1->Get(TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__%i", j+1, i));
            if (i==0 && j==0) continue;
	    else{
                h1_el_momentum_mc[0][0]->Add(h1_el_momentum_mc[i][j]);
            }
        }
    }
    
    h1_el_momentum_data[0] = (TH1F*)input2->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__5"));
    h1_el_momentum_data[1] = (TH1F*)input2->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__4"));
    h1_el_momentum_data[2] = (TH1F*)input2->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__3"));
    h1_el_momentum_data[3] = (TH1F*)input2->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__2"));
    h1_el_momentum_data[4] = (TH1F*)input2->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__1"));
    h1_el_momentum_data[5] = (TH1F*)input2->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__0"));

    // compile all the sectors into one histogram
    for( int i = 1; i < 6; i++) { 
        h1_el_momentum_data[5]->Add(h1_el_momentum_data[5-i]);
    }

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // format the histogram
    h1_el_momentum_mc[0][0]->Sumw2();
    h1_el_momentum_mc[0][0]->Rebin( 25);
 
    std::cout << "Normalizing factor = "<< h1_el_momentum_data[5]->Integral()/h1_el_momentum_mc[0][0]->Integral() << "\n";
    h1_el_momentum_mc[0][0]->Scale( h1_el_momentum_data[5]->Integral()/h1_el_momentum_mc[0][0]->Integral()); 
    h1_el_momentum_mc[0][0]->SetStats( 0); 
    h1_el_momentum_mc[0][0]->GetXaxis()->SetTitle("Electron Momentum [GeV/c]");
    h1_el_momentum_mc[0][0]->GetYaxis()->SetTitle("Events");
    h1_el_momentum_mc[0][0]->GetXaxis()->SetRangeUser( 0, 3);
    h1_el_momentum_mc[0][0]->GetXaxis()->CenterTitle( true);
    h1_el_momentum_mc[0][0]->GetYaxis()->CenterTitle( true);
    h1_el_momentum_mc[0][0]->SetTitle(TString(info1)+TString(cuts1));

    h1_el_momentum_data[5]->Sumw2();
    h1_el_momentum_data[5]->Rebin( 25);
    h1_el_momentum_data[5]->SetStats( 0); 
    h1_el_momentum_data[5]->GetXaxis()->SetRangeUser( 0, 3);
    h1_el_momentum_data[5]->SetMarkerStyle(2);
    h1_el_momentum_data[5]->SetMarkerColor(kBlack);

  
    // draw the histogram
    h1_el_momentum_mc[0][0]->Draw("HIST E");
    h1_el_momentum_data[5]->Draw("SAME E");
    
    // draw a legend for our plot
    TLegend *legend = new TLegend( 0.2, 0.725, 0.35, 0.825);
    legend->AddEntry( h1_el_momentum_mc[0][0], "GENIE SuSAv2");
    legend->AddEntry( h1_el_momentum_data[5], "CLAS data");
    legend->SetBorderSize( 0);
    legend->SetFillStyle( 0);
    legend->Draw();

    // crop the margins of the canvas
    c->SetLeftMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();

    // save this histogram as a PDF file
    std::string save_file_path ("invalid");


    save_file_path = "../output/h1_el_mom/"+TString(info1)+TString(cuts1)+"MottXSecEq1.pdf"; 
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));

}
