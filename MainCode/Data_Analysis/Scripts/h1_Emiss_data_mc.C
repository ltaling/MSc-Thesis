// File to plot histograms of the proton momentum for a simulation.
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
#include <THStack.h>

#include <iostream>
#include <vector>
#include <string.h>
#include <stdio.h>
//#include "myFunctions.cpp"
#include "parse_file.C"

void h1_Emiss_data_mc( int range){

        // declare variables that will be used in the plotting/formatting of histograms and file names
    std::string info1 ("");
    std::string cuts1 ("");
    std::string info2 ("");
    std::string cuts2 ("");

    // import the ROOT data file from the output of doing ./genie_analysis in the e4nu directory
    std::string file_name1;
    std::string file_name2;

    file_name1 = "/genie/app/users/nsteinbe/e4nu/CLAS/GENIE/Second_try/Noemi_run_Genie_1_C12_2.261000.root";
    file_name2 = "/genie/app/users/nsteinbe/e4nu/CLAS/DATA/Second_try/Noemi_run_Data__C12_2.261000.root";
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

    // declare variables that will be used in the plotting/formatting of histograms and file names
    std::string target ("invalid");
    std::string target_plot ("invalid");
    std::string beam_energy ("invalid");
    std::string beam_energy_plot ("invalid");
    std::string elSectors ("invalid");
    std::string elSectors_plot ("invalid");
    std::string deltaPhiEl ("invalid");
    std::string thetaEl_lb ("invalid");
    std::string thetaEl_ub ("invalid");
    std::string elMom_lb ("invalid");
    std::string elMom_lb_plot ("invalid");
    std::string protSectors ("invalid");
    std::string protSectors_plot ("invalid");
    std::string deltaPhiProt ("invalid");
    std::string thetaProt_lb ("invalid");
    std::string thetaProt_ub ("invalid");
    std::string protMom_lb ("invalid");
    std::string protMom_lb_plot ("invalid");

    // histogram initialization
    TH1F* h1_Emiss_mc[4];
    TH1F* h1_Emiss_mc_total;
    TH1F* h1_Emiss_InteractionBreakDown_inSector_mc[4][6];

      // histogram initialization
    TH1F* h1_Emiss_data[6];

    THStack *hs = new THStack("hs","Stacked 1D histograms");

    //if (run_type.find("Genie") !=  string::npos) { // GENIE MC is separated into sectors and interactions
	for (int i = 0; i < 1; i++) { // for all the interactions
	    for (int j = 0; j < 6; j++) { // for all the sectors
    		// extract the histrograms
    		h1_Emiss_InteractionBreakDown_inSector_mc[i][j] = (TH1F*)input1->Get(TString::Format("h1_Int_%i_Sect_%i_Em", i+1, j));
    	    
    	    //Add everything to the first sector
            if(j == 0) continue;
            h1_Emiss_InteractionBreakDown_inSector_mc[i][0]->Add(h1_Emiss_InteractionBreakDown_inSector_mc[i][j]);
	    }

	}

    h1_Emiss_mc_total = h1_Emiss_InteractionBreakDown_inSector_mc[0][0];


    for (int j = 0; j < 1; j++) {
        h1_Emiss_mc[j] = h1_Emiss_InteractionBreakDown_inSector_mc[j][0];
        if(j == 0) continue;
        h1_Emiss_mc_total->Add(h1_Emiss_mc[j]);

    }
	
    for (int i = 0; i < 6; i++) { // for all the sectors
		h1_Emiss_data[i] = (TH1F*)input2->Get( TString::Format("h1_Int_0_Sect_%i_Em", i));
    }


    // compile all the sectors into one histogram
    for( int i = 1; i < 6; i++) { 
       h1_Emiss_data[0]->Add( h1_Emiss_data[i]); 
    }

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    double norm = h1_Emiss_data[0]->Integral()/h1_Emiss_mc_total->Integral();
    std::cout << "normalizing factor = " << norm << "\n";

    const int color_options[4] = {kBlue, kAzure+10, kGreen-3, kViolet};
    // format the histograms
    for(int j = 0; j < 1; j++) {
        h1_Emiss_mc[j]->Sumw2();
        h1_Emiss_mc[j]->Scale(norm); // normalize the mc to the data
        h1_Emiss_mc[j]->SetLineColor(color_options[j]);
        hs->Add(h1_Emiss_mc[j]);
    }
    
      // format the histograms
    h1_Emiss_data[0]->Sumw2();
    h1_Emiss_data[0]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
    h1_Emiss_data[0]->GetXaxis()->SetTitle("Missing Energy [MeV]");
    h1_Emiss_data[0]->GetYaxis()->SetTitle("Events");
    h1_Emiss_data[0]->GetXaxis()->SetRangeUser( 0, 200);
    h1_Emiss_data[0]->GetXaxis()->CenterTitle( true);
    h1_Emiss_data[0]->GetYaxis()->CenterTitle( true);
    h1_Emiss_data[0]->SetTitle(TString(info2)+TString(cuts2));
    h1_Emiss_data[0]->SetMarkerStyle(2);
    h1_Emiss_data[0]->SetMarkerColor(kBlack);
    h1_Emiss_data[0]->GetYaxis()->SetRangeUser(0,2000);

//    h1_prot_momentum[0]->GetXaxis()->SetTitleSize( 0.06);
//    h1_prot_momentum[0]->GetYaxis()->SetTitleSize( 0.06);
    
    h1_Emiss_data[0]->Draw("HIST E");
    hs->Draw("SAME");

    // draw a legend for our plot
    TLegend *legend = new TLegend( 0.2, 0.725, 0.35, 0.825);
    legend->AddEntry( h1_Emiss_mc[0], "QE-CBF SF");
    //legend->AddEntry( h1_Emiss_mc[1], "MEC-GSuSAv2");
    //legend->AddEntry( h1_Emiss_mc[2], "RES-SuSAv2");
    //legend->AddEntry( h1_Emiss_mc[3], "DIS-SuSAv2");
    legend->AddEntry( h1_Emiss_data[0], "CLAS data");
    legend->SetBorderSize( 0);
    legend->SetFillStyle( 0);
    legend->Draw();


    // crop the margins of the canvas
    c->SetLeftMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();

    // save this histogram as a PDF file
    std::string save_file_path ("invalid");
    save_file_path = "../output/h1_prot_mom/"+TString(info1)+TString(cuts1); 
    //TFile *f = new TFile(TString(save_file_path),"RECREATE");
    //cout << "File saved at: " << save_file_path << "\n";
    //h1_prot_momentum[0]->Write();
    //f->Close();
    c->SaveAs(TString(save_file_path)+".pdf");
    // smithja: this section of code integrates the distribution and prints it to the screen
    /*
    int num_events = h1_prot_momentum[0]->Integral();

    cout << "// --- PROGRAM OUTPUT ------------------------------------------------------" << endl;
    cout << "There are " << num_events << " many events in the h1_prot_mom histogram." << endl;
    cout << "// -------------------------------------------------------------------------" << endl << endl; // separate program output
    */
}
