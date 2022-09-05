// File to plot histograms of the proton momentum for a simulation. The histograms generated
// with this macro are separated into the different interaction components: quasi-elastic (QE),
// mesonic current (MEC), resonance (RES), and deep inelastic scattering (DIS). The histograms
// generated with this macro are also separated into FSI and no FSI components. This functionality
// is indicated by "intAndFSIseparated" being in the macro's file name.
//
// Author: Jacob Smith (smithja)
// Creation Date: 10/19/2021

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
//#include "myFunctions.cpp"
#include "parse_file.C"

void h1_prot_mom_intAndFSISeparated( std::string file_name){
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
    TH1F* h1_prot_momentum[48];

    h1_prot_momentum[0] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector5_FSI_prot_mom_QE"));
    h1_prot_momentum[1] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector5_FSI_prot_mom_QE"));
    h1_prot_momentum[2] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector5_FSI_prot_mom_QE"));
    h1_prot_momentum[3] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector5_FSI_prot_mom_QE"));

    h1_prot_momentum[4] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector4_FSI_prot_mom_QE"));
    h1_prot_momentum[5] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector4_FSI_prot_mom_QE"));
    h1_prot_momentum[6] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector4_FSI_prot_mom_QE"));
    h1_prot_momentum[7] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector4_FSI_prot_mom_QE"));

    h1_prot_momentum[8] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector3_FSI_prot_mom_QE"));
    h1_prot_momentum[9] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector3_FSI_prot_mom_QE"));
    h1_prot_momentum[10] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector3_FSI_prot_mom_QE"));
    h1_prot_momentum[11] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector3_FSI_prot_mom_QE"));

    h1_prot_momentum[12] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector2_FSI_prot_mom_QE"));
    h1_prot_momentum[13] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector2_FSI_prot_mom_QE"));
    h1_prot_momentum[14] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector2_FSI_prot_mom_QE"));
    h1_prot_momentum[15] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector2_FSI_prot_mom_QE"));

    h1_prot_momentum[16] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector1_FSI_prot_mom_QE"));
    h1_prot_momentum[17] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector1_FSI_prot_mom_QE"));
    h1_prot_momentum[18] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector1_FSI_prot_mom_QE"));
    h1_prot_momentum[19] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector1_FSI_prot_mom_QE"));

    h1_prot_momentum[20] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector0_FSI_prot_mom_QE"));
    h1_prot_momentum[21] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector0_FSI_prot_mom_QE"));
    h1_prot_momentum[22] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector0_FSI_prot_mom_QE"));
    h1_prot_momentum[23] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector0_FSI_prot_mom_QE"));
// ^---- FSI | no FSI ----v
    h1_prot_momentum[24] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector5_noFSI_prot_mom_QE"));
    h1_prot_momentum[25] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector5_noFSI_prot_mom_QE"));
    h1_prot_momentum[26] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector5_noFSI_prot_mom_QE"));
    h1_prot_momentum[27] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector5_noFSI_prot_mom_QE"));

    h1_prot_momentum[28] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector4_noFSI_prot_mom_QE"));
    h1_prot_momentum[29] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector4_noFSI_prot_mom_QE"));
    h1_prot_momentum[30] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector4_noFSI_prot_mom_QE"));
    h1_prot_momentum[31] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector4_noFSI_prot_mom_QE"));
    
    h1_prot_momentum[32] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector3_noFSI_prot_mom_QE"));
    h1_prot_momentum[33] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector3_noFSI_prot_mom_QE"));
    h1_prot_momentum[34] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector3_noFSI_prot_mom_QE"));
    h1_prot_momentum[35] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector3_noFSI_prot_mom_QE"));

    h1_prot_momentum[36] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector2_noFSI_prot_mom_QE"));
    h1_prot_momentum[37] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector2_noFSI_prot_mom_QE"));
    h1_prot_momentum[38] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector2_noFSI_prot_mom_QE"));
    h1_prot_momentum[39] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector2_noFSI_prot_mom_QE"));

    h1_prot_momentum[40] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector1_noFSI_prot_mom_QE"));
    h1_prot_momentum[41] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector1_noFSI_prot_mom_QE"));
    h1_prot_momentum[42] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector1_noFSI_prot_mom_QE"));
    h1_prot_momentum[43] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector1_noFSI_prot_mom_QE"));

    h1_prot_momentum[44] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector0_noFSI_prot_mom_QE"));
    h1_prot_momentum[45] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector0_noFSI_prot_mom_QE"));
    h1_prot_momentum[46] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector0_noFSI_prot_mom_QE"));
    h1_prot_momentum[47] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector0_noFSI_prot_mom_QE"));

    // compile all the histograms into the four interaction types (QE, MEC, RES, and DIS for quasi-elastic,
    // mesonic current, resonance, and deep inelastic scattering, respectively.
    for( int i = 1; i < 6; i++) {
        h1_prot_momentum[20]->Add(h1_prot_momentum[20-(4*i)]);
	h1_prot_momentum[21]->Add(h1_prot_momentum[21-(4*i)]); 
	h1_prot_momentum[22]->Add(h1_prot_momentum[22-(4*i)]); 
	h1_prot_momentum[23]->Add(h1_prot_momentum[23-(4*i)]); 
    
	h1_prot_momentum[44]->Add(h1_prot_momentum[44-(4*i)]);
        h1_prot_momentum[45]->Add(h1_prot_momentum[45-(4*i)]);
        h1_prot_momentum[46]->Add(h1_prot_momentum[46-(4*i)]);
        h1_prot_momentum[47]->Add(h1_prot_momentum[47-(4*i)]);
    }

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // initialize things such as color options that will be used in plotting the histograms below
    const int color_options[8] = {kBlue, kAzure+10, kGreen-3, kViolet, kOrange+8, kRed+1, kGray, kGreen+3};
    const char * draw_options[8] = {"HIST","HIST SAME","HIST SAME", "HIST SAME", "HIST SAME", "HIST SAME", "HIST SAME", "HIST SAME"};
    const char * draw_options2[8] = {"E SAME", "E SAME", "E SAME", "E SAME", "E SAME", "E SAME", "E SAME", "E SAME"};

    double hist_max = 0.;
    double hist_max_FSI = 0.;
    double hist_max_noFSI = 0.;
    double hist_max_tmp = 0.;
    // format the histograms and then make them
    for( int i = 20; i < 24; i++){
	// formatting statements
        h1_prot_momentum[i]->Sumw2();
        h1_prot_momentum[i]->Rebin( 25);
//        h1_prot_momentum[i]->Scale( 1./h1_prot_momentum[i]->Integral(),"width"); // normalize the histogram to 1
        h1_prot_momentum[i]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
        h1_prot_momentum[i]->SetLineColor( color_options[i-20]);
        h1_prot_momentum[i]->GetXaxis()->SetTitle("Proton Momentum [GeV/c]");
        h1_prot_momentum[i]->GetYaxis()->SetTitle("Events");
        h1_prot_momentum[i]->GetXaxis()->SetRangeUser( 0, 3);
        h1_prot_momentum[i]->GetXaxis()->CenterTitle( true);
        h1_prot_momentum[i]->GetYaxis()->CenterTitle( true);
//        h1_prot_momentum[i]->GetXaxis()->SetTitleSize( 0.06);
//        h1_prot_momentum[i]->GetYaxis()->SetTitleSize( 0.06);

        h1_prot_momentum[i+24]->Sumw2();
        h1_prot_momentum[i+24]->Rebin( 25);
//        h1_prot_momentum[i+24]->Scale( 1./h1_prot_momentum[i+24]->Integral(),"width"); // normalize the histogram to 1
        h1_prot_momentum[i+24]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
        h1_prot_momentum[i+24]->SetLineColor( color_options[i-20+4]);
        h1_prot_momentum[i+24]->GetXaxis()->SetTitle("Proton Momentum [GeV/c]");
        h1_prot_momentum[i+24]->GetYaxis()->SetTitle("Events");
        h1_prot_momentum[i+24]->GetXaxis()->SetRangeUser( 0, 3);
        h1_prot_momentum[i+24]->GetXaxis()->CenterTitle( true);
        h1_prot_momentum[i+24]->GetYaxis()->CenterTitle( true);
//        h1_prot_momentum[i+24]->GetXaxis()->SetTitleSize( 0.06);
//        h1_prot_momentum[i+24]->GetYaxis()->SetTitleSize( 0.06);


        h1_prot_momentum[i]->SetTitle(TString(info)+TString(cuts));
        h1_prot_momentum[i+24]->SetTitle(TString(info)+TString(cuts));

        // take the maximum of the four interactions with and without FSI to use later in plotting
        hist_max_FSI = h1_prot_momentum[i]->GetMaximum();
	hist_max_noFSI = h1_prot_momentum[i+24]->GetMaximum();
	hist_max_tmp = std::max( hist_max_FSI, hist_max_noFSI);
        if (hist_max_tmp > hist_max) { hist_max = hist_max_tmp; }
    }

    // draw the histograms
    for( int i = 20; i < 24; i++){
	h1_prot_momentum[i]->GetYaxis()->SetRangeUser( 0, hist_max*1.1); // smithja: hist_max mult. by 1.1 to give the graph some head room
        h1_prot_momentum[i]->Draw( draw_options[i-20]);
        h1_prot_momentum[i]->Draw( draw_options2[i-20]);
    }
    for( int i = 20; i < 24; i++){
        h1_prot_momentum[i+24]->GetYaxis()->SetRangeUser( 0, hist_max*1.1); // smithja: hist_max mult. by 1.1 to give the graph some head room
        h1_prot_momentum[i+24]->Draw( draw_options[i-16]);
        h1_prot_momentum[i+24]->Draw( draw_options2[i-16]);
    }

    // crop the margins of the canvas
    c->SetLeftMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->Update();

    // create and draw a legend for our histograms
    TLegend *legend = new TLegend( 0.6, 0.510, 1.2, 0.8);
    legend->AddEntry( h1_prot_momentum[20],"QE, FSI","l");
    legend->AddEntry( h1_prot_momentum[21],"MEC, FSI","l");
    legend->AddEntry( h1_prot_momentum[22],"RES, FSI","l");
    legend->AddEntry( h1_prot_momentum[23],"DIS, FSI","l");
    legend->AddEntry( h1_prot_momentum[44],"QE, no FSI","l");
    legend->AddEntry( h1_prot_momentum[45],"MEC, no FSI","l");
    legend->AddEntry( h1_prot_momentum[46],"RES, no FSI","l");
    legend->AddEntry( h1_prot_momentum[47],"DIS, no FSI","l");
    legend->SetBorderSize( 0);
    legend->SetFillStyle( 0);
    legend->Draw();

    // save this histogram as a PDF file
    std::string save_file_path ("invalid");
    save_file_path="../output/h1_prot_mom_intAndFSISeparated/"+TString(info)+TString(cuts)+"MottXSecEq1.pdf";
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));

    // smithja: this section of code integrates the distributions and prints
    //          it to the screen
    int num_QE_FSI_events = h1_prot_momentum[20]->Integral();
    int num_MEC_FSI_events = h1_prot_momentum[21]->Integral();
    int num_RES_FSI_events = h1_prot_momentum[22]->Integral();
    int num_DIS_FSI_events = h1_prot_momentum[23]->Integral();
    int num_QE_noFSI_events = h1_prot_momentum[44]->Integral();
    int num_MEC_noFSI_events = h1_prot_momentum[45]->Integral();
    int num_RES_noFSI_events = h1_prot_momentum[46]->Integral();
    int num_DIS_noFSI_events = h1_prot_momentum[47]->Integral();
    int num_FSI_events = num_QE_FSI_events + num_MEC_FSI_events + num_RES_FSI_events + num_DIS_FSI_events;
    int num_noFSI_events = num_QE_noFSI_events + num_MEC_noFSI_events + num_RES_noFSI_events + num_DIS_noFSI_events;
    int num_events = num_FSI_events + num_noFSI_events;

    cout << "// --- PROGRAM OUTPUT -------------------------------------------------------------------" << endl;
    cout << "There are " << num_QE_FSI_events << " many QE events with FSI in this histogram." << endl;
    cout << "There are " << num_MEC_FSI_events << " many MEC events with FSI in this histogram." << endl;
    cout << "There are " << num_RES_FSI_events << " many RES events with FSI in this histogram." << endl;
    cout << "There are " << num_DIS_FSI_events << " many DIS events with FSI in this histogram." << endl;
    cout << endl;

    cout << "There are " << num_QE_noFSI_events << " many QE events without FSI in this histogram." << endl;
    cout << "There are " << num_MEC_noFSI_events << " many MEC events without FSI in this histogram." << endl;
    cout << "There are " << num_RES_noFSI_events << " many RES events without FSI in this histogram." << endl;
    cout << "There are " << num_DIS_noFSI_events << " many DIS events without FSI in this histogram." << endl;
    cout << endl;

    cout << "There are " << num_FSI_events << " many events with FSI in this histogram." << endl;
    cout << "There are " << num_noFSI_events << " many events without FSI in this histogram." << endl;
    cout << "There are " << num_events << " many total events in the h1_prot_mom_intAndFSISeparated histogram." << endl;
    cout << "// --------------------------------------------------------------------------------------" << endl << endl; // separate program output
}
