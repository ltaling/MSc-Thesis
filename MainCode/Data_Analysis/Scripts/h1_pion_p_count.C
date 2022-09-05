// This macro makes the pion energy plots

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

#include <iostream>
#include <vector>
#include <string.h>
#include <stdio.h>

#include "parse_file.C"

void h1_pion_p_count (std::string file_name, bool pimin, bool pi0, bool piplu){
  // declare variables that will be used in the plotting/formatting of histograms and file names
  std::string info ("");
  std::string cuts ("");

  // import the ROOT data file from the output of doing ./genie_analysis in the e4nu directory
  TFile *input1 = TFile::Open( TString( file_name));

  // Grab the TList in the ROOT File for info on the run
  TList *Run_Info = (TList*)input1->Get("Run_Info");
  gStyle->SetOptStat(0);

  std::vector<std::string> parse = parse_file(Run_Info);
  std:string types = "";

  info = parse[0];
  cuts = parse[1];

  // Histogram initialisation
  TH1F* h1_p_stacked_npi[6]; // Bins of all pion types added together
  TH1F* h1_pi_p_0_npi[6]; //pi-
  TH1F* h1_pi_p_1_npi[6]; //pi0
  TH1F* h1_pi_p_2_npi[6]; //pi+
  TH1F* h1_toPlot_npi[4]; //things to actually plot

  TH1F* h1_p_stacked_1pi[6]; // Bins of all pion types added together
  TH1F* h1_pi_p_0_1pi[6]; //pi-
  TH1F* h1_pi_p_1_1pi[6]; //pi0
  TH1F* h1_pi_p_2_1pi[6]; //pi+
  TH1F* h1_toPlot_1pi[4]; //things to actually plot
  // Drawing histograms
  if (pimin){
    h1_pi_p_0_npi[0] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_1_type_0")); //QES
    h1_pi_p_0_npi[1] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_2_type_0")); //MEC
    h1_pi_p_0_npi[2] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_3_type_0")); //OtherRES
    h1_pi_p_0_npi[3] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_4_type_0")); //DIS
    h1_pi_p_0_npi[4] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_5_type_0")); //1232RES
    h1_pi_p_0_npi[5] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_0_type_0")); //Other
    h1_pi_p_0_1pi[0] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_1_type_0")); //QES
    h1_pi_p_0_1pi[1] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_2_type_0")); //MEC
    h1_pi_p_0_1pi[2] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_3_type_0")); //OtherRES
    h1_pi_p_0_1pi[3] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_4_type_0")); //DIS
    h1_pi_p_0_1pi[4] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_5_type_0")); //1232RES
    h1_pi_p_0_1pi[5] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_0_type_0")); //Other
    if (h1_p_stacked_npi[2] == nullptr){
      types = types + "#pi^{-}";
      h1_p_stacked_npi[0] = h1_pi_p_0_npi[0];
      h1_p_stacked_npi[1] = h1_pi_p_0_npi[1];
      h1_p_stacked_npi[2] = h1_pi_p_0_npi[2];
      h1_p_stacked_npi[3] = h1_pi_p_0_npi[3];
      h1_p_stacked_npi[4] = h1_pi_p_0_npi[4];
      h1_p_stacked_npi[5] = h1_pi_p_0_npi[5];
      h1_p_stacked_1pi[0] = h1_pi_p_0_1pi[0];
      h1_p_stacked_1pi[1] = h1_pi_p_0_1pi[1];
      h1_p_stacked_1pi[2] = h1_pi_p_0_1pi[2];
      h1_p_stacked_1pi[3] = h1_pi_p_0_1pi[3];
      h1_p_stacked_1pi[4] = h1_pi_p_0_1pi[4];
      h1_p_stacked_1pi[5] = h1_pi_p_0_1pi[5];
      }
    else{
      types = types + " & #pi^{-}";
      h1_p_stacked_npi[0]->Add(h1_pi_p_0_npi[0]);
      h1_p_stacked_npi[1]->Add(h1_pi_p_0_npi[1]);
      h1_p_stacked_npi[2]->Add(h1_pi_p_0_npi[2]);
      h1_p_stacked_npi[3]->Add(h1_pi_p_0_npi[3]);
      h1_p_stacked_npi[4]->Add(h1_pi_p_0_npi[4]);
      h1_p_stacked_npi[5]->Add(h1_pi_p_0_npi[5]);
      h1_p_stacked_1pi[0]->Add(h1_pi_p_0_1pi[0]);
      h1_p_stacked_1pi[1]->Add(h1_pi_p_0_1pi[1]);
      h1_p_stacked_1pi[2]->Add(h1_pi_p_0_1pi[2]);
      h1_p_stacked_1pi[3]->Add(h1_pi_p_0_1pi[3]);
      h1_p_stacked_1pi[4]->Add(h1_pi_p_0_1pi[4]);
      h1_p_stacked_1pi[5]->Add(h1_pi_p_0_1pi[5]);
    }
  }
  if (pi0){
    h1_pi_p_1_npi[0] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_1_type_1")); //QES
    h1_pi_p_1_npi[1] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_2_type_1")); //MEC
    h1_pi_p_1_npi[2] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_3_type_1")); //OtherRES
    h1_pi_p_1_npi[3] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_4_type_1")); //DIS
    h1_pi_p_1_npi[4] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_5_type_1")); //1232RES
    h1_pi_p_1_npi[5] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_0_type_1")); //Other
    h1_pi_p_1_1pi[0] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_1_type_1")); //QES
    h1_pi_p_1_1pi[1] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_2_type_1")); //MEC
    h1_pi_p_1_1pi[2] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_3_type_1")); //OtherRES
    h1_pi_p_1_1pi[3] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_4_type_1")); //DIS
    h1_pi_p_1_1pi[4] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_5_type_1")); //1232RES
    h1_pi_p_1_1pi[5] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_0_type_1")); //Other
    if (h1_p_stacked_npi[2] == nullptr){
      types = types + "#pi^{0}";
      h1_p_stacked_npi[0] = h1_pi_p_1_npi[0];
      h1_p_stacked_npi[1] = h1_pi_p_1_npi[1];
      h1_p_stacked_npi[2] = h1_pi_p_1_npi[2];
      h1_p_stacked_npi[3] = h1_pi_p_1_npi[3];
      h1_p_stacked_npi[4] = h1_pi_p_1_npi[4];
      h1_p_stacked_npi[5] = h1_pi_p_1_npi[5];
      h1_p_stacked_1pi[0] = h1_pi_p_1_1pi[0];
      h1_p_stacked_1pi[1] = h1_pi_p_1_1pi[1];
      h1_p_stacked_1pi[2] = h1_pi_p_1_1pi[2];
      h1_p_stacked_1pi[3] = h1_pi_p_1_1pi[3];
      h1_p_stacked_1pi[4] = h1_pi_p_1_1pi[4];
      h1_p_stacked_1pi[5] = h1_pi_p_1_1pi[5];
      }
    else{
      types = types + " & #pi^{0}";
      h1_p_stacked_npi[0]->Add(h1_pi_p_1_npi[0]);
      h1_p_stacked_npi[1]->Add(h1_pi_p_1_npi[1]);
      h1_p_stacked_npi[2]->Add(h1_pi_p_1_npi[2]);
      h1_p_stacked_npi[3]->Add(h1_pi_p_1_npi[3]);
      h1_p_stacked_npi[4]->Add(h1_pi_p_1_npi[4]);
      h1_p_stacked_npi[5]->Add(h1_pi_p_1_npi[5]);
      h1_p_stacked_1pi[0]->Add(h1_pi_p_1_1pi[0]);
      h1_p_stacked_1pi[1]->Add(h1_pi_p_1_1pi[1]);
      h1_p_stacked_1pi[2]->Add(h1_pi_p_1_1pi[2]);
      h1_p_stacked_1pi[3]->Add(h1_pi_p_1_1pi[3]);
      h1_p_stacked_1pi[4]->Add(h1_pi_p_1_1pi[4]);
      h1_p_stacked_1pi[5]->Add(h1_pi_p_1_1pi[5]);
    }
  }
  if (piplu){
    h1_pi_p_2_npi[0] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_1_type_2")); //QES
    h1_pi_p_2_npi[1] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_2_type_2")); //MEC
    h1_pi_p_2_npi[2] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_3_type_2")); //OtherRES
    h1_pi_p_2_npi[3] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_4_type_2")); //DIS
    h1_pi_p_2_npi[4] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_5_type_2")); //1232RES
    h1_pi_p_2_npi[5] = (TH1F*)input1->Get(TString::Format("h1_Npion_p_interaction_0_type_2")); //Other
    h1_pi_p_2_1pi[0] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_1_type_2")); //QES
    h1_pi_p_2_1pi[1] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_2_type_2")); //MEC
    h1_pi_p_2_1pi[2] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_3_type_2")); //OtherRES
    h1_pi_p_2_1pi[3] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_4_type_2")); //DIS
    h1_pi_p_2_1pi[4] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_5_type_2")); //1232RES
    h1_pi_p_2_1pi[5] = (TH1F*)input1->Get(TString::Format("h1_1pion_p_interaction_0_type_2")); //Other
    if (h1_p_stacked_npi[2] == nullptr){
      types = types + "#pi^{+}";
      h1_p_stacked_npi[0] = h1_pi_p_2_npi[0];
      h1_p_stacked_npi[1] = h1_pi_p_2_npi[1];
      h1_p_stacked_npi[2] = h1_pi_p_2_npi[2];
      h1_p_stacked_npi[3] = h1_pi_p_2_npi[3];
      h1_p_stacked_npi[4] = h1_pi_p_2_npi[4];
      h1_p_stacked_npi[5] = h1_pi_p_2_npi[5];
      h1_p_stacked_1pi[0] = h1_pi_p_2_1pi[0];
      h1_p_stacked_1pi[1] = h1_pi_p_2_1pi[1];
      h1_p_stacked_1pi[2] = h1_pi_p_2_1pi[2];
      h1_p_stacked_1pi[3] = h1_pi_p_2_1pi[3];
      h1_p_stacked_1pi[4] = h1_pi_p_2_1pi[4];
      h1_p_stacked_1pi[5] = h1_pi_p_2_1pi[5];
      }
    else{
      types = types + " & #pi^{+}";
      h1_p_stacked_npi[0]->Add(h1_pi_p_2_npi[0]);
      h1_p_stacked_npi[1]->Add(h1_pi_p_2_npi[1]);
      h1_p_stacked_npi[2]->Add(h1_pi_p_2_npi[2]);
      h1_p_stacked_npi[3]->Add(h1_pi_p_2_npi[3]);
      h1_p_stacked_npi[4]->Add(h1_pi_p_2_npi[4]);
      h1_p_stacked_npi[5]->Add(h1_pi_p_2_npi[5]);
      h1_p_stacked_1pi[0]->Add(h1_pi_p_2_1pi[0]);
      h1_p_stacked_1pi[1]->Add(h1_pi_p_2_1pi[1]);
      h1_p_stacked_1pi[2]->Add(h1_pi_p_2_1pi[2]);
      h1_p_stacked_1pi[3]->Add(h1_pi_p_2_1pi[3]);
      h1_p_stacked_1pi[4]->Add(h1_pi_p_2_1pi[4]);
      h1_p_stacked_1pi[5]->Add(h1_pi_p_2_1pi[5]);
    }
  }

  //Create the things to plot
  h1_toPlot_npi[0] = h1_p_stacked_npi[3]; //DIS

  h1_toPlot_npi[1] = h1_p_stacked_npi[4]; //1232RES

  h1_toPlot_npi[2] = (TH1F*)h1_p_stacked_npi[0]->Clone();//All other ones
  h1_toPlot_npi[2]->Add(h1_p_stacked_npi[1]);
  h1_toPlot_npi[2]->Add(h1_p_stacked_npi[2]);
  h1_toPlot_npi[2]->Add(h1_p_stacked_npi[5]);

  h1_toPlot_npi[3] = (TH1F*)h1_p_stacked_npi[0]->Clone(); // Total of everything
  h1_toPlot_npi[3]->Add(h1_p_stacked_npi[1]);
  h1_toPlot_npi[3]->Add(h1_p_stacked_npi[2]);
  h1_toPlot_npi[3]->Add(h1_p_stacked_npi[3]);
  h1_toPlot_npi[3]->Add(h1_p_stacked_npi[4]);
  h1_toPlot_npi[3]->Add(h1_p_stacked_npi[5]);

  h1_toPlot_1pi[0] = h1_p_stacked_1pi[3]; //DIS

  h1_toPlot_1pi[1] = h1_p_stacked_1pi[4]; //1232RES

  h1_toPlot_1pi[2] = (TH1F*)h1_p_stacked_1pi[0]->Clone();//All other ones
  h1_toPlot_1pi[2]->Add(h1_p_stacked_1pi[1]);
  h1_toPlot_1pi[2]->Add(h1_p_stacked_1pi[2]);
  h1_toPlot_1pi[2]->Add(h1_p_stacked_1pi[5]);

  h1_toPlot_1pi[3] = (TH1F*)h1_p_stacked_1pi[0]->Clone(); // Total of everything
  h1_toPlot_1pi[3]->Add(h1_p_stacked_1pi[1]);
  h1_toPlot_1pi[3]->Add(h1_p_stacked_1pi[2]);
  h1_toPlot_1pi[3]->Add(h1_p_stacked_1pi[3]);
  h1_toPlot_1pi[3]->Add(h1_p_stacked_1pi[4]);
  h1_toPlot_1pi[3]->Add(h1_p_stacked_1pi[5]);

  std::string particles = "";
  if (pimin){particles = particles + "Pimin";}
  if (pi0){particles = particles + "Pizero";}
  if (piplu){particles = particles + "Piplus";}

  // Makes the title look neater
  std::string accept = "";
  std::string a = "";
  std::string numProt = "";
  std::string energ = "";
  std::string nucleus = "";
  if (file_name.find("noAcc") != string::npos) {
    accept = accept + "Without Detector Acceptance";
    a = a + "noAcc";
  }
  else {
    accept = accept + "With Detector Acceptance";
    a = a + "Acc";
  }
  if (file_name.find("0P") != string::npos) {numProt = numProt + "0";}
  else if (file_name.find("-1P") != string::npos) {numProt = numProt + "n";}
  else if (file_name.find("2P") != string::npos) {numProt = numProt + "2";}
  else if (file_name.find("3P") != string::npos) {numProt = numProt + "3";}
  else if (file_name.find("1P") != string::npos) {numProt = numProt + "1";}
  if (file_name.find("1.161") != string::npos) {energ = energ + "1.161";}
  else if (file_name.find("2.261") != string::npos) {energ = energ + "2.261";}
  else if (file_name.find("4.461") != string::npos) {energ = energ + "4.461";}
  if (file_name.find("4He") != string::npos) {nucleus = "4He";}
  else if (file_name.find("56Fe") != string::npos) {nucleus = "56Fe";}
  else if (file_name.find("C12") != string::npos) {nucleus = "C12";}

  //Place Threshold position
  float line_pos;
  if (pi0) {line_pos = 0.3;}
  else {line_pos = 0.15;}

  //------------------------------------------------------------ N pion events plotting -------------------------------


  TCanvas* c1;
  c1 = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

  // initialize things such as color options that will be used in plotting the histograms below
  const int color_options[4] = {kBlue, kAzure+10, kGreen-3, kViolet};
  const char * draw_options[4] = {"HIST","HIST SAME","HIST SAME", "HIST SAME"};
  const char * draw_options2[4] = {"E SAME", "E SAME","E SAME", "E SAME"};

  double hist_max_n = 0.;
  double x_max_n[4] = {0.,0.,0.,0.};
  double maximumx_n = 0.;
  // format the histograms and then make them
  for( int i = 0; i < 4; i++){
  // formatting statements
     if (pi0 != 1){h1_toPlot_npi[i]->Rebin(2);} // Halve the number of bins when not pi0. Comment out (or delete) if unnecesary
     else {h1_toPlot_npi[i]->Rebin(4);} // Decrease the number of bins for pi0 (much fewer events than the others). Comment out (or delete) if unnecesary
     h1_toPlot_npi[i]->Sumw2();
  //       h1_el_energyentum_[i]->Scale( 1./h1_el_energyentum_[i]->Integral(),"width"); // normalize the histogram to 1
     h1_toPlot_npi[i]->SetStats(0); // get rid of the stats box that usually appears at the top right of plots
     h1_toPlot_npi[i]->SetLineColor( color_options[i]);
     h1_toPlot_npi[i]->GetXaxis()->SetTitle("Pion Momentum (GeV/c)");
     if (pi0 != 1){h1_toPlot_npi[i]->GetYaxis()->SetTitle("Pions per 20 MeV");}
     else {h1_toPlot_npi[i]->GetYaxis()->SetTitle("Pions per 40 MeV");}
     //h1_toPlot_npi[i]->GetXaxis()->SetRangeUser( 0, 1);
      //h1_toPlot_npi[i]->GetYaxis()->SetRangeUser( 0, 10);
     h1_toPlot_npi[i]->GetXaxis()->CenterTitle( true);
     h1_toPlot_npi[i]->GetYaxis()->CenterTitle( true);
  //        h1_el_momentum[i]->GetXaxis()->SetTitleSize( 0.06);
  //        h1_el_momentum[i]->GetYaxis()->SetTitleSize( 0.06);

     // Take your pick for the title -- top is for all information and bottom is for basic info
     //h1_toPlot_npi[i]->SetTitle(TString(energ)+" GeV " + TString(nucleus) + " (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));
     h1_toPlot_npi[i]->SetTitle("Momentum Distribution of (1enpi) Pions "+TString(accept) + " for "+TString(types));

      // take the maximum of the four interactions to use later in plotting
      if (h1_toPlot_npi[i]->GetMaximum() > hist_max_n) { hist_max_n = h1_toPlot_npi[i]->GetMaximum(); }

      //Find the maximum x range for each graph
      for (int j = h1_toPlot_npi[i]->GetNbinsX(); j-- > 0; ){ // Don't ask me why, but this is the best way to count backwards in a for-loop
        if (h1_toPlot_npi[i]->GetBinContent(j-1) >= 1. ) {break;}
        x_max_n[i] = h1_toPlot_npi[i]->GetBinCenter(j-1);
      }

      //Find maximum x-value in array
      for (int m = 0; m < 4; m++){
        if (x_max_n[m] > maximumx_n) {maximumx_n = x_max_n[m];}
      }

  }

  // draw the histograms
  for( int i = 0; i < 4; i++){
     h1_toPlot_npi[i]->GetYaxis()->SetRangeUser( 0, hist_max_n*1.1); // smithja: hist_max mult. by 1.1 to give the graph some head room
     h1_toPlot_npi[i]->GetXaxis()->SetRangeUser(0,TMath::Ceil(maximumx_n));
     h1_toPlot_npi[i]->Draw( draw_options[i]);
     h1_toPlot_npi[i]->Draw( draw_options2[i]);
  }

  // Draw a line at the threshold
  TLine* line_npi = new TLine(line_pos, 0, line_pos, 1.05*hist_max_n);
    line_npi->SetLineColor(kBlack);
    line_npi->SetLineStyle(1);
    line_npi->SetLineWidth(1);
    line_npi->Draw("SAME");

  // crop the margins of the canvas
  c1->SetLeftMargin( 0.14);
  c1->SetBottomMargin( 0.14);
  c1->Update();

  // smithja: this section of code integrates the distributions and prints
  //          it to the screen
  int num_DIS_events_above_npi;
  int num_1232_events_above_npi;
  int num_Other_events_above_npi;
  int num_Total_events_above_npi;

  int num_DIS_events_below_npi;
  int num_1232_events_below_npi;
  int num_Other_events_below_npi;
  int num_Total_events_below_npi;

  if (pi0){
  num_DIS_events_above_npi = h1_toPlot_npi[0]->Integral(h1_toPlot_npi[0]->FindBin(0.3),4);
  num_1232_events_above_npi = h1_toPlot_npi[1]->Integral(h1_toPlot_npi[1]->FindBin(0.3),4);
  num_Other_events_above_npi = h1_toPlot_npi[2]->Integral(h1_toPlot_npi[2]->FindBin(0.3),4);
  num_Total_events_above_npi = h1_toPlot_npi[3]->Integral(h1_toPlot_npi[3]->FindBin(0.3),4);
}

  else{
  num_DIS_events_above_npi = h1_toPlot_npi[0]->Integral(h1_toPlot_npi[0]->FindBin(0.15),4);
  num_1232_events_above_npi = h1_toPlot_npi[1]->Integral(h1_toPlot_npi[1]->FindBin(0.15),4);
  num_Other_events_above_npi = h1_toPlot_npi[2]->Integral(h1_toPlot_npi[2]->FindBin(0.15),4);
  num_Total_events_above_npi = h1_toPlot_npi[3]->Integral(h1_toPlot_npi[3]->FindBin(0.15),4);
  }

  if (pi0){
  num_DIS_events_below_npi = h1_toPlot_npi[0]->Integral(0,h1_toPlot_npi[0]->FindBin(0.3));
  num_1232_events_below_npi = h1_toPlot_npi[1]->Integral(0,h1_toPlot_npi[1]->FindBin(0.3));
  num_Other_events_below_npi = h1_toPlot_npi[2]->Integral(0,h1_toPlot_npi[2]->FindBin(0.3));
  num_Total_events_below_npi = h1_toPlot_npi[3]->Integral(0,h1_toPlot_npi[3]->FindBin(0.3));
}

  else{
  num_DIS_events_below_npi = h1_toPlot_npi[0]->Integral(0,h1_toPlot_npi[0]->FindBin(0.15));
  num_1232_events_below_npi = h1_toPlot_npi[1]->Integral(0,h1_toPlot_npi[1]->FindBin(0.15));
  num_Other_events_below_npi = h1_toPlot_npi[2]->Integral(0,h1_toPlot_npi[2]->FindBin(0.15));
  num_Total_events_below_npi = h1_toPlot_npi[3]->Integral(0,h1_toPlot_npi[3]->FindBin(0.15));
  }
  //int num_events = num_QE_events + num_MEC_events + num_RES_events + num_DIS_events;

  // create and draw a legend for our histograms
  TLegend *legend_npi_above = new TLegend( 0.62, 0.59, 1.07, 0.74);
  legend_npi_above->SetHeader("# of Pions Above Threshold", "l");
  legend_npi_above->AddEntry( h1_toPlot_npi[0],Form("DIS (%d)",num_DIS_events_above_npi),"l");
  legend_npi_above->AddEntry( h1_toPlot_npi[1],Form("#Delta(1232) (%d)",num_1232_events_above_npi),"l");
  legend_npi_above->AddEntry( h1_toPlot_npi[2],Form("Other (%d)",num_Other_events_above_npi),"l");
  legend_npi_above->AddEntry( h1_toPlot_npi[3],Form("Total (%d)",num_Total_events_above_npi),"l");
  legend_npi_above->SetBorderSize( 0);
  legend_npi_above->SetFillStyle( 0);
  legend_npi_above->Draw();

  // create and draw a legend for our histograms
  TLegend *legend_npi_below = new TLegend( 0.62, 0.75, 1.07, 0.9);
  legend_npi_below->SetHeader("# of Pions Below Threshold", "l");
  legend_npi_below->AddEntry( h1_toPlot_npi[0],Form("DIS (%d)",num_DIS_events_below_npi),"l");
  legend_npi_below->AddEntry( h1_toPlot_npi[1],Form("#Delta(1232) (%d)",num_1232_events_below_npi),"l");
  legend_npi_below->AddEntry( h1_toPlot_npi[2],Form("Other (%d)",num_Other_events_below_npi),"l");
  legend_npi_below->AddEntry( h1_toPlot_npi[3],Form("Total (%d)",num_Total_events_below_npi),"l");
  legend_npi_below->SetBorderSize( 0);
  legend_npi_below->SetFillStyle( 0);
  legend_npi_below->Draw();

  // save this histogram as a PDF file
  std::string save_file_path_npi ("invalid");
  save_file_path_npi = "../output/h1_pion_momcount/"+TString(particles)+"momentum_NpiEvents"+TString(a)+".pdf";

  c1->SaveAs( TString::Format("%s", save_file_path_npi.c_str()));

  //--------------------------------------------------- 1 pion events plotting ---------------------------------------------

  TCanvas* c2;
  c2 = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

  double hist_max_1 = 0.;
  double x_max_1[4] = {0.,0.,0.,0.};
  double maximumx_1 = 0.;
  // format the histograms and then make them
  for( int i = 0; i < 4; i++){
  // formatting statements
     if (pi0 != 1){h1_toPlot_1pi[i]->Rebin(2);} // Halve the number of bins when not pi0. Comment out (or delete) if unnecesary
     else {h1_toPlot_1pi[i]->Rebin(4);} // Decrease the number of bins for pi0 (much fewer events than the others). Comment out (or delete) if unnecesary
     h1_toPlot_1pi[i]->Sumw2();
  //       h1_el_energyentum_[i]->Scale( 1./h1_el_energyentum_[i]->Integral(),"width"); // normalize the histogram to 1
     h1_toPlot_1pi[i]->SetStats(0); // get rid of the stats box that usually appears at the top right of plots
     h1_toPlot_1pi[i]->SetLineColor( color_options[i]);
     h1_toPlot_1pi[i]->GetXaxis()->SetTitle("Pion Momentum (GeV/c)");
     if (pi0 != 1){h1_toPlot_1pi[i]->GetYaxis()->SetTitle("Pions per 20 MeV");}
     else {h1_toPlot_1pi[i]->GetYaxis()->SetTitle("Pions per 40 MeV");}
     //h1_toPlot_1pi[i]->GetXaxis()->SetRangeUser( 0, 1);
      //h1_toPlot_1pi[i]->GetYaxis()->SetRangeUser( 0, 10);
     h1_toPlot_1pi[i]->GetXaxis()->CenterTitle( true);
     h1_toPlot_1pi[i]->GetYaxis()->CenterTitle( true);
  //        h1_el_momentum[i]->GetXaxis()->SetTitleSize( 0.06);
  //        h1_el_momentum[i]->GetYaxis()->SetTitleSize( 0.06);

     // Take your pick for the title -- top is for all information and bottom is for basic info
     //h1_toPlot_1pi[i]->SetTitle(TString(energ)+" GeV " + TString(nucleus) + " (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));
     h1_toPlot_1pi[i]->SetTitle("Momentum Distribution of (1e1pi) Pions "+TString(accept) + " for "+TString(types));

      // take the maximum of the four interactions to use later in plotting
      if (h1_toPlot_1pi[i]->GetMaximum() > hist_max_1) { hist_max_1 = h1_toPlot_1pi[i]->GetMaximum(); }

      //Find the maximum x range for each graph
      for (int j = h1_toPlot_1pi[i]->GetNbinsX(); j-- > 0; ){ // Don't ask me why, but this is the best way to count backwards in a for-loop
        if (h1_toPlot_1pi[i]->GetBinContent(j-1) >= 1. ) {break;}
        x_max_1[i] = h1_toPlot_1pi[i]->GetBinCenter(j-1);
      }

      //Find maximum x-value in array
      for (int m = 0; m < 4; m++){
        if (x_max_1[m] > maximumx_1) {maximumx_1 = x_max_1[m];}
      }

  }

  // draw the histograms
  for( int i = 0; i < 4; i++){
     h1_toPlot_1pi[i]->GetYaxis()->SetRangeUser( 0, hist_max_1*1.1); // smithja: hist_max mult. by 1.1 to give the graph some head room
     h1_toPlot_1pi[i]->GetXaxis()->SetRangeUser(0,TMath::Ceil(maximumx_1));
     h1_toPlot_1pi[i]->Draw( draw_options[i]);
     h1_toPlot_1pi[i]->Draw( draw_options2[i]);
  }

  // Draw a line at the threshold
  TLine* line_1pi = new TLine(line_pos, 0, line_pos, 1.05*hist_max_1);
    line_1pi->SetLineColor(kBlack);
    line_1pi->SetLineStyle(1);
    line_1pi->SetLineWidth(1);
    line_1pi->Draw("SAME");

  // crop the margins of the canvas
  c2->SetLeftMargin( 0.14);
  c2->SetBottomMargin( 0.14);
  c2->Update();

  // smithja: this section of code integrates the distributions and prints
  //          it to the screen
  int num_DIS_events_above_1pi;
  int num_1232_events_above_1pi;
  int num_Other_events_above_1pi;
  int num_Total_events_above_1pi;

  int num_DIS_events_below_1pi;
  int num_1232_events_below_1pi;
  int num_Other_events_below_1pi;
  int num_Total_events_below_1pi;

  if (pi0){
  num_DIS_events_above_1pi = h1_toPlot_1pi[0]->Integral(h1_toPlot_1pi[0]->FindBin(0.3)+1,4);
  num_1232_events_above_1pi = h1_toPlot_1pi[1]->Integral(h1_toPlot_1pi[1]->FindBin(0.3)+1,4);
  num_Other_events_above_1pi = h1_toPlot_1pi[2]->Integral(h1_toPlot_1pi[2]->FindBin(0.3)+1,4);
  num_Total_events_above_1pi = h1_toPlot_1pi[3]->Integral(h1_toPlot_1pi[3]->FindBin(0.3)+1,4);
}

  else{
  num_DIS_events_above_1pi = h1_toPlot_1pi[0]->Integral(h1_toPlot_1pi[0]->FindBin(0.15)+1,4);
  num_1232_events_above_1pi = h1_toPlot_1pi[1]->Integral(h1_toPlot_1pi[1]->FindBin(0.15)+1,4);
  num_Other_events_above_1pi = h1_toPlot_1pi[2]->Integral(h1_toPlot_1pi[2]->FindBin(0.15)+1,4);
  num_Total_events_above_1pi = h1_toPlot_1pi[3]->Integral(h1_toPlot_1pi[3]->FindBin(0.15)+1,4);
  }

  if (pi0){
  num_DIS_events_below_1pi = h1_toPlot_1pi[0]->Integral(0,h1_toPlot_1pi[0]->FindBin(0.3));
  num_1232_events_below_1pi = h1_toPlot_1pi[1]->Integral(0,h1_toPlot_1pi[0]->FindBin(0.3));
  num_Other_events_below_1pi = h1_toPlot_1pi[2]->Integral(0,h1_toPlot_1pi[0]->FindBin(0.3));
  num_Total_events_below_1pi = h1_toPlot_1pi[3]->Integral(0,h1_toPlot_1pi[0]->FindBin(0.3));
}

  else{
  num_DIS_events_below_1pi = h1_toPlot_1pi[0]->Integral(0,h1_toPlot_1pi[0]->FindBin(0.15));
  num_1232_events_below_1pi = h1_toPlot_1pi[1]->Integral(0,h1_toPlot_1pi[1]->FindBin(0.15));
  num_Other_events_below_1pi = h1_toPlot_1pi[2]->Integral(0,h1_toPlot_1pi[2]->FindBin(0.15));
  num_Total_events_below_1pi = h1_toPlot_1pi[3]->Integral(0,h1_toPlot_1pi[3]->FindBin(0.15));
  }
  //int num_events = num_QE_events + num_MEC_events + num_RES_events + num_DIS_events;

  // create and draw a legend for our histograms
  TLegend *legend_1pi_above = new TLegend( 0.62, 0.75, 1.07, 0.9);
  legend_1pi_above->SetHeader("# of Pions Above Threshold", "l");
  legend_1pi_above->AddEntry( h1_toPlot_1pi[0],Form("DIS (%d)",num_DIS_events_above_1pi),"l");
  legend_1pi_above->AddEntry( h1_toPlot_1pi[1],Form("#Delta(1232) (%d)",num_1232_events_above_1pi),"l");
  legend_1pi_above->AddEntry( h1_toPlot_1pi[2],Form("Other (%d)",num_Other_events_above_1pi),"l");
  legend_1pi_above->AddEntry( h1_toPlot_1pi[3],Form("Total (%d)",num_Total_events_above_1pi),"l");
  legend_1pi_above->SetBorderSize( 0);
  legend_1pi_above->SetFillStyle( 0);
  legend_1pi_above->Draw();

  // create and draw a legend for our histograms
  TLegend *legend_1pi_below = new TLegend( 0.62, 0.59, 1.07, 0.74);
  legend_1pi_below->SetHeader("# of Pions Below Threshold", "l");
  legend_1pi_below->AddEntry( h1_toPlot_1pi[0],Form("DIS (%d)",num_DIS_events_below_1pi),"l");
  legend_1pi_below->AddEntry( h1_toPlot_1pi[1],Form("#Delta(1232) (%d)",num_1232_events_below_1pi),"l");
  legend_1pi_below->AddEntry( h1_toPlot_1pi[2],Form("Other (%d)",num_Other_events_below_1pi),"l");
  legend_1pi_below->AddEntry( h1_toPlot_1pi[3],Form("Total (%d)",num_Total_events_below_1pi),"l");
  legend_1pi_below->SetBorderSize( 0);
  legend_1pi_below->SetFillStyle( 0);
  legend_1pi_below->Draw();

  // save this histogram as a PDF file
  std::string save_file_path_1pi ("invalid");
  save_file_path_1pi = "../output/h1_pion_momcount/"+TString(particles)+"momentum_1piEvents"+TString(a)+".pdf";

  c2->SaveAs( TString::Format("%s", save_file_path_1pi.c_str()));
}
