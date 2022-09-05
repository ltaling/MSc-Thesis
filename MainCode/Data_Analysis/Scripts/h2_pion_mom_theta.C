// This simple macro creates the 2D histograms for the momentum and theta of the signal inducing pions
// Specifically, 1e0p1pi events.
// This is heavily based on, and in some places copied from, the code by Noah and Graham that also produce plots.
// My sincerest apologies for the dramatic decrease in quality.

// These include statements were here when I received Graham's code. Delete/Modify them at
// your own peril.
#include <TFile.h>
#include <TH2D.h>
#include <TF2.h>
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

void h2_pion_mom_theta (std::string file_name, bool pimin, bool pi0, bool piplu){

  // declare variables that will be used in the plotting/formatting of histograms and file names
  std::string info ("");
  std::string cuts ("");

  // import the ROOT data file from the output of doing ./genie_analysis in the e4nu directory
  TFile *input1 = TFile::Open( TString( file_name));

  // Grab the TList in the ROOT File for info on the run
  TList *Run_Info = (TList*)input1->Get("Run_Info");

  std::vector<std::string> parse = parse_file(Run_Info);
  std:string types = "";

  info = parse[0];
  cuts = parse[1];

  // Histogram initialisation
  TH2F* h2_stacked[5]; // Bins of all pion types added together
  TH2F* h2_pi_mom_theta_0[5]; //pi-
  TH2F* h2_pi_mom_theta_1[5]; //pi0
  TH2F* h2_pi_mom_theta_2[5]; //pi+
  // Drawing histograms
  if (pimin){
    h2_pi_mom_theta_0[0] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_1_type_0")); //QES
    h2_pi_mom_theta_0[1] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_2_type_0")); //MEC
    h2_pi_mom_theta_0[2] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_3_type_0")); //RES
    h2_pi_mom_theta_0[3] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_4_type_0")); //DIS
    h2_pi_mom_theta_0[4] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_0_type_0")); //Other
    if (h2_stacked[2] == nullptr){
      types = types + "min";
      h2_stacked[0] = h2_pi_mom_theta_0[0];
      h2_stacked[1] = h2_pi_mom_theta_0[1];
      h2_stacked[2] = h2_pi_mom_theta_0[2];
      h2_stacked[3] = h2_pi_mom_theta_0[3];
      h2_stacked[4] = h2_pi_mom_theta_0[4];
      }
    else{
      types = types + " & min ";
      h2_stacked[0]->Add(h2_pi_mom_theta_0[0]);
      h2_stacked[1]->Add(h2_pi_mom_theta_0[1]);
      h2_stacked[2]->Add(h2_pi_mom_theta_0[2]);
      h2_stacked[3]->Add(h2_pi_mom_theta_0[3]);
      h2_stacked[4]->Add(h2_pi_mom_theta_0[4]);
    }
  }
  if (pi0){
    h2_pi_mom_theta_1[0] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_1_type_1")); //QES
    h2_pi_mom_theta_1[1] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_2_type_1")); //MEC
    h2_pi_mom_theta_1[2] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_3_type_1")); //RES
    h2_pi_mom_theta_1[3] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_4_type_1")); //DIS
    h2_pi_mom_theta_1[4] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_0_type_1")); //Other
    if (h2_stacked[2] == nullptr){
      types = types + "0";
      h2_stacked[0] = h2_pi_mom_theta_1[0];
      h2_stacked[1] = h2_pi_mom_theta_1[1];
      h2_stacked[2] = h2_pi_mom_theta_1[2];
      h2_stacked[3] = h2_pi_mom_theta_1[3];
      h2_stacked[4] = h2_pi_mom_theta_1[4];
      }
    else{
      types = types + " & 0";
      h2_stacked[0]->Add(h2_pi_mom_theta_1[0]);
      h2_stacked[1]->Add(h2_pi_mom_theta_1[1]);
      h2_stacked[2]->Add(h2_pi_mom_theta_1[2]);
      h2_stacked[3]->Add(h2_pi_mom_theta_1[3]);
      h2_stacked[4]->Add(h2_pi_mom_theta_1[4]);
    }
  }
  if (piplu){
    h2_pi_mom_theta_2[0] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_1_type_2")); //QES
    h2_pi_mom_theta_2[1] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_2_type_2")); //MEC
    h2_pi_mom_theta_2[2] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_3_type_2")); //RES
    h2_pi_mom_theta_2[3] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_4_type_2")); //DIS
    h2_pi_mom_theta_2[4] = (TH2F*)input1->Get(TString::Format("h2_pion_mom_theta_interaction_0_type_2")); //Other
    if (h2_stacked[2] == nullptr){
      types = types + "plus";
      h2_stacked[0] = h2_pi_mom_theta_2[0];
      h2_stacked[1] = h2_pi_mom_theta_2[1];
      h2_stacked[2] = h2_pi_mom_theta_2[2];
      h2_stacked[3] = h2_pi_mom_theta_2[3];
      h2_stacked[4] = h2_pi_mom_theta_2[4];
      }
    else{
      types = types + " & plus";
      h2_stacked[0]->Add(h2_pi_mom_theta_2[0]);
      h2_stacked[1]->Add(h2_pi_mom_theta_2[1]);
      h2_stacked[2]->Add(h2_pi_mom_theta_2[2]);
      h2_stacked[3]->Add(h2_pi_mom_theta_2[3]);
      h2_stacked[4]->Add(h2_pi_mom_theta_2[4]);
    }
  }

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
  else if (file_name.find("1P") != string::npos) {numProt = numProt + "1";}
  else if (file_name.find("2P") != string::npos) {numProt = numProt + "2";}
  else if (file_name.find("3P") != string::npos) {numProt = numProt + "3";}
  else if (file_name.find("nP") != string::npos) {numProt = numProt + "n";}
  if (file_name.find("1.161") != string::npos) {energ = energ + "1.161";}
  else if (file_name.find("2.261") != string::npos) {energ = energ + "2.261";}
  else if (file_name.find("4.461") != string::npos) {energ = energ + "4.461";}
  if (file_name.find("4He") != string::npos) {nucleus = "4He";}
  else if (file_name.find("56Fe") != string::npos) {nucleus = "56Fe";}
  else if (file_name.find("C12") != string::npos) {nucleus = "C12";}

  TCanvas* c;
  c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

  // initialize things such as color options that will be used in plotting the histograms below
  const int color_options[5] = {kBlue, kAzure+10, kGreen-3, kViolet, kBlue+2};
  const char * draw_options[5] = {"HIST","HIST SAME","HIST SAME", "HIST SAME","HIST SAME"};
  //const char * draw_options2[5] = {"E SAME", "E SAME","E SAME", "E SAME","E SAME"};

  double hist_max = 0.;
  double x_max[5] = {0.,0.,0.,0.,0.};
  double maximumx = 0.;
  // format the histograms and then make them
  for( int i = 0; i < 5; i++){
// formatting statements
      //h2_theta_[i]->Sumw2();
//        h1_el_momentum[i]->Scale( 1./h1_el_momentum[i]->Integral(),"width"); // normalize the histogram to 1
      h2_stacked[i]->SetStats(0); // get rid of the stats box that usually appears at the top right of plots
      h2_stacked[i]->SetMarkerColor( color_options[i]);
      h2_stacked[i]->SetFillColor( color_options[i]);
      h2_stacked[i]->GetXaxis()->SetTitle("Pion Momentum (GeV/c)");
      h2_stacked[i]->GetYaxis()->SetTitle("Theta (degrees)");
      //h2_stacked[i]->GetXaxis()->SetRangeUser( 0, 1);
      //h2_stacked[i]->GetYaxis()->SetRangeUser( 0, 10);
      h2_stacked[i]->GetXaxis()->CenterTitle( true);
      h2_stacked[i]->GetYaxis()->CenterTitle( true);
//        h1_el_momentum[i]->GetXaxis()->SetTitleSize( 0.06);
//        h1_el_momentum[i]->GetYaxis()->SetTitleSize( 0.06);


      h2_stacked[i]->SetTitle(TString(energ)+" Gev " + TString(nucleus) + " (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));

      //Find the maximum x range for each graph
      for (int j = h2_stacked[i]->GetNbinsX(); j-- > 0; ){ // Don't ask me why, but this is the best way to count backwards in a for-loop
        double yvals = 0;
        for (int k = h2_stacked[i]->GetNbinsY(); k-- > 0; ){
          yvals = yvals + h2_stacked[i]->GetBinContent(j-1,k);
        }
        if (yvals >= 1. ) {break;}
        x_max[i] = h2_stacked[i]->GetXaxis()->GetBinCenter(j-1);
      }

      //Find maximum x-value in array
      for (int m = 0; m < 5; m++){
        if (x_max[m] > maximumx) {maximumx = x_max[m];}
      }

      // take the maximum of the four interactions to use later in plotting
      //if (h2_stacked[i]->GetMaximum() > hist_max) { hist_max = h2_stacked[i]->GetMaximum(); }
  }

  // draw the histograms
  for( int i = 0; i < 5; i++){
//h2_stacked[i]->GetYaxis()->SetRangeUser( 0, hist_max*1.1); // smithja: hist_max mult. by 1.1 to give the graph some head room
      h2_stacked[i]->Draw( draw_options[i]);
      h2_stacked[i]->GetXaxis()->SetRangeUser(0,TMath::Ceil(maximumx));
      //h2_stacked[i]->Draw( draw_options2[i]);
  }

  // crop the margins of the canvas
  c->SetLeftMargin( 0.14);
  c->SetBottomMargin( 0.14);
  c->Update();

  // smithja: this section of code integrates the distributions and prints
  //          it to the screen
  int num_QE_events =h2_stacked[0]->Integral();
  int num_MEC_events =h2_stacked[1]->Integral();
  int num_RES_events =h2_stacked[2]->Integral();
  int num_DIS_events =h2_stacked[3]->Integral();
  int num_events = num_QE_events + num_MEC_events + num_RES_events + num_DIS_events;

  // create and draw a legend for our histograms
  TLegend *legend = new TLegend( 0.65, 0.75, 1.1, 0.9);
  legend->AddEntry( h2_stacked[0],Form("QE (%d)",num_QE_events),"f");
  legend->AddEntry( h2_stacked[1],Form("MEC (%d)",num_MEC_events),"f");
  legend->AddEntry( h2_stacked[2],Form("RES (%d)",num_RES_events),"f");
  legend->AddEntry( h2_stacked[3],Form("DIS (%d)",num_DIS_events),"f");
  legend->SetBorderSize( 0);
  legend->SetFillStyle( 0);
  legend->Draw();

  // save this histogram as a PDF file
  std::string save_file_path ("invalid");
  save_file_path = "../output/h2_pion_mom_theta/"+TString(energ)+TString(nucleus)+TString(cuts)+TString(particles)+TString(numProt)+"P"+TString(a)+".pdf";

  c->SaveAs( TString::Format("%s", save_file_path.c_str()));

/*
  cout << "// --- PROGRAM OUTPUT ------------------------------------------------------" << endl;
  cout << "There are " << num_QE_events << " many QE events in this histogram." << endl;
  cout << "There are " << num_MEC_events << " many MEC events in this histogram." << endl;
  cout << "There are " << num_RES_events << " many RES events in this histogram." << endl;
  cout << "There are " << num_DIS_events << " many DIS events in this histogram." << endl;
  cout << "There are " << num_events << " many total events in the h2_el_mom_insSeparated histogram." << endl;
  cout << "// -------------------------------------------------------------------------" << endl << endl; // separate program output
  */
}
