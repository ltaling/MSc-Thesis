// This macro draws the 1 dimensional histogram of outgoing electron angle (theta) with respect to the beam

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

void h1_pion_electrontheta (std::string file_name, bool pimin, bool pi0, bool piplu){
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
  TH1F* h1_stacked[6]; // Bins of all pion types added together
  TH1F* h1_pi_electrontheta_0[6]; //pi-
  TH1F* h1_pi_electrontheta_1[6]; //pi0
  TH1F* h1_pi_electrontheta_2[6]; //pi+
  TH1F* h1_toPlot[4]; //things to actally plot
  // Drawing histograms
  if (pimin){
    h1_pi_electrontheta_0[0] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_1_type_0")); //QES
    h1_pi_electrontheta_0[1] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_2_type_0")); //MEC
    h1_pi_electrontheta_0[2] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_3_type_0")); //OtherRES
    h1_pi_electrontheta_0[3] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_4_type_0")); //DIS
    h1_pi_electrontheta_0[4] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_5_type_0")); //1232RES
    h1_pi_electrontheta_0[5] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_0_type_0")); //Other
    if (h1_stacked[2] == nullptr){
      types = types + "#pi^{-}";
      h1_stacked[0] = h1_pi_electrontheta_0[0];
      h1_stacked[1] = h1_pi_electrontheta_0[1];
      h1_stacked[2] = h1_pi_electrontheta_0[2];
      h1_stacked[3] = h1_pi_electrontheta_0[3];
      h1_stacked[4] = h1_pi_electrontheta_0[4];
      h1_stacked[5] = h1_pi_electrontheta_0[5];
      }
    else{
      types = types + " & #pi^{-} ";
      h1_stacked[0]->Add(h1_pi_electrontheta_0[0]);
      h1_stacked[1]->Add(h1_pi_electrontheta_0[1]);
      h1_stacked[2]->Add(h1_pi_electrontheta_0[2]);
      h1_stacked[3]->Add(h1_pi_electrontheta_0[3]);
      h1_stacked[4]->Add(h1_pi_electrontheta_0[4]);
      h1_stacked[5]->Add(h1_pi_electrontheta_0[5]);
    }
  }
  if (pi0){
    h1_pi_electrontheta_1[0] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_1_type_1")); //QES
    h1_pi_electrontheta_1[1] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_2_type_1")); //MEC
    h1_pi_electrontheta_1[2] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_3_type_1")); //OtherRES
    h1_pi_electrontheta_1[3] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_4_type_1")); //DIS
    h1_pi_electrontheta_1[4] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_5_type_1")); //1232RES
    h1_pi_electrontheta_1[5] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_0_type_1")); //Other
    if (h1_stacked[2] == nullptr){
      types = types + "#pi^{0}";
      h1_stacked[0] = h1_pi_electrontheta_1[0];
      h1_stacked[1] = h1_pi_electrontheta_1[1];
      h1_stacked[2] = h1_pi_electrontheta_1[2];
      h1_stacked[3] = h1_pi_electrontheta_1[3];
      h1_stacked[4] = h1_pi_electrontheta_1[4];
      h1_stacked[5] = h1_pi_electrontheta_1[5];
      }
    else{
      types = types + " & #pi^{0}";
      h1_stacked[0]->Add(h1_pi_electrontheta_1[0]);
      h1_stacked[1]->Add(h1_pi_electrontheta_1[1]);
      h1_stacked[2]->Add(h1_pi_electrontheta_1[2]);
      h1_stacked[3]->Add(h1_pi_electrontheta_1[3]);
      h1_stacked[4]->Add(h1_pi_electrontheta_1[4]);
      h1_stacked[5]->Add(h1_pi_electrontheta_1[5]);
    }
  }
  if (piplu){
    h1_pi_electrontheta_2[0] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_1_type_2")); //QES
    h1_pi_electrontheta_2[1] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_2_type_2")); //MEC
    h1_pi_electrontheta_2[2] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_3_type_2")); //OtherRES
    h1_pi_electrontheta_2[3] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_4_type_2")); //DIS
    h1_pi_electrontheta_2[4] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_5_type_2")); //1232RES
    h1_pi_electrontheta_2[5] = (TH1F*)input1->Get(TString::Format("h1_pion_electrontheta_interaction_0_type_2")); //Other
    if (h1_stacked[2] == nullptr){
      types = types + "#pi^{+}";
      h1_stacked[0] = h1_pi_electrontheta_2[0];
      h1_stacked[1] = h1_pi_electrontheta_2[1];
      h1_stacked[2] = h1_pi_electrontheta_2[2];
      h1_stacked[3] = h1_pi_electrontheta_2[3];
      h1_stacked[4] = h1_pi_electrontheta_2[4];
      h1_stacked[5] = h1_pi_electrontheta_2[5];
      }
    else{
      types = types + " & #pi^{+}";
      h1_stacked[0]->Add(h1_pi_electrontheta_2[0]);
      h1_stacked[1]->Add(h1_pi_electrontheta_2[1]);
      h1_stacked[2]->Add(h1_pi_electrontheta_2[2]);
      h1_stacked[3]->Add(h1_pi_electrontheta_2[3]);
      h1_stacked[4]->Add(h1_pi_electrontheta_2[4]);
      h1_stacked[5]->Add(h1_pi_electrontheta_2[5]);
    }
  }

  //Create the things to plot
  h1_toPlot[0] = h1_stacked[3]; //DIS

  h1_toPlot[1] = h1_stacked[4]; //1232RES

  h1_toPlot[2] = (TH1F*)h1_stacked[0]->Clone();//All other ones
  h1_toPlot[2]->Add(h1_stacked[1]);
  h1_toPlot[2]->Add(h1_stacked[2]);
  h1_toPlot[2]->Add(h1_stacked[5]);

  h1_toPlot[3] = (TH1F*)h1_stacked[0]->Clone(); // Total of everything
  h1_toPlot[3]->Add(h1_stacked[1]);
  h1_toPlot[3]->Add(h1_stacked[2]);
  h1_toPlot[3]->Add(h1_stacked[3]);
  h1_toPlot[3]->Add(h1_stacked[4]);
  h1_toPlot[3]->Add(h1_stacked[5]);


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

  TCanvas* c;
  c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

  // initialize things such as color options that will be used in plotting the histograms below
  const int color_options[5] = {kBlue, kAzure+10, kGreen-3, kViolet};
  const char * draw_options[5] = {"HIST","HIST SAME","HIST SAME", "HIST SAME"};
  const char * draw_options2[5] = {"E SAME", "E SAME","E SAME", "E SAME"};

  double hist_max = 0.;
  // format the histograms and then make them
  for( int i = 0; i < 4; i++){
  // formatting statements
      if (pi0 != 1){h1_toPlot[i]->Rebin(2);} // Halve the number of bins when not pi0. Comment out (or delete) if unnecesary
      else {h1_toPlot[i]->Rebin(4);} // Decrease the number of bins for pi0 (much fewer events than the others). Comment out (or delete) if unnecesary
      h1_toPlot[i]->Sumw2();
  //       h1_el_phientum[i]->Scale( 1./h1_el_momentum[i]->Integral(),"width"); // normalize the histogram to 1
      h1_toPlot[i]->SetStats(0); // get rid of the stats box that usually appears at the top right of plots
      h1_toPlot[i]->SetLineColor( color_options[i]);
      h1_toPlot[i]->GetXaxis()->SetTitle("#theta (Degrees)");
      if (pi0 != 1){h1_toPlot[i]->GetYaxis()->SetTitle("Events per 2^{#circ}");}
      else {h1_toPlot[i]->GetYaxis()->SetTitle("Events per 4^{#circ}");}
      //h1_toPlot[i]->GetXaxis()->SetRangeUser( 0, 180);
      //h1_toPlot[i]->GetYaxis()->SetRangeUser( 0, 10);
      h1_toPlot[i]->GetXaxis()->CenterTitle( true);
      h1_toPlot[i]->GetYaxis()->CenterTitle( true);
  //        h1_el_momentum[i]->GetXaxis()->SetTitleSize( 0.06);
  //        h1_el_momentum[i]->GetYaxis()->SetTitleSize( 0.06);


        //h1_toPlot[i]->SetTitle(TString(energ)+" GeV " + TString(nucleus) + " (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));
        h1_toPlot[i]->SetTitle("Electron #theta Distribution of (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));

      // take the maximum of the four interactions to use later in plotting
      if (h1_toPlot[i]->GetMaximum() > hist_max) { hist_max = h1_toPlot[i]->GetMaximum(); }
  }

  // draw the histograms
  for( int i = 0; i < 4; i++){
      h1_toPlot[i]->GetYaxis()->SetRangeUser( 0, hist_max*1.1); // smithja: hist_max mult. by 1.1 to give the graph some head room
      h1_toPlot[i]->GetXaxis()->SetRangeUser( 0, 180);
      h1_toPlot[i]->Draw( draw_options[i]);
      h1_toPlot[i]->Draw( draw_options2[i]);
  }

  // crop the margins of the canvas
  c->SetLeftMargin( 0.14);
  c->SetBottomMargin( 0.14);
  c->Update();


  // smithja: this section of code integrates the distributions and prints
  //          it to the screen
  int num_DIS_events = h1_toPlot[0]->Integral();
  int num_1232_events = h1_toPlot[1]->Integral();
  int num_Other_events = h1_toPlot[2]->Integral();
  int num_Total_events = h1_toPlot[3]->Integral();
  //int num_events = num_QE_events + num_MEC_events + num_RES_events + num_DIS_events;

  // create and draw a legend for our histograms
  TLegend *legend = new TLegend( 0.62, 0.75, 1.07, 0.9);
  legend->AddEntry( h1_toPlot[0],Form("DIS (%d)",num_DIS_events),"l");
  legend->AddEntry( h1_toPlot[1],Form("#Delta(1232) (%d)",num_1232_events),"l");
  legend->AddEntry( h1_toPlot[2],Form("Other (%d)",num_Other_events),"l");
  legend->AddEntry( h1_toPlot[3],Form("Total (%d)",num_Total_events),"l");
  legend->SetBorderSize( 0);
  legend->SetFillStyle( 0);
  legend->Draw();

  // save this histogram as a PDF file
  std::string save_file_path ("invalid");
  save_file_path = "../output/h1_pion_electrontheta/"+TString(energ)+TString(nucleus)+TString(cuts)+TString(particles)+TString(numProt)+"P"+TString(a)+".pdf";

  c->SaveAs( TString::Format("%s", save_file_path.c_str()));
}
