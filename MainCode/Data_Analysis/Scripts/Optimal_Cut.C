// Finds optimal high cut on W, Q^2 and x to purify the 1232 RES signal
// This code is an amalgamation of multiple programs and is thus very unreadable
// If you wish to modify this code, I recommend you make a new program and view the original files, which are much easier to read
// The plotting comes from h1_pion_mom.C (or any other h1_pion plotting macro)
// The optimal cut algorithm comes from Nathaniel Boyle's code, with the same name as this file, on his github (phronmist)
// Then it has been copied 3 times to accomodate for the 3 different cuts to be made

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

void Optimal_Cut(std::string file_name, bool pimin, bool pi0, bool piplus)

{
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
  //W
  TH1F* h1_w_stacked[6]; // Bins of all pion types added together
  TH1F* h1_pi_w_0[6]; //pi-
  TH1F* h1_pi_w_1[6]; //pi0
  TH1F* h1_pi_w_2[6]; //pi+
  TH1F* h1_toPlot_w[4]; //things to actually plot
  //Q^2
  TH1F* h1_q2_stacked[6]; // Bins of all pion types added together
  TH1F* h1_pi_q2_0[6]; //pi-
  TH1F* h1_pi_q2_1[6]; //pi0
  TH1F* h1_pi_q2_2[6]; //pi+
  TH1F* h1_toPlot_q2[4]; //things to actually plot
  //Bjorken x
  TH1F* h1_x_stacked[6]; // Bins of all pion types added together
  TH1F* h1_pi_x_0[6]; //pi-
  TH1F* h1_pi_x_1[6]; //pi0
  TH1F* h1_pi_x_2[6]; //pi+
  TH1F* h1_toPlot_x[4]; //things to actually plot


  // Drawing histograms
  if (pimin){
    h1_pi_w_0[0] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_1_type_0")); //QES
    h1_pi_w_0[1] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_2_type_0")); //MEC
    h1_pi_w_0[2] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_3_type_0")); //OtherRES
    h1_pi_w_0[3] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_4_type_0")); //DIS
    h1_pi_w_0[4] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_5_type_0")); //1232RES
    h1_pi_w_0[5] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_0_type_0")); //Other
    h1_pi_q2_0[0] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_1_type_0")); //QES
    h1_pi_q2_0[1] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_2_type_0")); //MEC
    h1_pi_q2_0[2] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_3_type_0")); //OtherRES
    h1_pi_q2_0[3] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_4_type_0")); //DIS
    h1_pi_q2_0[4] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_5_type_0")); //1232RES
    h1_pi_q2_0[5] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_0_type_0")); //Other
    h1_pi_x_0[0] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_1_type_0")); //QES
    h1_pi_x_0[1] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_2_type_0")); //MEC
    h1_pi_x_0[2] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_3_type_0")); //OtherRES
    h1_pi_x_0[3] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_4_type_0")); //DIS
    h1_pi_x_0[4] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_5_type_0")); //1232RES
    h1_pi_x_0[5] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_0_type_0")); //Other
    if (h1_w_stacked[2] == nullptr){
      types = types + "#pi^{-}";
      h1_w_stacked[0] = h1_pi_w_0[0];
      h1_w_stacked[1] = h1_pi_w_0[1];
      h1_w_stacked[2] = h1_pi_w_0[2];
      h1_w_stacked[3] = h1_pi_w_0[3];
      h1_w_stacked[4] = h1_pi_w_0[4];
      h1_w_stacked[5] = h1_pi_w_0[5];
      h1_q2_stacked[0] = h1_pi_q2_0[0];
      h1_q2_stacked[1] = h1_pi_q2_0[1];
      h1_q2_stacked[2] = h1_pi_q2_0[2];
      h1_q2_stacked[3] = h1_pi_q2_0[3];
      h1_q2_stacked[4] = h1_pi_q2_0[4];
      h1_q2_stacked[5] = h1_pi_q2_0[5];
      h1_x_stacked[0] = h1_pi_x_0[0];
      h1_x_stacked[1] = h1_pi_x_0[1];
      h1_x_stacked[2] = h1_pi_x_0[2];
      h1_x_stacked[3] = h1_pi_x_0[3];
      h1_x_stacked[4] = h1_pi_x_0[4];
      h1_x_stacked[5] = h1_pi_x_0[5];
      }
    else{
      types = types + " & #pi^{-}";
      h1_w_stacked[0]->Add(h1_pi_w_0[0]);
      h1_w_stacked[1]->Add(h1_pi_w_0[1]);
      h1_w_stacked[2]->Add(h1_pi_w_0[2]);
      h1_w_stacked[3]->Add(h1_pi_w_0[3]);
      h1_w_stacked[4]->Add(h1_pi_w_0[4]);
      h1_w_stacked[5]->Add(h1_pi_w_0[5]);
      h1_q2_stacked[0]->Add(h1_pi_q2_0[0]);
      h1_q2_stacked[1]->Add(h1_pi_q2_0[1]);
      h1_q2_stacked[2]->Add(h1_pi_q2_0[2]);
      h1_q2_stacked[3]->Add(h1_pi_q2_0[3]);
      h1_q2_stacked[4]->Add(h1_pi_q2_0[4]);
      h1_q2_stacked[5]->Add(h1_pi_q2_0[5]);
      h1_x_stacked[0]->Add(h1_pi_x_0[0]);
      h1_x_stacked[1]->Add(h1_pi_x_0[1]);
      h1_x_stacked[2]->Add(h1_pi_x_0[2]);
      h1_x_stacked[3]->Add(h1_pi_x_0[3]);
      h1_x_stacked[4]->Add(h1_pi_x_0[4]);
      h1_x_stacked[5]->Add(h1_pi_x_0[5]);
    }
  }
  if (pi0){
    h1_pi_w_1[0] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_1_type_1")); //QES
    h1_pi_w_1[1] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_2_type_1")); //MEC
    h1_pi_w_1[2] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_3_type_1")); //OtherRES
    h1_pi_w_1[3] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_4_type_1")); //DIS
    h1_pi_w_1[4] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_5_type_1")); //1232RES
    h1_pi_w_1[5] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_0_type_1")); //Other
    h1_pi_q2_1[0] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_1_type_1")); //QES
    h1_pi_q2_1[1] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_2_type_1")); //MEC
    h1_pi_q2_1[2] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_3_type_1")); //OtherRES
    h1_pi_q2_1[3] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_4_type_1")); //DIS
    h1_pi_q2_1[4] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_5_type_1")); //1232RES
    h1_pi_q2_1[5] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_0_type_1")); //Other
    h1_pi_x_1[0] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_1_type_1")); //QES
    h1_pi_x_1[1] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_2_type_1")); //MEC
    h1_pi_x_1[2] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_3_type_1")); //OtherRES
    h1_pi_x_1[3] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_4_type_1")); //DIS
    h1_pi_x_1[4] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_5_type_1")); //1232RES
    h1_pi_x_1[5] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_0_type_1")); //Other
    if (h1_w_stacked[2] == nullptr){
      types = types + "#pi^{0}";
      h1_w_stacked[0] = h1_pi_w_1[0];
      h1_w_stacked[1] = h1_pi_w_1[1];
      h1_w_stacked[2] = h1_pi_w_1[2];
      h1_w_stacked[3] = h1_pi_w_1[3];
      h1_w_stacked[4] = h1_pi_w_1[4];
      h1_w_stacked[5] = h1_pi_w_1[5];
      h1_q2_stacked[0] = h1_pi_q2_1[0];
      h1_q2_stacked[1] = h1_pi_q2_1[1];
      h1_q2_stacked[2] = h1_pi_q2_1[2];
      h1_q2_stacked[3] = h1_pi_q2_1[3];
      h1_q2_stacked[4] = h1_pi_q2_1[4];
      h1_q2_stacked[5] = h1_pi_q2_1[5];
      h1_x_stacked[0] = h1_pi_x_1[0];
      h1_x_stacked[1] = h1_pi_x_1[1];
      h1_x_stacked[2] = h1_pi_x_1[2];
      h1_x_stacked[3] = h1_pi_x_1[3];
      h1_x_stacked[4] = h1_pi_x_1[4];
      h1_x_stacked[5] = h1_pi_x_1[5];
      }
    else{
      types = types + " & #pi^{0}";
      h1_w_stacked[0]->Add(h1_pi_w_1[0]);
      h1_w_stacked[1]->Add(h1_pi_w_1[1]);
      h1_w_stacked[2]->Add(h1_pi_w_1[2]);
      h1_w_stacked[3]->Add(h1_pi_w_1[3]);
      h1_w_stacked[4]->Add(h1_pi_w_1[4]);
      h1_w_stacked[5]->Add(h1_pi_w_1[5]);
      h1_q2_stacked[0]->Add(h1_pi_q2_1[0]);
      h1_q2_stacked[1]->Add(h1_pi_q2_1[1]);
      h1_q2_stacked[2]->Add(h1_pi_q2_1[2]);
      h1_q2_stacked[3]->Add(h1_pi_q2_1[3]);
      h1_q2_stacked[4]->Add(h1_pi_q2_1[4]);
      h1_q2_stacked[5]->Add(h1_pi_q2_1[5]);
      h1_x_stacked[0]->Add(h1_pi_x_1[0]);
      h1_x_stacked[1]->Add(h1_pi_x_1[1]);
      h1_x_stacked[2]->Add(h1_pi_x_1[2]);
      h1_x_stacked[3]->Add(h1_pi_x_1[3]);
      h1_x_stacked[4]->Add(h1_pi_x_1[4]);
      h1_x_stacked[5]->Add(h1_pi_x_1[5]);
    }
  }
  if (piplus){
    h1_pi_w_2[0] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_1_type_2")); //QES
    h1_pi_w_2[1] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_2_type_2")); //MEC
    h1_pi_w_2[2] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_3_type_2")); //OtherRES
    h1_pi_w_2[3] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_4_type_2")); //DIS
    h1_pi_w_2[4] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_5_type_2")); //1232RES
    h1_pi_w_2[5] = (TH1F*)input1->Get(TString::Format("h1_pion_w_interaction_0_type_2")); //Other
    h1_pi_q2_2[0] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_1_type_2")); //QES
    h1_pi_q2_2[1] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_2_type_2")); //MEC
    h1_pi_q2_2[2] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_3_type_2")); //OtherRES
    h1_pi_q2_2[3] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_4_type_2")); //DIS
    h1_pi_q2_2[4] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_5_type_2")); //1232RES
    h1_pi_q2_2[5] = (TH1F*)input1->Get(TString::Format("h1_pion_q2_interaction_0_type_2")); //Other
    h1_pi_x_2[0] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_1_type_2")); //QES
    h1_pi_x_2[1] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_2_type_2")); //MEC
    h1_pi_x_2[2] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_3_type_2")); //OtherRES
    h1_pi_x_2[3] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_4_type_2")); //DIS
    h1_pi_x_2[4] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_5_type_2")); //1232RES
    h1_pi_x_2[5] = (TH1F*)input1->Get(TString::Format("h1_pion_x_interaction_0_type_2")); //Other
    if (h1_w_stacked[2] == nullptr){
      types = types + "#pi^{+}";
      h1_w_stacked[0] = h1_pi_w_2[0];
      h1_w_stacked[1] = h1_pi_w_2[1];
      h1_w_stacked[2] = h1_pi_w_2[2];
      h1_w_stacked[3] = h1_pi_w_2[3];
      h1_w_stacked[4] = h1_pi_w_2[4];
      h1_w_stacked[5] = h1_pi_w_2[5];
      h1_q2_stacked[0] = h1_pi_q2_2[0];
      h1_q2_stacked[1] = h1_pi_q2_2[1];
      h1_q2_stacked[2] = h1_pi_q2_2[2];
      h1_q2_stacked[3] = h1_pi_q2_2[3];
      h1_q2_stacked[4] = h1_pi_q2_2[4];
      h1_q2_stacked[5] = h1_pi_q2_2[5];
      h1_x_stacked[0] = h1_pi_x_2[0];
      h1_x_stacked[1] = h1_pi_x_2[1];
      h1_x_stacked[2] = h1_pi_x_2[2];
      h1_x_stacked[3] = h1_pi_x_2[3];
      h1_x_stacked[4] = h1_pi_x_2[4];
      h1_x_stacked[5] = h1_pi_x_2[5];
      }
    else{
      types = types + " & #pi^{+}";
      h1_w_stacked[0]->Add(h1_pi_w_2[0]);
      h1_w_stacked[1]->Add(h1_pi_w_2[1]);
      h1_w_stacked[2]->Add(h1_pi_w_2[2]);
      h1_w_stacked[3]->Add(h1_pi_w_2[3]);
      h1_w_stacked[4]->Add(h1_pi_w_2[4]);
      h1_w_stacked[5]->Add(h1_pi_w_2[5]);
      h1_q2_stacked[0]->Add(h1_pi_q2_2[0]);
      h1_q2_stacked[1]->Add(h1_pi_q2_2[1]);
      h1_q2_stacked[2]->Add(h1_pi_q2_2[2]);
      h1_q2_stacked[3]->Add(h1_pi_q2_2[3]);
      h1_q2_stacked[4]->Add(h1_pi_q2_2[4]);
      h1_q2_stacked[5]->Add(h1_pi_q2_2[5]);
      h1_x_stacked[0]->Add(h1_pi_x_2[0]);
      h1_x_stacked[1]->Add(h1_pi_x_2[1]);
      h1_x_stacked[2]->Add(h1_pi_x_2[2]);
      h1_x_stacked[3]->Add(h1_pi_x_2[3]);
      h1_x_stacked[4]->Add(h1_pi_x_2[4]);
      h1_x_stacked[5]->Add(h1_pi_x_2[5]);
    }
  }

    //Create the things to plot
    h1_toPlot_w[0] = h1_w_stacked[3]; //DIS
    h1_toPlot_w[1] = h1_w_stacked[4]; //1232RES
    h1_toPlot_w[2] = (TH1F*)h1_w_stacked[0]->Clone();//All other ones
    h1_toPlot_w[2]->Add(h1_w_stacked[1]);
    h1_toPlot_w[2]->Add(h1_w_stacked[2]);
    h1_toPlot_w[2]->Add(h1_w_stacked[5]);
    h1_toPlot_w[3] = (TH1F*)h1_w_stacked[0]->Clone(); // Total of everything
    h1_toPlot_w[3]->Add(h1_w_stacked[1]);
    h1_toPlot_w[3]->Add(h1_w_stacked[2]);
    h1_toPlot_w[3]->Add(h1_w_stacked[3]);
    h1_toPlot_w[3]->Add(h1_w_stacked[4]);
    h1_toPlot_w[3]->Add(h1_w_stacked[5]);

    h1_toPlot_q2[0] = h1_q2_stacked[3]; //DIS
    h1_toPlot_q2[1] = h1_q2_stacked[4]; //1232RES
    h1_toPlot_q2[2] = (TH1F*)h1_q2_stacked[0]->Clone();//All other ones
    h1_toPlot_q2[2]->Add(h1_q2_stacked[1]);
    h1_toPlot_q2[2]->Add(h1_q2_stacked[2]);
    h1_toPlot_q2[2]->Add(h1_q2_stacked[5]);
    h1_toPlot_q2[3] = (TH1F*)h1_q2_stacked[0]->Clone(); // Total of everything
    h1_toPlot_q2[3]->Add(h1_q2_stacked[1]);
    h1_toPlot_q2[3]->Add(h1_q2_stacked[2]);
    h1_toPlot_q2[3]->Add(h1_q2_stacked[3]);
    h1_toPlot_q2[3]->Add(h1_q2_stacked[4]);
    h1_toPlot_q2[3]->Add(h1_q2_stacked[5]);

    h1_toPlot_x[0] = h1_x_stacked[3]; //DIS
    h1_toPlot_x[1] = h1_x_stacked[4]; //1232RES
    h1_toPlot_x[2] = (TH1F*)h1_x_stacked[0]->Clone();//All other ones
    h1_toPlot_x[2]->Add(h1_x_stacked[1]);
    h1_toPlot_x[2]->Add(h1_x_stacked[2]);
    h1_toPlot_x[2]->Add(h1_x_stacked[5]);
    h1_toPlot_x[3] = (TH1F*)h1_x_stacked[0]->Clone(); // Total of everything
    h1_toPlot_x[3]->Add(h1_x_stacked[1]);
    h1_toPlot_x[3]->Add(h1_x_stacked[2]);
    h1_toPlot_x[3]->Add(h1_x_stacked[3]);
    h1_toPlot_x[3]->Add(h1_x_stacked[4]);
    h1_toPlot_x[3]->Add(h1_x_stacked[5]);

    //Creating the names of each thing
    std::string particles = "";
    if (pimin){particles = particles + "Pimin";}
    if (pi0){particles = particles + "Pizero";}
    if (piplus){particles = particles + "Piplus";}


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

    int steps_w = h1_toPlot_w[1]->GetNbinsX();
    int steps_q2 = h1_toPlot_q2[1]->GetNbinsX();
    int steps_x = h1_toPlot_x[1]->GetNbinsX();

    //Initialize arrays for plotting
    float pur_x_eff_w[steps_w]; //Purity times Efficiency  (y value of plot)
    float Eff_new_w[steps_w];   //Efficiency (y value of plot)
    float Pur_new_w[steps_w];   //Purity (y value of plot)
    float pur_x_eff_q2[steps_q2]; //Purity times Efficiency  (y value of plot)
    float Eff_new_q2[steps_q2];   //Efficiency (y value of plot)
    float Pur_new_q2[steps_q2];   //Purity (y value of plot)
    float pur_x_eff_x[steps_x]; //Purity times Efficiency  (y value of plot)
    float Eff_new_x[steps_x];   //Efficiency (y value of plot)
    float Pur_new_x[steps_x];   //Purity (y value of plot)


    float resid_total_w = h1_toPlot_w[1]->Integral(); // # total weighted Delta 1232 events from histogram
    float resid_total_q2 = h1_toPlot_q2[1]->Integral();
    float resid_total_x = h1_toPlot_x[1]->Integral();

    float w_center[steps_w];
    float q2_center[steps_q2];
    float x_center[steps_x];

    // ------------------------------------------------- Start of finding purity and efficiency -------------------------------------------------

    //Calculate purity and efficiency for W
    for (int i = 0; i < steps_w; ++i)
    {
        float resid_cut_new_w = h1_toPlot_w[1]->Integral(0,i); // calculate Delta 1232 events within new cut
        Eff_new_w[i] = resid_cut_new_w / resid_total_w; // calculate Efficiency
        float n_total_cut_new_w = h1_toPlot_w[3]->Integral(0,i);    //Save new signal + background as a float
        w_center[i] = h1_toPlot_w[1]->GetBinCenter(i);

        if (n_total_cut_new_w == 0) //avoid dividing by zero
            Pur_new_w[i] = 0;
        else
            Pur_new_w[i] = resid_cut_new_w / n_total_cut_new_w; // calculate purity
        pur_x_eff_w[i] = Eff_new_w[i] * Pur_new_w[i]; // calculate new element
    }

    //Calculate purity and efficiency for Q^2 -- placing a high cut
    for (int i = 0; i < steps_q2; ++i)
    {
        float resid_cut_new_q2 = h1_toPlot_q2[1]->Integral(0,i); // calculate Delta 1232 events within new cut
        Eff_new_q2[i] = resid_cut_new_q2 / resid_total_q2; // calculate Efficiency
        float n_total_cut_new_q2 = h1_toPlot_q2[3]->Integral(0,i);    //Save new signal + background as a float
        q2_center[i] = h1_toPlot_q2[1]->GetBinCenter(i);

        if (n_total_cut_new_q2 == 0) //avoid dividing by zero
            Pur_new_q2[i] = 0;
        else
            Pur_new_q2[i] = resid_cut_new_q2 / n_total_cut_new_q2; // calculate purity
        pur_x_eff_q2[i] = Eff_new_q2[i] * Pur_new_q2[i]; // calculate new element
    }

    //Calculate purity and efficiency for Bjorken x -- placing a low cut
    for (int i = steps_x; i-- > 0; ) //Counting down in the for loop
    {
        float resid_cut_new_x = h1_toPlot_x[1]->Integral(i,steps_x); // calculate Delta 1232 events within new cut
        Eff_new_x[i] = resid_cut_new_x / resid_total_x; // calculate Efficiency
        float n_total_cut_new_x = h1_toPlot_x[3]->Integral(i,steps_x);    //Save new signal + background as a float
        x_center[i] = h1_toPlot_x[1]->GetBinCenter(i);

        if (n_total_cut_new_x == 0) //avoid dividing by zero
            Pur_new_x[i] = 0;
        else
            Pur_new_x[i] = resid_cut_new_x / n_total_cut_new_x; // calculate purity
        pur_x_eff_x[i] = Eff_new_x[i] * Pur_new_x[i]; // calculate new element
    }

    // --------------------------------------------------- End of finding purity and efficiency -------------------------------

    // --------------------------------------------------- Start of plotting --------------------------------------------------

    // create graph of purity and efficiency as a function of W,q2 and x

    TGraph* gr_w = new TGraph(steps_w, w_center, pur_x_eff_w);// Purity x Efficiency graph
    TGraph* gr2_w = new TGraph(steps_w, w_center, Eff_new_w); //Efficiency Graph
    TGraph* gr3_w = new TGraph(steps_w, w_center, Pur_new_w); //Purity Graph
    gr_w->GetXaxis()->SetTitle("High Cut on W"); //Label x Axis
    gr_w->GetXaxis()->CenterTitle(true); //Center x Axis
    gr_w->GetYaxis()->SetTitle("Value"); //Label y Axis
    gr_w->GetYaxis()->CenterTitle(true); //Center y Axis

    TGraph* gr_q2 = new TGraph(steps_q2, q2_center, pur_x_eff_q2);// Purity x Efficiency graph
    TGraph* gr2_q2 = new TGraph(steps_q2, q2_center, Eff_new_q2); //Efficiency Graph
    TGraph* gr3_q2 = new TGraph(steps_q2, q2_center, Pur_new_q2); //Purity Graph
    gr_q2->GetXaxis()->SetTitle("High Cut on Q^{2}"); //Label x Axis
    gr_q2->GetXaxis()->CenterTitle(true); //Center x Axis
    gr_q2->GetYaxis()->SetTitle("Value"); //Label y Axis
    gr_q2->GetYaxis()->CenterTitle(true); //Center y Axis

    TGraph* gr_x = new TGraph(steps_x, x_center, pur_x_eff_x);// Purity x Efficiency graph
    TGraph* gr2_x = new TGraph(steps_x, x_center, Eff_new_x); //Efficiency Graph
    TGraph* gr3_x = new TGraph(steps_x, x_center, Pur_new_x); //Purity Graph
    gr_x->GetXaxis()->SetTitle("Low Cut on Bjorken x"); //Label x Axis
    gr_x->GetXaxis()->CenterTitle(true); //Center x Axis
    gr_x->GetYaxis()->SetTitle("Value"); //Label y Axis
    gr_x->GetYaxis()->CenterTitle(true); //Center y Axis

    //set color of graphs

    gr_w->SetLineColor(kRed);
    //gr_w->SetMarkerSize(0);
    gr2_w->SetLineColor(kBlue);
    //gr2_w->SetMarkerSize(0);
    gr3_w->SetLineColor(kGreen);
    //gr3_w->SetMarkerSize(0);

    gr_q2->SetLineColor(kRed);
    //gr_q2->SetMarkerSize(0);
    gr2_q2->SetLineColor(kBlue);
    //gr2_q2->SetMarkerSize(0);
    gr3_q2->SetLineColor(kGreen);
    //gr3_q2->SetMarkerSize(0);

    gr_x->SetLineColor(kRed);
    //gr_x->SetMarkerSize(0);
    gr2_x->SetLineColor(kBlue);
    //gr2_x->SetMarkerSize(0);
    gr3_x->SetLineColor(kGreen);
    //gr3_x->SetMarkerSize(0);

    //Create canvas

    TCanvas* c1;
    c1 = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    TMultiGraph* mg_w = new TMultiGraph();
    mg_w->Add(gr_w);
    mg_w->Add(gr2_w);
    mg_w->Add(gr3_w);
    mg_w->SetTitle("Optimal W cut of (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));
    //mg_w->SetTitle("Optimal Cut on W for 4.461 GeV Electron Beam on {}^{12}C; High Cut on W; Value");
    //mg->GetHistogram()->GetXaxis()->CenterTitle( true);
    //mg->GetHistogram()->GetYaxis()->CenterTitle( true);


    //find where purity x efficiency is max
    const int N_w = sizeof(pur_x_eff_w) / sizeof(int);
    int max_index_w = distance(pur_x_eff_w, max_element(pur_x_eff_w, pur_x_eff_w + N_w));

    // Create solid line at max value of pur x eff
    TLine* l_w = new TLine(h1_toPlot_w[1]->GetBinCenter(max_index_w), 0, h1_toPlot_w[1]->GetBinCenter(max_index_w), 1);
    l_w->SetLineColor(kBlack);
    //l->SetLineStyle(1);
    l_w->SetLineWidth(2);
    mg_w->Draw("AL");
    l_w->Draw("SAME");

    // crop the margins of the canvas
    c1->SetLeftMargin( 0.14);
    c1->SetBottomMargin( 0.14);
    c1->Update();

    //Print optimal values in legend as strings

    std::string textW = "Optimal Cut: W =  ";//print value of optimal W
    textW += to_string(h1_toPlot_w[1]->GetBinCenter(max_index_w));

    std::string textPur_w = "Purity (Optimal Purity: "; //print value of optimal purity
        textPur_w += to_string(Pur_new_w[max_index_w]);
        textPur_w += ")";

        std::string textPxE_w = "Purity x Efficiency (Optimal Value: ";//print value of optimal purity times efficiency
        textPxE_w += to_string(pur_x_eff_w[max_index_w]);
        textPxE_w += ")";

        std::string textEff_w = "Efficiency (Optimal Efficiency: ";//print value of optimal efficiency
        textEff_w += to_string(Eff_new_w[max_index_w]);
        textEff_w += ")";

    // Create legend

    auto legend_w = new TLegend(0.53, 0.35, 0.9, 0.6);
    //legend_w->SetHeader("Legend", "C");
    legend_w->AddEntry(gr_w, Form("%s",textPxE_w.c_str()), "l");
    legend_w->AddEntry(gr2_w, Form("%s",textEff_w.c_str()), "l");
    legend_w->AddEntry(gr3_w, Form("%s",textPur_w.c_str()), "l");
    legend_w->AddEntry(l_w, Form("%s",textW.c_str()), "l");
    legend_w->SetTextSize(0.02);
    legend_w->Draw();

    // save this histogram as a PDF file
    std::string save_file_path_w ("invalid");
    save_file_path_w = "../output/Optimal_cuts/W_cuts_"+TString(energ)+TString(nucleus)+TString(cuts)+TString(particles)+TString(numProt)+"P"+TString(a)+".pdf";

    c1->SaveAs( TString::Format("%s", save_file_path_w.c_str()));

    //Create new canvas and graph for Q^2

    TCanvas* c2;
    c2 = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    TMultiGraph* mg_q2 = new TMultiGraph();
    mg_q2->Add(gr_q2);
    mg_q2->Add(gr2_q2);
    mg_q2->Add(gr3_q2);
    mg_q2->SetTitle("Optimal Q^{2} cut of (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));
    //mg_q2->SetTitle("Optimal Cut on W for 4.461 GeV Electron Beam on {}^{12}C; High Cut on W; Value");
    //mg->GetHistogram()->GetXaxis()->CenterTitle( true);
    //mg->GetHistogram()->GetYaxis()->CenterTitle( true);


    //find where purity x efficiency is max
    const int N_q2 = sizeof(pur_x_eff_q2) / sizeof(int);
    int max_index_q2 = distance(pur_x_eff_q2, max_element(pur_x_eff_q2, pur_x_eff_q2 + N_q2));

    // Create solid line at max value of pur x eff
    TLine* l_q2 = new TLine(h1_toPlot_q2[1]->GetBinCenter(max_index_q2), 0, h1_toPlot_q2[1]->GetBinCenter(max_index_q2), 1);
    l_q2->SetLineColor(kBlack);
    //l->SetLineStyle(1);
    l_q2->SetLineWidth(2);
    mg_q2->Draw("AL");
    l_q2->Draw("SAME");

    // crop the margins of the canvas
    c2->SetLeftMargin( 0.14);
    c2->SetBottomMargin( 0.14);
    c2->Update();

    //Print optimal values in legend as strings

    std::string textq2 = "Optimal Cut: Q^{2} =  ";//print value of optimal Q^2
    textq2 += to_string(h1_toPlot_q2[1]->GetBinCenter(max_index_q2));

    std::string textPur_q2 = "Purity (Optimal Purity: "; //print value of optimal purity
        textPur_q2 += to_string(Pur_new_q2[max_index_q2]);
        textPur_q2 += ")";

        std::string textPxE_q2 = "Purity x Efficiency (Optimal Value: ";//print value of optimal purity times efficiency
        textPxE_q2 += to_string(pur_x_eff_q2[max_index_q2]);
        textPxE_q2 += ")";

        std::string textEff_q2 = "Efficiency (Optimal Efficiency: ";//print value of optimal efficiency
        textEff_q2 += to_string(Eff_new_q2[max_index_q2]);
        textEff_q2 += ")";

    // Create legend

    auto legend_q2 = new TLegend(0.17, 0.65, 0.47, 0.9);
    legend_q2->SetHeader("Legend", "C");
    legend_q2->AddEntry(gr_q2, Form("%s",textPxE_q2.c_str()), "l");
    legend_q2->AddEntry(gr2_q2, Form("%s",textEff_q2.c_str()), "l");
    legend_q2->AddEntry(gr3_q2, Form("%s",textPur_q2.c_str()), "l");
    legend_q2->AddEntry(l_q2, Form("%s",textq2.c_str()), "l");
    legend_q2->Draw();

    // save this histogram as a PDF file
    std::string save_file_path_q2 ("invalid");
    save_file_path_q2 = "../output/Optimal_cuts/q2_cuts_"+TString(energ)+TString(nucleus)+TString(cuts)+TString(particles)+TString(numProt)+"P"+TString(a)+".pdf";

    c2->SaveAs( TString::Format("%s", save_file_path_q2.c_str()));


    //Drawing graph for Bjorken x
    TCanvas* c3;
    c3 = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    TMultiGraph* mg_x = new TMultiGraph();
    mg_x->Add(gr_x);
    mg_x->Add(gr2_x);
    mg_x->Add(gr3_x);
    mg_x->SetTitle("Optimal Bjorken x cut of (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));
    //mg_x->SetTitle("Optimal Cut on W for 4.461 GeV Electron Beam on {}^{12}C; High Cut on W; Value");
    //mg->GetHistogram()->GetXaxis()->CenterTitle( true);
    //mg->GetHistogram()->GetYaxis()->CenterTitle( true);


    //find where purity x efficiency is max
    const int N_x = sizeof(pur_x_eff_x) / sizeof(int);
    int max_index_x = distance(pur_x_eff_x, max_element(pur_x_eff_x, pur_x_eff_x + N_x));

    // Create solid line at max value of pur x eff
    TLine* l_x = new TLine(h1_toPlot_x[1]->GetBinCenter(max_index_x), 0, h1_toPlot_x[1]->GetBinCenter(max_index_x), 1);
    l_x->SetLineColor(kBlack);
    //l->SetLineStyle(1);
    l_x->SetLineWidth(2);
    mg_x->Draw("AL");
    l_x->Draw("SAME");

    // crop the margins of the canvas
    c3->SetLeftMargin( 0.14);
    c3->SetBottomMargin( 0.14);
    c3->Update();

    //Print optimal values in legend as strings

    std::string textx = "Optimal Cut: x =  ";//print value of optimal W
    textx += to_string(h1_toPlot_x[1]->GetBinCenter(max_index_x));

    std::string textPur_x = "Purity (Optimal Purity: "; //print value of optimal purity
        textPur_x += to_string(Pur_new_x[max_index_x]);
        textPur_x += ")";

        std::string textPxE_x = "Purity x Efficiency (Optimal Value: ";//print value of optimal purity times efficiency
        textPxE_x += to_string(pur_x_eff_x[max_index_x]);
        textPxE_x += ")";

        std::string textEff_x = "Efficiency (Optimal Efficiency: ";//print value of optimal efficiency
        textEff_x += to_string(Eff_new_x[max_index_x]);
        textEff_x += ")";

    // Create legend

    auto legend_x = new TLegend(0.53, 0.65, 0.9, 0.9);
    //legend_x->SetHeader("Legend", "C");
    legend_x->AddEntry(gr_x, Form("%s",textPxE_x.c_str()), "l");
    legend_x->AddEntry(gr2_x, Form("%s",textEff_x.c_str()), "l");
    legend_x->AddEntry(gr3_x, Form("%s",textPur_x.c_str()), "l");
    legend_x->AddEntry(l_x, Form("%s",textx.c_str()), "l");
    legend_x->SetTextSize(0.02);
    legend_x->Draw();

    // save this histogram as a PDF file
    std::string save_file_path_x ("invalid");
    save_file_path_x = "../output/Optimal_cuts/x_cuts_"+TString(energ)+TString(nucleus)+TString(cuts)+TString(particles)+TString(numProt)+"P"+TString(a)+".pdf";

    c3->SaveAs( TString::Format("%s", save_file_path_x.c_str()));
    //-------------------------------------------------- End of plotting --------------------------------------------------------------

    // ------------------------------------------------- Start of finding range and printing ------------------------------------------

    //Calculate the values of w,Q^2 and x between which the purity x efficiency is 75% its maximum value
    int w_range[2] = {0,0}; //lower[0] and upper[1] index of W
    int q2_range[2] = {0,0};
    int x_range[2] = {0,0};

    //Looping over to find the ranges
    for (int i = 0; i < steps_w; i++)
    {
      if (pur_x_eff_w[i] > (0.25*pur_x_eff_w[max_index_w]) && w_range[0] == 0) {w_range[0] = i;}
      else if (pur_x_eff_w[i] > (0.25*pur_x_eff_w[max_index_w]) && w_range[0] != 0) {w_range[1] = i;}
    }

    for (int i = 0; i < steps_q2; i++)
    {
      if (pur_x_eff_q2[i] > (0.25*pur_x_eff_q2[max_index_q2]) && q2_range[0] == 0) {q2_range[0] = i;}
      else if (pur_x_eff_q2[i] > (0.25*pur_x_eff_q2[max_index_q2]) && q2_range[0] != 0) {q2_range[1] = i;}
    }

    for (int i = 0; i < steps_x; i++)
    {
      if (pur_x_eff_x[i] > (0.25*pur_x_eff_x[max_index_x]) && x_range[0] == 0) {x_range[0] = i;}
      else if (pur_x_eff_x[i] > (0.25*pur_x_eff_x[max_index_x]) && x_range[0] != 0) {x_range[1] = i;}
    }

    //Print out the values
    cout << "Minimum W = " << w_center[w_range[0]] << endl;
    cout << "Maximum W = " << w_center[w_range[1]] << endl;
    cout << "-----------------" << endl;
    cout << "Minimum Q^2 = " << q2_center[q2_range[0]] << endl;
    cout << "Maximum Q^2 = " << q2_center[q2_range[1]] << endl;
    cout << "-----------------" << endl;
    cout << "Minimum x = " << x_center[x_range[0]] << endl;
    cout << "Maximum x = " << x_center[x_range[1]] << endl;





}
