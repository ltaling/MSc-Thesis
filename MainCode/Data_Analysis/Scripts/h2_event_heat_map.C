// This macro creates a heat map of the number of events
// The heat map I made looks terrible so I simply printed out the values to the screen and made the plots in Python
// You can either do the same, or you can improve this code so they look good in ROOT

#include <TFile.h>
#include <TH2I.h>
#include <TF2.h>
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

void h2_event_heat_map (std::string file_name, int type){
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

  TH2I* h2_heatmap;
  TH2I* h2_heatmap_temp[6];
  std::string title = "";

  h2_heatmap_temp[0] = (TH2I*)input1->Get(TString::Format("h2_threshold_passing_1")); //QES
  h2_heatmap_temp[1] = (TH2I*)input1->Get(TString::Format("h2_threshold_passing_2")); //MEC
  h2_heatmap_temp[2] = (TH2I*)input1->Get(TString::Format("h2_threshold_passing_3")); //OtherRES
  h2_heatmap_temp[3] = (TH2I*)input1->Get(TString::Format("h2_threshold_passing_4")); //DIS
  h2_heatmap_temp[4] = (TH2I*)input1->Get(TString::Format("h2_threshold_passing_5")); //1232RES
  h2_heatmap_temp[5] = (TH2I*)input1->Get(TString::Format("h2_threshold_passing_0")); //Other

  if (type == 1){
    h2_heatmap = h2_heatmap_temp[3]; //DIS
    title = title + "DIS";
  }

  else if (type == 2){
    h2_heatmap = h2_heatmap_temp[4]; //1232RES
    title = title + "#Delta(1232) Resonance";
  }

  else if (type == 3){
    h2_heatmap = (TH2I*)h2_heatmap_temp[0]->Clone();//All other ones
    h2_heatmap->Add(h2_heatmap_temp[1]);
    h2_heatmap->Add(h2_heatmap_temp[2]);
    h2_heatmap->Add(h2_heatmap_temp[5]);
    title = title + "Other";
  }

  else if (type == 4){
    h2_heatmap = (TH2I*)h2_heatmap_temp[0]->Clone(); // Total of everything
    h2_heatmap->Add(h2_heatmap_temp[1]);
    h2_heatmap->Add(h2_heatmap_temp[2]);
    h2_heatmap->Add(h2_heatmap_temp[3]);
    h2_heatmap->Add(h2_heatmap_temp[4]);
    h2_heatmap->Add(h2_heatmap_temp[5]);
    title = title + "Total";
  }

  h2_heatmap->Print("all"); // Print the all values to the screen
  // Note: Underflow bins are also printed, so any bin coordinate with a 0 in it (e.g. [0,1] or [2,0] etc.) will have a 0 entry and should be ignored

  TCanvas* c;
  c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

  h2_heatmap->SetStats(0); // get rid of the stats box that usually appears at the top right of plots
  h2_heatmap->GetXaxis()->SetTitle("Number of Pions Above Threshold");
  h2_heatmap->GetYaxis()->SetTitle("Number of Pions Below Threshold");
  //h2_heatmap->GetXaxis()->SetRangeUser( 0, 1);
  //h2_heatmap->GetYaxis()->SetRangeUser( 0, 10);
  h2_heatmap->GetXaxis()->CenterTitle( true);
  h2_heatmap->GetYaxis()->CenterTitle( true);
  TAxis* ax = h2_heatmap->GetXaxis();
  TAxis* ay = h2_heatmap->GetYaxis();
  ax->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  ax->ChangeLabel(2,-1,-1,-1,-1,-1,"1");
  ax->ChangeLabel(3,-1,-1,-1,-1,-1,"2");
  ax->ChangeLabel(4,-1,-1,-1,-1,-1,"3+");
  ay->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  ay->ChangeLabel(2,-1,-1,-1,-1,-1,"1");
  ay->ChangeLabel(3,-1,-1,-1,-1,-1,"2");
  ay->ChangeLabel(4,-1,-1,-1,-1,-1,"3+");
  h2_heatmap->LabelsDeflate("X");
  h2_heatmap->LabelsDeflate("Y");

//        h1_el_momentum[i]->GetXaxis()->SetTitleSize( 0.06);
//        h1_el_momentum[i]->GetYaxis()->SetTitleSize( 0.06);


  //h2_heatmap->SetTitle(TString(energ)+" Gev " + TString(nucleus) + " (1e"+TString(numProt)+"p1pi) Events "+TString(accept) + " for "+TString(types));
  h2_heatmap->SetTitle(TString(title));
  gStyle->SetPalette(kRust);
  h2_heatmap->Draw("COLZ");

  // save this histogram as a PDF file
  std::string save_file_path ("invalid");
  save_file_path = "../output/h2_pion_eventcount/Heatmap_"+TString(title)+".pdf";

  c->SaveAs( TString::Format("%s", save_file_path.c_str()));

}
