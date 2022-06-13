#include "TCanvas.h"
#include "TH1.h" 
#include "TLegend.h" 
#include <iostream>
#include <stdio.h>
#include "Constants.h"

using namespace std;

TTree* tree(){
  // This function opens the file
  gStyle->SetOptStat(false);  // Disable status box

  TFile* inputFile = new TFile("e4nu_test.root");  // Open ROOT file as input
  inputFile->ls();

  TTree* tree = (TTree*)inputFile->Get("gst");    // Load the ntuple from the file.
    
  return tree;
}

float PlotOptimalCut(TTree* tree){
  // Create cuts (default is 1 pion (plus or minus)) NOTE: we are not including pion 0 (nfpi0) in this case unless you want to
  string qel_cut = "qel  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";    // QE cuts
  string dis_cut = "dis  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";    // DIS cuts
  string mec_cut = "mec  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";    // MEC cuts
  string res_cut = "res  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";    // RES cuts
  string resid_0_cut = "res && resid == 0  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";  // Delta 1232 Cuts
  string resid_1_cut = "res && resid == 1  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";  // kS11_1535 RES cuts
  string resid_2_cut = "res && resid == 2  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";  // kD13_1520 RES cuts
  string resid_3_cut = "res && resid == 3  && (nfpip + nfpim  == 1) && nfpi0 == 0 ";  // kS11_1650 RES cuts

  //Initialize W cut parameters
  float Max_W_cut = 1.8;  // Value of W to stop calculating
  float Min_W_cut = 1;    // Value of W to begin calculating
  float step_size = 0.01; // Step size of W for calculation
  int steps = (Max_W_cut - Min_W_cut) * (1 / step_size) + 1;  // Number of steps between max W and min W (notice this is an integer and is subject to rounding)

  //Initialize arrays for plotting
  float pur_x_eff[steps]; //Purity times Efficiency  (y value of plot)
  float Eff_new[steps];   //Efficiency (y value of plot)
  float Pur_new[steps];   //Purity (y value of plot)
  float W[steps];         //Value of W (x value of plot)

  //Fill in W array with appropriate values (1.1, 1.2 etc.)
  for (int i = 0; i < steps; ++i)
    {
      W[i] = Min_W_cut + (i * step_size);
    }

  float resid_total = tree->GetEntries(resid_0_cut.c_str()); // # total Delta 1232 events from root file

  //Calculate purity and efficiency for loop
  for (int i = 0; i < steps; ++i)
    {
      string W_0_cut = Form("resid ==0  && (nfpip + nfpim == 1) && nfpi0 == 0 && W<=%f", W[i]); //define cut on W for every step for signal
      float resid_cut_new = tree->GetEntries(W_0_cut.c_str()); // calculate Delta 1232 events within new cut
      Eff_new[i] = resid_cut_new / resid_total; // calculate Efficiency
      string signal_background = Form("(nfpip + nfpim == 1) && nfpi0 == 0 &&  W<=%f", W[i]);// calculate signal background
      float n_total_cut_new = tree->GetEntries(signal_background.c_str());    //Save new signal background as a float

      if (n_total_cut_new == 0) //avoid dividing by zero
	Pur_new[i] = 0;
      else
	Pur_new[i] = resid_cut_new / n_total_cut_new; // calculate purity
      pur_x_eff[i] = Eff_new[i] * Pur_new[i]; // calculate new element
    }

  //find where purity x efficiency is max
  const int N = sizeof(pur_x_eff) / sizeof(int);
  int max_index = distance(pur_x_eff, max_element(pur_x_eff, pur_x_eff + N));
    
  return W[max_index]; //Save best W as float 
}

void draw(TTree* tree, float Wopt){

  //convert float to string
  string Woptimum = to_string(Wopt);
  
  // Set variables for different scattering
  Bool_t qel, mec, dis, res;
  tree->SetBranchAddress("qel",&qel);
  tree->SetBranchAddress("mec",&mec);
  tree->SetBranchAddress("dis",&dis);
  tree->SetBranchAddress("res",&res);

  //Outgoing kinematics
  Double_t el;
  tree->SetBranchAddress("El", &el);

  // Cuts for only 1 pion in the final state (for neutrino files consider adding "cc" in these cuts as well)
  string qel_cut = "qel  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1";
  string dis_cut = "dis  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1";
  string mec_cut = "mec  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1";
  string res_cut = "res  && resid != 0 && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1"; //note we are not including nfpi0
  string resd_cut = "res && resid == 0  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1";// Delta 1232 RES events
  //Now for energy cuts included
  string qel_Ecut = Form("qel  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1 && W <= (%s)",Woptimum.c_str());
  string dis_Ecut = Form("dis  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1 && W <= (%s)",Woptimum.c_str());
  string mec_Ecut = Form("mec  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1&& W <= (%s)",Woptimum.c_str());
  string res_Ecut = Form("res  && resid != 0 && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1&& W <= (%s)",Woptimum.c_str()); //note we are not including nfpi0
  string resd_Ecut = Form("res && resid == 0  && (nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1&& W <= (%s)",Woptimum.c_str());// Delta 1232 RES events

  // Create the histograms
  TH1D *hcal_qel = new TH1D("cal_qel", "Calorimetric energy reconstruction in QE events", 100, 0, 2);
  TH1D *hcal_res = new TH1D("cal_res", "Calorimetric energy reconstruction in other RES events", 100, 0, 2);
  TH1D *hcal_dres = new TH1D("cal_dres", "Calorimetric energy reconstruction in #Delta 1232 RES events", 100, 0, 2);
  TH1D *hcal_dis = new TH1D("cal_dis", "Calorimetric energy reconstruction in DIS events", 100, 0, 2);
  TH1D *hcal_mec = new TH1D("cal_mec", "Calorimetric energy reconstruction in MEC events", 100, 0, 2);
  TH1D *hcal_Eqel = new TH1D("Ecal_qel", "Calorimetric energy reconstruction in QE events with cuts", 100, 0, 2);
  TH1D *hcal_Eres = new TH1D("Ecal_res", "Calorimetric energy reconstruction in other RES events with cuts", 100, 0, 2);
  TH1D *hcal_Edres = new TH1D("Ecal_dres", "Calorimetric energy reconstruction in #Delta 1232 RES events with cuts", 100, 0, 2);
  TH1D *hcal_Edis = new TH1D("Ecal_dis", "Calorimetric energy reconstruction in DIS events with cuts", 100, 0, 2);
  TH1D *hcal_Emec = new TH1D("Ecal_mec", "Calorimetric energy reconstruction in MEC events with cuts", 100, 0, 2);


  //Cut for the events we care about (No energy cut first)
  tree->Draw(">>elist", "(nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1", "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  //With energy cut here
  tree->Draw(">>Eelist", Form("(nfpip + nfpim  == 1) && nfpi0 == 0 && nfp == 1 && W < (%s)", Woptimum.c_str()), "entrylist");
  TEntryList *Eelist = (TEntryList*)gDirectory->Get("Eelist");

  /* THIS DOESN'T WORK???
  // Find maximum number of nf
  Int_t nf;
  tree->SetBranchAddress("nf",&nf);
  int N = 0;
  for(Long64_t j = 0; j < elist->GetN(); j++){
  tree->GetEntry(j);
  if(nf > N){
  N = nf;
  }
  }
  */

  //Set up hadron labelling
  Int_t pdg[13], resid, nf;
  Double_t E[13]; //Hadron energy
  tree->SetBranchAddress("nf",&nf);
  tree->SetBranchAddress("pdgf",&pdg);
  tree->SetBranchAddress("Ef",&E);
  tree->SetBranchAddress("resid",&resid);
  Double_t kp, epi; //proton kinetic energy, pion energy

  //Initialise the required variables over a for loop with no energy cut
  for (Long64_t i = 0; i < elist->GetN(); i++){

    tree->GetEntry(elist->GetEntry(i)); //Get the required entry

    //Get the outgoing electron kinetic energy
    double_t ke = el - ELECTRON_MASS; 

    //Loop through the final state hadrons
    for (Long64_t k = 0; k < nf; k++){

      //Initialise hadronic kinematic variables
      // Double_t kp, epi; //proton kinetic energy, pion energy

      if (pdg[k] == 2212){ //If it is a proton
	kp = E[k] - PROTON_MASS; //Set the proton kinetic energy
      }

      if (pdg[k] == 211 || pdg[k] == -211){
	epi = E[k];
      }
    }

    //Calculate the reconstructed electron energy
    Double_t Ee = ke + kp + epi + C12_bind_en;  //Scattering off C-12 nuclei 

    //Add it to histograms depending on interaction type:
    if (qel){
      hcal_qel->Fill(Ee);
    }

    else if (mec){
      hcal_mec->Fill(Ee);
    }

    else if (dis){
      hcal_dis->Fill(Ee);
    }

    else if (res && resid == 0){
      hcal_dres->Fill(Ee);
    }

    else if (res && resid != 0){
      hcal_res->Fill(Ee);
    }

    else{}
  }

  //Initialise the required variables over a for loop with energy cut
  for (Long64_t i = 0; i < Eelist->GetN(); i++){

    tree->GetEntry(Eelist->GetEntry(i)); //Get the required entry

    //Get the outgoing electron kinetic energy
    double_t ke = el - ELECTRON_MASS; 

    //Loop through the final state hadrons
    for (Long64_t k = 0; k < nf; k++){

      //Initialise hadronic kinematic variables
      // Double_t kp, epi; //proton kinetic energy, pion energy

      if (pdg[k] == 2212){ //If it is a proton
	kp = E[k] - PROTON_MASS; //Set the proton kinetic energy
      }

      if (pdg[k] == 211 || pdg[k] == -211){
	epi = E[k];
      }
    }

    //Calculate the reconstructed electron energy
    Double_t Ee = ke + kp + epi + C12_bind_en;  //Scattering off C-12 nuclei 

    //Add it to histograms depending on interaction type:
    if (qel){
      hcal_Eqel->Fill(Ee);
    }

    else if (mec){
      hcal_Emec->Fill(Ee);
    }

    else if (dis){
      hcal_Edis->Fill(Ee);
    }

    else if (res && resid == 0){
      hcal_Edres->Fill(Ee);
    }

    else if (res && resid != 0){
      hcal_Eres->Fill(Ee);
    }

    else{}
  }

  //Stack the background histograms
  THStack* hcal_background = new THStack("hs_background", "Reconstructed Electron Energy for {}^{12}C 2.2GeV Electron #Delta Resonances; E (GeV); Events / 50 MeV");
  THStack* hcal_Ebackground = new THStack("hs_Ebackground", "Reconstructed Electron Energy for {}^{12}C 2.2GeV Electron #Delta Resonances; E (GeV) With Cuts; Events / 50 MeV");

  hcal_background->Add(hcal_mec);
  hcal_mec->SetFillColorAlpha(kOrange,0.2);
  hcal_background->Add(hcal_qel);
  hcal_qel->SetFillColorAlpha(kPink,0.2);
  hcal_background->Add(hcal_res);
  hcal_res->SetFillColorAlpha(kGreen,0.2);
  hcal_background->Add(hcal_dis);
  hcal_dis->SetFillColorAlpha(kBlue,0.2);

  hcal_Ebackground->Add(hcal_Emec);
  hcal_Emec->SetFillColorAlpha(kOrange,0.2);
  hcal_Ebackground->Add(hcal_Eqel);
  hcal_Eqel->SetFillColorAlpha(kPink,0.2);
  hcal_Ebackground->Add(hcal_Eres);
  hcal_Eres->SetFillColorAlpha(kGreen,0.2);
  hcal_Ebackground->Add(hcal_Edis);
  hcal_Edis->SetFillColorAlpha(kBlue,0.2);

  //Create canvas for the regular histograms
  TCanvas* c_cal = new TCanvas("c_cal", "stacked hists", 10, 10, 1920, 1080);
  hcal_dres->SetLineColor(kRed);
  hcal_dres->SetLineWidth(4);
  hcal_dres->SetTitle("Reconstructed Electron Energy for {}^{12}C 2.2GeV Electron #Delta Resonances; E (GeV); Events / 50 MeV");
  hcal_dres->SetMaximum(250);
  hcal_dres->Draw("HIST");
  hcal_background->Draw("HIST SAME");

  //Create Legend
  auto hcal_legend = new TLegend(0.68, 0.68, 0.9, 0.9);
  hcal_legend->SetHeader("Events", "C");
  hcal_legend->AddEntry(hcal_dres, Form("#Delta 1232 Resonances (%lld)", tree->GetEntries(resd_cut.c_str())), "l");
  hcal_legend->AddEntry(hcal_qel, Form("QE (%lld)", tree->GetEntries(qel_cut.c_str())), "f");
  hcal_legend->AddEntry(hcal_res, Form("All other RES (%lld)", tree->GetEntries(res_cut.c_str())), "f");
  hcal_legend->AddEntry(hcal_mec, Form("MEC (%lld)", tree->GetEntries(mec_cut.c_str())), "f");
  hcal_legend->AddEntry(hcal_dis, Form("DIS (%lld)", tree->GetEntries(dis_cut.c_str())), "f");
  hcal_legend->Draw();

  //Save canvas
  c_cal->SaveAs("Calorimetric_allW.png");

  //Create canvas for cut histogram
  TCanvas* c_Ecal = new TCanvas("c_Ecal", "stacked hists", 10, 10, 1920, 1080);
  hcal_Edres->SetLineColor(kRed);
  hcal_Edres->SetLineWidth(4);
  hcal_Edres->SetTitle("Reconstructed Electron Energy for {}^{12}C 2.2GeV Electron #Delta Resonances With Energy Cut; E (GeV); Events / 50 MeV");
  hcal_Edres->SetMaximum(250);
  hcal_Edres->Draw("HIST");
  hcal_Ebackground->Draw("HIST SAME");


  auto hcal_Elegend = new TLegend(0.68, 0.68, 0.9, 0.9);
  hcal_Elegend->SetHeader("Events", "C");
  hcal_Elegend->AddEntry(hcal_Edres, Form("#Delta 1232 Resonances (%lld)", tree->GetEntries(resd_Ecut.c_str())), "l");
  hcal_Elegend->AddEntry(hcal_Eqel, Form("QE (%lld)", tree->GetEntries(qel_Ecut.c_str())), "f");
  hcal_Elegend->AddEntry(hcal_Eres, Form("All other RES (%lld)", tree->GetEntries(res_Ecut.c_str())), "f");
  hcal_Elegend->AddEntry(hcal_Emec, Form("MEC (%lld)", tree->GetEntries(mec_Ecut.c_str())), "f");
  hcal_Elegend->AddEntry(hcal_Edis, Form("DIS (%lld)", tree->GetEntries(dis_Ecut.c_str())), "f");
  hcal_Elegend->Draw();

  c_Ecal->SaveAs("Calorimetric_cutW.png");
}

int Reconstruction(){
  TTree* t1 = tree();
  float Woptimal = PlotOptimalCut(t1);
  draw(t1, Woptimal);
  return 0;
}

