// This file takes in GENIE Monte Carlo (MC) data and the CLAS experimental
// data for the He-4, C-12, and Fe-56 targets with a 2.261 GeV beam energy 
// and calculates the proton transparency (ratio of exclusive cut case 
// events to inclusive cut case events) as a function of target. 
//
// Make sure to edit the code not just at the top of this script. I have
// separated all of the y-axis data arrays into their respective ranges
// (as well as by GENIE MC and CLAS data) manually. I did not want to 
// fool around with how to project an N-dimensional array down to a 
// 1-dimensional array that is required for the TGraph constructor.
// Future users of this code may know of a better way to structure the
// data and how it is processed, and there may be a simple way to 
// project an N-dimensional array to a 1-dimensional array––like I said
// I did not try to even play around with this.
//
// Author: Jacob Smith (smithja)
// Date of Creation: 11/08/2021 

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
#include <TGraph.h> // smithja: added to plot the data arrays

#include <iostream>
#include <vector>
#include <stdio.h>

void prot_trans_vs_target_exclPlus(){
	const Int_t N_TARGETS = 3; // number of targets we process
	const Int_t N_RANGES = 3; // number of Ranges (see similarities in the thetaEl cuts)
	const Int_t N_FILES = N_TARGETS * N_RANGES * 2 * 2; // number of files to process; factors of two are because we consider two types of data and each transparency ratio comes from two runs, namely the inclusive run and the exclusive/exclusive+ run
        const Double_t nucleon_num[N_TARGETS] = {4, 12, 56} // x ordinates to be plotted
 
	const Int_t N_SECTORS = 6; // number of sectors
	const Int_t N_INT = 4; // number of interactions

	// this is the file path where all the GENIE MC files will be
	std::string mc_file_path ("/pnfs/genie/persistent/users/smithja/e4nu_GENIE_simulation_output_files/");
	std::string data_file_path ("/pnfs/genie/persistent/users/smithja/e4nu_CLAS_simulation_output_files/");

        TH1F* h1_prot_mom_sectors_interactions[N_FILES][N_SECTORS][N_INT]; // array of TH1Fs to store all of the unprocessed histograms
	
	TFile *input;
	std::string file_names[ N_FILES];
	// CAUTION: If you want to add more files, make sure to keep the file groupings below.
	//          If you change this ordering of the input files, you will have to propgate 
	//          the changes in the for loops, etc., below.
	// The grouping scheme is as follows:
	// 1 ) Split the files into halves with the MC data being first and the CLAS data being second
	// 2 ) Split those halves into sections according to the Range used (this is most clearly seen by the groupings in thetaEl); we used three
	// 3 ) For each of those Ranges, go through all of the targets (we had He-4, C-12, and Fe-56 at the time I wrote this)
	// 4 ) Each of those targets should have two simulations, the inclusive case and the exclusive+ case, and they should appear in that order
	file_names[0] =	mc_file_path + "4He_2261/genie_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_MottXSecEq1.root";
	file_names[1] =	mc_file_path + "4He_2261/genie_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_protSectors453_deltaPhiProtEq45pt00_thetaProt40pt00to80pt00_protMomLBEq0pt60_MottXSecEq1.root";
        file_names[6] = mc_file_path + "C12_2261/genie_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt80_MottXSecEq1.root";
        file_names[7] = mc_file_path + "C12_2261/genie_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt80_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to80pt00_protMomLBEq0pt60_MottXSecEq1.root";
        file_names[12] = mc_file_path + "56Fe_2261/genie_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_MottXSecEq1.root";
        file_names[13] = mc_file_path + "56Fe_2261/genie_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_protSectors453_deltaPhiProtEq45pt00_thetaProt47pt00to70pt00_protMomLBEq0pt60_MottXSecEq1.root";
        file_names[2] = mc_file_path + "4He_2261/genie_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_MottXSecEq1.root";
        file_names[3] = mc_file_path + "4He_2261/genie_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt38pt00to62pt00_protMomLBEq0pt90_MottXSecEq1.root";
	file_names[8] = mc_file_path + "C12_2261/genie_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_MottXSecEq1.root";
        file_names[9] = mc_file_path + "C12_2261/genie_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to65pt00_protMomLBEq0pt90_MottXSecEq1.root";
        file_names[14] = mc_file_path + "56Fe_2261/genie_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_MottXSecEq1.root";
        file_names[15] = mc_file_path + "56Fe_2261/genie_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt40pt00to55pt00_protMomLBEq1pt00_MottXSecEq1.root";
        file_names[4] = mc_file_path + "4He_2261/genie_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_MottXSecEq1.root";
        file_names[5] = mc_file_path + "4He_2261/genie_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt30pt00to50pt00_protMomLBEq1pt25_MottXSecEq1.root";
        file_names[10] = mc_file_path + "C12_2261/genie_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_MottXSecEq1.root";
        file_names[11] = mc_file_path + "C12_2261/genie_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt30pt00to50pt00_protMomLBEq1pt25_MottXSecEq1.root";
        file_names[16] = mc_file_path + "56Fe_2261/genie_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_MottXSecEq1.root";
        file_names[17] = mc_file_path + "56Fe_2261/genie_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to45pt00_protMomLBEq1pt30_MottXSecEq1.root";
	file_names[18] = data_file_path + "4He_2261/data_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90.root";
        file_names[19] = data_file_path + "4He_2261/data_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_protSectors453_deltaPhiProtEq45pt00_thetaProt40pt00to80pt00_protMomLBEq0pt60.root";
        file_names[24] = data_file_path + "C12_2261/data_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt80.root";
        file_names[25] = data_file_path + "C12_2261/data_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt80_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to80pt00_protMomLBEq0pt60.root"; 
        file_names[30] = data_file_path + "56Fe_2261/data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90.root";
        file_names[31] = data_file_path + "56Fe_2261/data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_protSectors453_deltaPhiProtEq45pt00_thetaProt47pt00to70pt00_protMomLBEq0pt60.root";
	file_names[20] = data_file_path + "4He_2261/data_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65.root";
        file_names[21] = data_file_path + "4He_2261/data_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt38pt00to62pt00_protMomLBEq0pt90.root";
        file_names[26] = data_file_path + "C12_2261/data_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65.root";
        file_names[27] = data_file_path + "C12_2261/data_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to65pt00_protMomLBEq0pt90.root"; 
        file_names[32] = data_file_path + "56Fe_2261/data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65.root";
        file_names[33] = data_file_path + "56Fe_2261/data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt40pt00to55pt00_protMomLBEq1pt00.root";
	file_names[22] = data_file_path + "4He_2261/data_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40.root";
        file_names[23] = data_file_path + "4He_2261/data_target4He_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt30pt00to50pt00_protMomLBEq1pt25.root";
        file_names[28] = data_file_path + "C12_2261/data_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40.root";
        file_names[29] = data_file_path + "C12_2261/data_targetC12_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt30pt00to50pt00_protMomLBEq1pt25.root";
        file_names[34] = data_file_path + "56Fe_2261/data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40.root";
        file_names[35] = data_file_path + "56Fe_2261/data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to45pt00_protMomLBEq1pt30.root";

	char sector_str [20];
	char interaction_str [20];
	// create the proton momentum distributions from the ROOT file data and store
	// those histograms in the h1_prot_momentum TH1F array
	for (int file_i = 0; file_i < N_FILES; file_i++){ // for all of the files in file_names[] 
		input = TFile::Open( TString::Format("%s", file_names[file_i].c_str()));
		
		// extract the ROOT histograms for each sector/interaction from the ith file
		for (int sector = 0; sector < N_SECTORS; sector++) {
			for (int interaction = 0; interaction < N_INT; interaction++) {
				sprintf(sector_str, "%i", sector);
				sprintf(interaction_str, "%i", interaction+1);
				h1_prot_mom_sectors_interactions[file_i][sector][interaction] = (TH1F*)input->Get( TString::Format("h1_%s_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%s", interaction_str, sector_str));
			}
		}

		// compile the sector/interaction plots
		for (int sector = 0; sector < N_SECTORS; sector++) {
			for (int interaction = 0; interaction < N_INT; interaction++){
				if (sector == 0 && interaction == 0) { continue; } // skip the first iteration because we are adding everything to that histogram
				h1_prot_mom_sectors_interactions[file_i][0][0]->Add( h1_prot_mom_sectors_interactions[file_i][sector][interaction]); 
			}
		} 
	}

	// NOTE: be sure to add any new data arrays here
	Double_t mc_trans_range1 [ N_TARGETS];
	Double_t mc_trans_range2 [ N_TARGETS]; 
	Double_t mc_trans_range3 [ N_TARGETS];
	Double_t data_trans_range1 [ N_TARGETS];
	Double_t data_trans_range2 [ N_TARGETS];
	Double_t data_trans_range3 [ N_TARGETS];

	Double_t incl_integ;
	Double_t exclPlus_integ;
	Double_t ratio;
	for (int idx = 0; idx < (N_FILES / 2); idx++){
		incl_integ = h1_prot_mom_sectors_interactions[2*idx][0][0]->Integral();
		exclPlus_integ = h1_prot_mom_sectors_interactions[2*idx+1][0][0]->Integral();
		ratio = exclPlus_integ / incl_integ;
		
		// NOTE: be sure to add more if statements to accomodate any new data
		if (idx < N_TARGETS) { mc_trans_range1[idx%N_TARGETS] = ratio;
		else if (idx >= N_TARGETS && idx < 2*N_TARGETS) { mc_trans_range2[idx%N_TARGETS] = ratio;
		else if (idx >= 2*N_TARGETS && idx < 3*N_TARGETS) { mc_trans_range3[idx%N_TARGETS] = ratio;
		else if (idx >= 3*N_TARGETS && idx < 4*N_TARGETS) { data_trans_range1[idx%N_TARGETS] = ratio;
		else if (idx >= 4*N_TARGETS && idx < 5*N_TARGETS) { data_trans_range2[idx%N_TARGETS] = ratio;
		else if (idx >= 5*N_TARGETS && idx < 6*N_TARGETS) { data_trans_range3[idx%N_TARGETS] = ratio;
	} // in retrospect this is probably not the most efficent and future-proof way I could have done this loop; somebody could improve on this loop if they changed the grouping structure for the file_names

	// convert the data arrays into TGraph objects
	// NOTE: make sure to add new TGraph objects for
	//       new data arrays
	TGraphErrors* mc_range1_graph = new TGraph( N_TARGETS, nucleon_num, mc_trans_range1); // no error arguments passed for MC
        TGraphErrors* mc_range2_graph = new TGraph( N_TARGETS, nucleon_num, mc_trans_range3);
        TGraphErrors* mc_range3_graph = new TGraph( N_TARGETS, nucleon_num, mc_trans_range3);
        TGraphErrors* data_range1_graph = new TGraph( N_TARGETS, nucleon_num, data_trans_range1, ); NEED TO SOMEHOW CALCULATE ERRORS FOR THESE PLOTS. EASIEST WILL PROBABLY BE SQRT( INTEG) FROM FOR LOOP ABOVE
        TGraphErrors* data_range2_graph = new TGraph( N_TARGETS, nucleon_num, data_trans_range3, );
        TGraphErrors* data_range3_graph = new TGraph( N_TARGETS, nucleon_num, data_trans_range3, );

        // put all of the TGraph objects into a storage 
        // array so we can process them eficently in a for loop
        // NOTE: make sure to add any TGraphs to these arrays
        // 	 so they are formatted/processed correctly
        TGraph* mc_graph_array[N_RANGES] = { mc_trans_range1, mc_trans_range2, mc_trans_range3};
	TGraphErrors* data_graph_array[N_RANGES] = { data_trans_range1, data_trans_range2, data_trans_range3}
	TGraph

	// format all the TGraph objects, both data and MC
        const int color_options[N_RANGES] = {kBlue, kBlack, kRed};	
	for (int idx = 0; idx < (N_RANGES); idx++){ 
		// format the data
                mc_graph_array[idx]->GetYaxis()->SetRangeUser( 0, 1.1);
		mc_graph_array[idx]->GetXaxis()->SetRangeUser( nucleon_num[0], nucleon_num[ (sizeof(nucleon_num) / sizeof(nucleon_num[0])) - 1]); // second arg is the last element in nucleon_num
		mc_graph_array[idx]->SetLineColor( color_options[idx]);
		mc_graph_array[idx]->SetStats( 0);
		mc_graph_array[idx]->GetXaxis()->SetTitle("Nucleon Number");
		mc_graph_array[idx]->GetYaxis()->SetTitle("Proton Transparency (Exclusive+ / Inclusive)");
		mc_graph_array[idx]->GetXaxis()->CenterTitle( true);
                mc_graph_array[idx]->GetYaxis()->CenterTitle( true);
		mc_graph_array[idx]->SetTitle("Proton Transparency vs. Nucleon Number (GENIE & CLAS)");

		data_graph_array[idx]->GetYaxis()->SetRangeUser( 0, 1.1);
                data_graph_array[idx]->GetXaxis()->SetRangeUser( nucleon_num[0], nucleon_num[ (sizeof(nucleon_num) / sizeof(nucleon_num[0])) - 1]); // second arg is the last element in nucleon_num
		data_graph_array[idx]->SetMarkerColor( color_options[idx]);
                data_graph_array[idx]->SetMarkerStyle( 4);
                data_graph_array[idx]->SetStats( 0);
		data_graph_array[idx]->GetXaxis()->SetTitle("Nucleon Number");
		data_graph_array[idx]->GetYaxis()->SetTitle("Proton Transparency (Exclusive+ / Inclusive)");
                data_graph_array[idx]->GetXaxis()->CenterTitle( true);
                data_graph_array[idx]->GetYaxis()->CenterTitle( true);
                data_graph_array[idx]->SetTitle("Proton Transparency vs. Nucleon Number (GENIE & CLAS)");
	}	
        gStyle->SetErrorX( 0);

    	// create a canvas on which to draw the histograms
    	TCanvas* c;
    	c = new TCanvas( TString::Format("c"), TString::Format("c"), 800, 600);

        // draw the histograms
        NEED TO REDO DRAWING OPTIONS SINCE WE'RE DEALING WITH TGraphErrors OBJECTS INSTEAD OF TH1 OBJECTS
        const char* mc_draw_options[3] = {"HIST", "HIST SAME", "HIST SAME"};
//        const char* mc_draw_options2[3] = {"E SAME", "E SAME", "E SAME"};
	const char* data_draw_options[3] = {};
//	const char* data_draw_options2[3] = {};
	for (int idx = 0; idx < (N_RANGES); idx++) {
		mc_graph_array[idx]->Draw( mc_draw_options[idx]);
		data_graph_array[idx]->Draw( data_draw_options[idx]);
	}

	NEED TO UPDATE THE LEGEND
	// create and draw a legend for our histograms
    	TLegend *legend = new TLegend( 0.20, 0.7, 0.325, 0.825);
    	legend->AddEntry( h1_prot_mom_sectors_interactions[1][0][0], "Range 1");
    	legend->AddEntry( h1_prot_mom_sectors_interactions[3][0][0], "Range 2");
    	legend->AddEntry( h1_prot_mom_sectors_interactions[5][0][0], "Range 3");
    	legend->SetBorderSize( 0);
    	legend->SetFillStyle( 0);
    	legend->Draw();
	
	// save the histogram
	c->SetLeftMargin( 0.15);
	c->SetBottomMargin( 0.15);
	c->Update();
	c->SaveAs(TString::Format("../output/trans_vs_target/trans_vs_target_beamEnergy2261.pdf"));
}
