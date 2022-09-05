// This file takes in GENIE Monte Carlo (MC) data for the Fe-56
// targets with a 2.261 GeV beam energy and calculates the proton 
// transparency (ratio of exclusive cut case events to inclusive cut case 
// events) as a function of proton momentum. Note this script take 
// the exclusive+ case as the numerator in the transparency ratio. 
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

#include <iostream>
#include <vector>
#include <stdio.h>


void prot_trans_vs_prot_mom_target56Fe_beamEnergy2261_exclusivePlus(){
	// You will need to change these variables if you change the amount of data you run
	// over or if you want to more than just GENIE MC and CLAS data.
	const int N_RANGES = 3; // number of Ranges; see similarities in thetaEl arguments in file_names below
	const int N_FILES = N_RANGES * 2 * 2; // number of files to process; 
					      // first factor of two because there is a numerator and denominator for the transparency ratios
					      // second factor of two because we consider both GENIE MC and CLAS data

        // You will need to change these if you want to change how things are plotted.
        // NOTE: be sure to add new color_options if you consider new targets
        // NOTE: be sure to add new draw_options if you consider new data
	const int color_options[N_RANGES] = {kBlue, kBlack, kRed}; // colors of each Range plotted
        const char* mc_draw_options[N_RANGES] = {"HIST", "HIST SAME", "HIST SAME"};
        const char* data_draw_options[N_RANGES] = {"E1 SAME", "E1 SAME", "E1 SAME"};
	const double xaxis_lb = 0.0; // x-axis lower bound for plotting
        const double xaxis_ub = 1.6; // x-axis upper bound for plotting
        const double yaxis_lb = 0.0; // y-axis lower bound for plotting
        const double yaxis_ub = 1.5; // y-axis upper bound for plotting
        const int rebin_num = 50; // how many bins you want to put into each one of your MC bins in the final processed plot
				     // NOTE: the histograms I use here are made in the genie_analysis.C script in the main
				     //       e4nu directory, and at the time of writing this those histograms ranged from
				     //       0 GeV to 6 GeV and had a total of 6000 bins (1000 bins per GeV)
	const double leg_x1 = 0.2; // x-axis lower bound for the legend
	const double leg_x2 = 0.4; // x-axis upper bound for the legend
	const double leg_y1 = 0.65; // y-axis lower bound for the legend
	const double leg_y2 = 0.875; // y-axis upper bound for the legend
        gStyle->SetErrorX(0);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Users don't have to read all the code below this dashed line unless they want to do a deep dive into the code, or if they want to debug, etc.

        // You probably won't be changing N_SECTORS unless you use an experiment other than
        // CLAS, which will most likely make this script obsolete anyway. Since I build all
        // of my distributions from histograms that have everything separated into different
        
        const int N_SECTORS = 6; // number of sectors
        const int N_INT = 4; // number of interactions

        // create a canvas on which to draw the histograms
        TCanvas* c;
        c = new TCanvas( TString::Format("c"), TString::Format("c"), 800, 600);

	// these are the file paths where the GENIE MC files and CLAS data files will be
	std::string mc_file_path ("/pnfs/genie/persistent/users/smithja/e4nu_GENIE_simulation_output_files/56Fe_2261/");
	std::string data_file_path ("/pnfs/genie/persistent/users/smithja/e4nu_CLAS_data_output_files/56Fe_2261/");

	TFile *input;
	std::string file_names[ N_FILES];
	file_names[0] =	mc_file_path + "genie_target56Fe_beamEnergy2261_50000000Events_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_MottXSecEq1.root";
	file_names[1] =	mc_file_path + "genie_target56Fe_beamEnergy2261_50000000Events_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_protSectors453_deltaPhiProtEq45pt00_thetaProt47pt00to70pt00_protMomLBEq0pt60_MottXSecEq1.root";
	file_names[2] = mc_file_path + "genie_target56Fe_beamEnergy2261_50000000Events_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_MottXSecEq1.root";
        file_names[3] = mc_file_path + "genie_target56Fe_beamEnergy2261_50000000Events_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt40pt00to55pt00_protMomLBEq1pt00_MottXSecEq1.root";
        file_names[4] = mc_file_path + "genie_target56Fe_beamEnergy2261_50000000Events_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_MottXSecEq1.root";
        file_names[5] = mc_file_path + "genie_target56Fe_beamEnergy2261_50000000Events_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to45pt00_protMomLBEq1pt30_MottXSecEq1.root"; 
        file_names[6] = data_file_path + "data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90.root";
        file_names[7] = data_file_path + "data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl20pt00to25pt00_elMomLBEq1pt90_protSectors453_deltaPhiProtEq45pt00_thetaProt47pt00to70pt00_protMomLBEq0pt60.root";
        file_names[8] = data_file_path + "data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65.root";
        file_names[9] = data_file_path + "data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl30pt00to35pt00_elMomLBEq1pt65_protSectors453_deltaPhiProtEq45pt00_thetaProt40pt00to55pt00_protMomLBEq1pt00.root";
        file_names[10] = data_file_path + "data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40.root";
        file_names[11] = data_file_path + "data_target56Fe_beamEnergy2261_elSectors126_deltaPhiElEq12pt00_thetaEl40pt00to45pt00_elMomLBEq1pt40_protSectors453_deltaPhiProtEq45pt00_thetaProt35pt00to45pt00_protMomLBEq1pt30.root";	

        TH1F* h1_prot_mom_sectors_interactions_mc[N_FILES/2][N_SECTORS][N_INT]; // array of TH1Fs to store all of the unprocessed MC histograms
	TH1F* h1_prot_mom_sectors_data[N_FILES/2][N_SECTORS]; // array of TH1Fs to store all of the unprocessed data histograms
	// read and format the data files
	for (int file_i = 0; file_i < N_FILES; file_i++){ // for all of the files in file_names[] 
		input = TFile::Open( TString::Format("%s", file_names[file_i].c_str()));
	
		// extract the ROOT histograms for each sector/interaction from the ith file
		for (int sector = 0; sector < N_SECTORS; sector++) {
			if (file_i < (N_FILES / 2)) { // MC has its data separated into interactions
				for (int interaction = 0; interaction < N_INT; interaction++) {
					h1_prot_mom_sectors_interactions_mc[file_i][sector][interaction] = (TH1F*)input->Get( TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", interaction+1, sector));
				}
			}
			else { // CLAS data, on the other hand, cannot separate by interaction and thus has its data stored all in the interaction = 0 file
				h1_prot_mom_sectors_data[file_i%(N_FILES/2)][sector] = (TH1F*)input->Get( TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", sector));
			}
		}

		// compile the sector/interaction plots
		for (int sector = 0; sector < N_SECTORS; sector++) {
			if (file_i < (N_FILES/2)) { // MC has its data separated into interactions
				for (int interaction = 0; interaction < N_INT; interaction++){
					if (sector == 0 && interaction == 0) { continue; } // skip the first iteration because we are adding everything to that histogram
					h1_prot_mom_sectors_interactions_mc[file_i][0][0]->Add( h1_prot_mom_sectors_interactions_mc[file_i][sector][interaction]); 
//					h1_prot_mom_sectors_interactions_mc[file_i][0][0]->Sumw2(); // update the errors as we compile; this should produce a lot of warning statements
				}
			}
			else { // CLAS data, on the other hand, cannot so we do not need to compile based on interactions
				if (sector == 0) { continue; } // skip the first sector because we are adding everything to that histogram
				h1_prot_mom_sectors_data[file_i%(N_FILES/2)][0]->Add( h1_prot_mom_sectors_data[file_i%(N_FILES/2)][sector]);
//				h1_prot_mom_sectors_data[file_i%(N_FILES/2)][0]->Sumw2(); // update the errors as we compile; this should produce a lot of warning statements
			}
		} 
	}
	Double_t mc_conglomerate_trans[N_RANGES]; // transparency values for each range summed up, MC
	Double_t mc_conglomerate_trans_err[N_RANGES];
        Double_t data_conglomerate_trans[N_RANGES]; // transparency values for each range summed up, data
        Double_t data_conglomerate_trans_err[N_RANGES];
	Double_t mc_avg_prot_mom[N_RANGES]; // average of the proton momentum to be used as the x coordinate for plotting conglomerate transparencies, MC
	Double_t data_avg_prot_mom[N_RANGES]; // average of the proton momentum to be used as the x coordinate for plotting conglomerate transparencies, data
	
	for (int idx = 0; idx < N_RANGES; idx++) { 

		// format the data in the histograms so that both have the same axis bounds and number of bins
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetXaxis()->SetRangeUser( xaxis_lb, xaxis_ub);
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->Rebin( rebin_num); 
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->Sumw2();
		h1_prot_mom_sectors_interactions_mc[2*idx][0][0]->GetXaxis()->SetRangeUser( xaxis_lb, xaxis_ub);
                h1_prot_mom_sectors_interactions_mc[2*idx][0][0]->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);
                h1_prot_mom_sectors_interactions_mc[2*idx][0][0]->Rebin( rebin_num); 
		h1_prot_mom_sectors_interactions_mc[2*idx][0][0]->Sumw2();
		
		h1_prot_mom_sectors_data[2*idx+1][0]->GetXaxis()->SetRangeUser( xaxis_lb, xaxis_ub);
                h1_prot_mom_sectors_data[2*idx+1][0]->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);
                h1_prot_mom_sectors_data[2*idx+1][0]->Rebin( rebin_num); 
                h1_prot_mom_sectors_data[2*idx+1][0]->Sumw2();
                h1_prot_mom_sectors_data[2*idx][0]->GetXaxis()->SetRangeUser( xaxis_lb, xaxis_ub);
                h1_prot_mom_sectors_data[2*idx][0]->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);
                h1_prot_mom_sectors_data[2*idx][0]->Rebin( rebin_num); 
                h1_prot_mom_sectors_data[2*idx][0]->Sumw2();
		
		// calculate the proton transparency histograms from the respective proton momentum histograms
		h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->Divide( h1_prot_mom_sectors_interactions_mc[2*idx][0][0]);
		h1_prot_mom_sectors_data[2*idx+1][0]->Divide( h1_prot_mom_sectors_data[2*idx][0]);

		// calculate the proton transparency for each range as a conglomerate
		mc_conglomerate_trans[idx] = h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->IntegralAndError(0,h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetNbinsX(), mc_conglomerate_trans_err[idx]); 
		//std::cout << "Nbins = " << h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetNbinsX() << "\n";

		data_conglomerate_trans[idx] = h1_prot_mom_sectors_data[2*idx+1][0]->IntegralAndError(0,h1_prot_mom_sectors_data[2*idx+1][0]->GetNbinsX(), data_conglomerate_trans_err[idx]);

		Int_t mc_lowbin = h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->FindFirstBinAbove(0,1);
		Int_t mc_highbin = h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->FindLastBinAbove(0,1);

		Int_t data_lowbin = h1_prot_mom_sectors_data[2*idx+1][0]->FindFirstBinAbove(0,1);
		Int_t data_highbin = h1_prot_mom_sectors_data[2*idx+1][0]->FindLastBinAbove(0,1);

		int dataN = data_highbin - data_lowbin + 1;
		int mcN = mc_highbin - mc_lowbin +1;

		std::cout << "dataN = " << dataN << "\n";
		std::cout << "mcN = " << mcN << "\n";

		mc_conglomerate_trans[idx] /= mcN;
		mc_conglomerate_trans_err[idx] /= mcN;

		data_conglomerate_trans[idx] /= dataN;
		data_conglomerate_trans_err[idx] /= dataN;

		std::cout << "mc_conglomerate_trans[idx]: " << mc_conglomerate_trans[idx] << " +/- " << mc_conglomerate_trans_err[idx] <<  ", data_conglomerate_trans[idx]: " << data_conglomerate_trans[idx] << " +/- " << data_conglomerate_trans_err[idx] << "\n";

		// Calculate the mean of each range to use as the x coordinate when plotting conglomerate transparencies.
                // Note that this is done after we have created the transparency plots as a function of proton momentum
                // because we would otherwise have two average proton momenta, one for the exclusive+ (numerator) case and
                // one for the inclusive (denominator) case. Rather than trying to reconcile those two averages, we find
                // the mean of the transparency v. proton momentum histograms.
                // If you want to change this proton momentum mean so that you consider the exclusive+ and inclusive
                // histograms separately, make sure you calculate the mean values before doing Divide().
                mc_avg_prot_mom[idx] = h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetMean();
		data_avg_prot_mom[idx] = h1_prot_mom_sectors_data[2*idx+1][0]->GetMean();
		
		// format the final histogram to make it look pretty
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->SetLineColor( color_options[idx]);
		h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetXaxis()->SetRangeUser( xaxis_lb-0.1, xaxis_ub+0.1); // give the x-axis bounds some room
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->SetStats( 0);
		h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->Sumw2();
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetXaxis()->SetTitle("Proton Momentum [GeV/c]");
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetYaxis()->SetTitle("Proton Transparency (Exclusive+ / Inclusive)");
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetXaxis()->CenterTitle( true);
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->GetYaxis()->CenterTitle( true);
                h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->SetTitle("Fe-56 with 2.261 GeV Beam Energy (GENIE & CLAS)");

                h1_prot_mom_sectors_data[2*idx+1][0]->SetMarkerColor( color_options[idx]);
                h1_prot_mom_sectors_data[2*idx+1][0]->SetMarkerStyle( kFullCircle);
                h1_prot_mom_sectors_data[2*idx+1][0]->SetLineColor( kBlack);
                h1_prot_mom_sectors_data[2*idx+1][0]->GetXaxis()->SetRangeUser( xaxis_lb-0.1, xaxis_ub+0.1); // give the x-axis bounds some room
                h1_prot_mom_sectors_data[2*idx+1][0]->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);
                h1_prot_mom_sectors_data[2*idx+1][0]->SetStats( 0);
                h1_prot_mom_sectors_data[2*idx+1][0]->Sumw2(); 
                h1_prot_mom_sectors_data[2*idx+1][0]->GetXaxis()->SetTitle("Proton Momentum [GeV/c]");
                h1_prot_mom_sectors_data[2*idx+1][0]->GetYaxis()->SetTitle("Proton Transparency (Exclusive+ / Inclusive)");
                h1_prot_mom_sectors_data[2*idx+1][0]->GetXaxis()->CenterTitle( true);
                h1_prot_mom_sectors_data[2*idx+1][0]->GetYaxis()->CenterTitle( true);
                h1_prot_mom_sectors_data[2*idx+1][0]->SetTitle("Fe-56 with 2.261 GeV Beam Energy (GENIE & CLAS)");
	}
	
	// draw the histograms
	for (int idx = 0; idx < N_RANGES; idx++) {
		h1_prot_mom_sectors_interactions_mc[2*idx+1][0][0]->Draw( mc_draw_options[idx]);
                h1_prot_mom_sectors_data[2*idx+1][0]->Draw( data_draw_options[idx]);
	}	

	// draw a legend for our histogram
        TLegend *legend = new TLegend( leg_x1, leg_y1, leg_x2, leg_y2);
        legend->AddEntry( h1_prot_mom_sectors_interactions_mc[1][0][0], "GENIE, Range 1");
        legend->AddEntry( h1_prot_mom_sectors_interactions_mc[3][0][0], "GENIE, Range 2");
        legend->AddEntry( h1_prot_mom_sectors_interactions_mc[5][0][0], "GENIE, Range 3");
        legend->AddEntry( h1_prot_mom_sectors_data[1][0], "CLAS, Range 1");
        legend->AddEntry( h1_prot_mom_sectors_data[3][0], "CLAS, Range 2");
        legend->AddEntry( h1_prot_mom_sectors_data[5][0], "CLAS, Range 3");
        legend->SetBorderSize( 0);
        legend->SetFillStyle( 0);
	legend->Draw();
	
	// save the histogram
	c->SetLeftMargin( 0.15);
	c->SetBottomMargin( 0.15);
	c->Update();
	c->SaveAs(TString::Format("../output/trans_vs_prot_mom/trans_vs_prot_mom_target56Fe_beamEnergy2261_exclusivePlus.pdf"));

	// CONGLOMERATE HISTOGRAM CODE BELOW --------------------------------------------
	TGraph *g1_mc_conglomerate_trans = new TGraph( N_RANGES, mc_avg_prot_mom, mc_conglomerate_trans);
	TGraph *g1_data_conglomerate_trans = new TGraph( N_RANGES, data_avg_prot_mom, data_conglomerate_trans);

	g1_mc_conglomerate_trans->SetMarkerStyle( kFullSquare);
	g1_mc_conglomerate_trans->SetMarkerColor( kGreen+2);
	//g1_mc_conglomerate_trans->SetLineColor( kBlack);
	g1_mc_conglomerate_trans->GetXaxis()->SetRangeUser( xaxis_lb-0.1, xaxis_ub+0.1); // give the x-axis bounds some room
	g1_mc_conglomerate_trans->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);
	g1_data_conglomerate_trans->SetMarkerStyle(kFullCircle);
	g1_data_conglomerate_trans->SetMarkerColor(kBlack);
	g1_data_conglomerate_trans->GetXaxis()->SetTitle("Proton momentum (GeV)");
	g1_data_conglomerate_trans->GetYaxis()->SetTitle("Transparency");
	g1_data_conglomerate_trans->SetTitle("2.261 GeV e- + 56Fe");

	TCanvas* c2;
        c2 = new TCanvas( TString::Format("c2"), TString::Format("c2"), 800, 600);
        c2->cd();
        g1_data_conglomerate_trans->Draw("ap");
        g1_mc_conglomerate_trans->Draw("p same");
        TFile *output = TFile::Open("../output/trans_vs_prot_mom/56Fe_2261.root", "RECREATE");
        output->cd();
        g1_data_conglomerate_trans->Write("data_graph");
        g1_mc_conglomerate_trans->Write("genie_graph");


        //c2->SaveAs(TString::Format("../output/trans_vs_prot_mom/trans_vs_prot_mom_target56Fe_beamEnergy2261_exclusivePlus_conglom.pdf"));

}
