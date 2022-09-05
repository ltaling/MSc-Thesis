// This file takes in GENIE Monte Carlo (MC) data for the C-12
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

double weighted_average(TH1F* h1) {
	double num_sum = 0.0;
	double den_sum = 0.0;
	for (int i = 0; i < h1->GetNbinsX(); i++) {
		if (h1->GetBinContent(i) == 0) continue;
		num_sum += h1->GetBinCenter(i)*h1->GetBinContent(i);
		den_sum += h1->GetBinContent(i);
	}

	return num_sum/den_sum;
}

void prot_trans_vs_prot_mom_targetC12_beamEnergy2261_exclusivePlus(){
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

	// These correction factors are obtained from looking at the number of struck protons and neutrons in each sample
	// In SuSAv2, the number of protons and neutrons was incorrectly set to be equal
	double corrRange[3] = {4900./(4900. + 1800.), 1830./(1830. + 830.), 395./(395. + 205.)};

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
	std::string mc_file_path ("/genie/app/users/nsteinbe/e4nu/CLAS/GENIE/Second_try/");
	std::string data_file_path ("/genie/app/users/nsteinbe/e4nu/CLAS/DATA/Second_try/");

	TFile *input_mc_incl;
	TFile *input_mc_excl;
	TFile *input_data_incl;
	TFile *input_data_excl;

	std::string mc_incl[N_RANGES];
        std::string mc_excl[N_RANGES];

        std::string data_incl[N_RANGES];
        std::string data_excl[N_RANGES];

        mc_incl[0] = mc_file_path + "FullSample_Range1_Genie_1_C12_2.261000.root";
        mc_incl[1] = mc_file_path + "FullSample_Range2_Genie_1_C12_2.261000.root";
        mc_incl[2] = mc_file_path + "FullSample_Range3_Genie_1_C12_2.261000.root";
        mc_excl[0] = mc_file_path + "Excl_Plus_Range1_Genie_1_C12_2.261000.root";
        mc_excl[1] = mc_file_path + "Excl_Plus_Range2_Genie_1_C12_2.261000.root";
        mc_excl[2] = mc_file_path + "Excl_Plus_Range3_Genie_1_C12_2.261000.root";

        data_incl[0] = data_file_path + "FullSample_Range1_Data__C12_2.261000.root";
        data_incl[1] = data_file_path + "FullSample_Range2_Data__C12_2.261000.root";
        data_incl[2] = data_file_path + "FullSample_Range3_Data__C12_2.261000.root";
        data_excl[0] = data_file_path + "Excl_Plus_Range1_Data__C12_2.261000.root";
        data_excl[1] = data_file_path + "Excl_Plus_Range2_Data__C12_2.261000.root";
        data_excl[2] = data_file_path + "Excl_Plus_Range3_Data__C12_2.261000.root";


        //exclusive histograms: proton momenta
        TH1F* h1_prot_mom_sectors_interactions_mc[N_RANGES][N_SECTORS][N_INT];
        
        //TH1F* h1_prot_mom_ranges_mc[N_RANGES];
        //TH1F* h1_prot_mom_sectors_interactions_mc[N_SECTORS][N_INT]; 
	

	TH1F* h1_prot_mom_sectors_data[N_RANGES][N_SECTORS]; 


	//Inclusive histograms: electron momenta
	TH1F* h1_el_mom_sectors_interactions_mc[N_RANGES][N_SECTORS][N_INT];
	//TH1F* h1_el_mom_ranges_mc[N_RANGES];
	//TH1F* h1_el_mom_sectors_interactions_mc[N_SECTORS][N_INT];

	TH1F *h1_el_mom_sectors_data[N_RANGES][N_SECTORS];


	// loop over all files and extract histograms and add them up
	for (int file_i = 0; file_i < N_RANGES; file_i++){ 
		input_mc_incl = TFile::Open( TString::Format("%s", mc_incl[file_i].c_str()));
		input_data_incl = TFile::Open( TString::Format("%s", data_incl[file_i].c_str()));
		input_mc_excl = TFile::Open( TString::Format("%s", mc_excl[file_i].c_str()));
		input_data_excl = TFile::Open( TString::Format("%s", data_excl[file_i].c_str()));
	
		// loop over sectors
		for (int sector = 0; sector < N_SECTORS; sector++) {

			h1_prot_mom_sectors_data[file_i][sector] = (TH1F*)input_data_excl->Get( TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", sector));

			h1_el_mom_sectors_data[file_i][sector] = (TH1F*)input_data_incl->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__%i", sector));
			//h1_el_mom_sectors_data[file_i][sector] = (TH1F*)input_data_incl->Get(TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", sector));


			// loop over interactions for MC
			for (int interaction = 0; interaction < N_INT; interaction++) {
				h1_prot_mom_sectors_interactions_mc[file_i][sector][interaction] = (TH1F*)input_mc_excl->Get( TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", interaction+1, sector));

				h1_el_mom_sectors_interactions_mc[file_i][sector][interaction] = (TH1F*)input_mc_incl->Get(TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__%i", interaction+1, sector));

				//h1_el_mom_sectors_interactions_mc[file_i][sector][interaction] = (TH1F*)input_mc_incl->Get(TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", interaction+1, sector));


				if(sector == 0 && interaction == 0) continue;
				else { 
					h1_prot_mom_sectors_interactions_mc[file_i][0][0]->Add(h1_prot_mom_sectors_interactions_mc[file_i][sector][interaction]);
					h1_el_mom_sectors_interactions_mc[file_i][0][0]->Add(h1_el_mom_sectors_interactions_mc[file_i][sector][interaction]);
				}
			}

			if(sector == 0) continue;
			else {
				h1_prot_mom_sectors_data[file_i][0]->Add(h1_prot_mom_sectors_data[file_i][sector]);
				h1_el_mom_sectors_data[file_i][0]->Add(h1_el_mom_sectors_data[file_i][sector]);
			}
		}		
	}


	
	double mc_conglomerate_trans[N_RANGES]; // transparency values for each range summed up, MC
        double data_conglomerate_trans[N_RANGES]; // transparency values for each range summed up, data

	double mc_conglomerate_trans_err[N_RANGES]; // transparency values for each range summed up, MC
        double data_conglomerate_trans_err[N_RANGES]; // transparency values for each range summed up, data

	double mc_avg_prot_mom[N_RANGES]; // average of the proton momentum to be used as the x coordinate for plotting conglomerate transparencies, MC
	double data_avg_prot_mom[N_RANGES]; // average of the proton momentum to be used as the x coordinate for plotting conglomerate transparencies, data


	for (int idx = 0; idx < N_RANGES; idx++) { 
	        h1_prot_mom_sectors_interactions_mc[idx][0][0]->Rebin( rebin_num); 
	        h1_prot_mom_sectors_interactions_mc[idx][0][0]->Sumw2();

	        h1_el_mom_sectors_interactions_mc[idx][0][0]->Rebin(rebin_num);
		h1_el_mom_sectors_interactions_mc[idx][0][0]->Sumw2();

	        h1_prot_mom_sectors_data[idx][0]->Rebin( rebin_num); 
	        h1_prot_mom_sectors_data[idx][0]->Sumw2();

	        h1_el_mom_sectors_data[idx][0]->Rebin(rebin_num);
	        h1_el_mom_sectors_data[idx][0]->Sumw2();

	        mc_avg_prot_mom[idx] = weighted_average(h1_prot_mom_sectors_interactions_mc[idx][0][0]);
	        data_avg_prot_mom[idx] = weighted_average(h1_prot_mom_sectors_data[idx][0]);

	        int nbins_mc_num = h1_prot_mom_sectors_interactions_mc[idx][0][0]->GetNbinsX();
	        int nbins_data_num = h1_prot_mom_sectors_data[idx][0]->GetNbinsX();
	        int nbins_mc_denom = h1_el_mom_sectors_interactions_mc[idx][0][0]->GetNbinsX();
	        int nbins_data_denom = h1_el_mom_sectors_data[idx][0]->GetNbinsX();

	        double mc_num_err;
        	double mc_denom_err;
        	double data_num_err;
        	double data_denom_err;

	        // calculate the proton transparency for each range by integrating over exclusive and inclusive distributions and taking the ratios
	        double mc_num = h1_prot_mom_sectors_interactions_mc[idx][0][0]->IntegralAndError(1,nbins_mc_num,mc_num_err);
	        std::cout << "mc num: " << mc_num << "+/- " << mc_num_err << "\n";
	        double mc_denom = h1_el_mom_sectors_interactions_mc[idx][0][0]->IntegralAndError(1,nbins_mc_denom,mc_denom_err);
	        std::cout << "mc denom: " << mc_denom << "+/- " << mc_denom_err << "\n";
	        double data_num = h1_prot_mom_sectors_data[idx][0]->IntegralAndError(1,nbins_data_num,data_num_err);
	        std::cout << "data num: " << data_num << "+/- " << data_num_err << "\n";
	        double data_denom = h1_el_mom_sectors_data[idx][0]->IntegralAndError(1,nbins_data_denom,data_denom_err);
	        std::cout << "data denom: " << data_denom << "+/- " << data_denom_err << "\n";
		data_denom *= corrRange[idx];

		mc_conglomerate_trans[idx] = mc_num/mc_denom;


		data_conglomerate_trans[idx] = data_num/data_denom;

		mc_conglomerate_trans_err[idx] = mc_conglomerate_trans[idx]*sqrt(pow(mc_num_err/mc_num,2) + pow(mc_denom_err/mc_denom,2));

		data_conglomerate_trans_err[idx] = data_conglomerate_trans[idx]*sqrt(pow(data_num_err/data_num,2) + pow(data_denom_err/data_denom,2));

		//mc_conglomerate_trans[idx] *= corrRange[idx];

		std::cout << "MC" << "\n";
		std::cout << "P proton = " << mc_avg_prot_mom[idx]<<", T =  " << mc_conglomerate_trans[idx] << "+/- " << mc_conglomerate_trans_err[idx] << "\n";
		std::cout << "Data" << "\n";
		std::cout << "P proton = " << data_avg_prot_mom[idx]<<", T =  " << data_conglomerate_trans[idx] << "+/- " << data_conglomerate_trans_err[idx]<< "\n";

	}

	// CONGLOMERATE HISTOGRAM CODE BELOW --------------------------------------------
	TGraph *g1_mc_conglomerate_trans = new TGraph( N_RANGES, mc_avg_prot_mom, mc_conglomerate_trans);
	TGraph *g1_data_conglomerate_trans = new TGraph( N_RANGES, data_avg_prot_mom, data_conglomerate_trans);

	g1_mc_conglomerate_trans->SetMarkerStyle( kFullSquare);
	g1_mc_conglomerate_trans->SetMarkerColor( kGreen+2);
	g1_mc_conglomerate_trans->SetLineColor( kBlack);
	g1_mc_conglomerate_trans->GetXaxis()->SetRangeUser( xaxis_lb-0.1, xaxis_ub+0.1); // give the x-axis bounds some room
	g1_mc_conglomerate_trans->GetYaxis()->SetRangeUser( yaxis_lb, yaxis_ub);

	TCanvas* c2;
        c2 = new TCanvas( TString::Format("c2"), TString::Format("c2"), 800, 600);
        c2->cd();
        g1_mc_conglomerate_trans->Draw();
        g1_data_conglomerate_trans->Draw("same");

        TFile *output = TFile::Open("../output/trans_vs_prot_mom/C12_2261.root", "RECREATE");
        output->cd();
        g1_data_conglomerate_trans->Write("data_graph");
        g1_mc_conglomerate_trans->Write("genie_graph");
        //c2->SaveAs(TString::Format("../output/trans_vs_prot_mom/trans_vs_prot_mom_targetC12_beamEnergy2261_exclusivePlus_conglom.pdf"));
}
