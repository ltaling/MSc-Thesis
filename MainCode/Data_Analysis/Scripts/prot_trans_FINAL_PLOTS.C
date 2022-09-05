// Author: Jacob Smith (smithja)
// Date of Creation: 12/03/2021 

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

void prot_trans_FINAL_PLOTS(){
	// manually entered data
	//
	// I was in a pinch when I created this file,
	// so I externally calculated the transparencies
	// and the associated errors according to the
	// following formulas:
	//
	// Transparency := T = E / I, where E is the exclusive+
	// (more accurately known as the true QE) case and I is
	// the inclusive case. Note that the exclusive+/true QE
	// case can be swapped out for just the exclusive case
	// or something else entirely if desired.
	//
	// The error on transparency is just adding the relevant
	// errors and partial derivatives up in quadrature given
	// that T = E / I.
	//
	// I take the weighted average (weighted by the number
	// of inclusive events) of each target's three
	// electron theta ranges to obtain the arrays without 
	// nuclear target specification. Weighting by the
	// number of inclusive events makes it so that the
	// error on the final transparencies is of the same
	// form as each individual target's transparency with
	// E_{total} = E_{Range 1} + E_{Range 2} + E_{Range 3},
        // I_{total} = I_{Range 1} + I_{Range 2} + I_{Range 3},
        // and T_{total} = E_{total} / I_{total}.
        //
        // You can find this data as well as its explicit calculation in the presentation at the following link:
        // https://docs.google.com/presentation/d/1f9JfIsC4lVHNc8YdF0F5NkQ-yupqg623cDoX-jhXPjc/edit?usp=sharing
        // The relevant quantities were on slides 36, 71, and 72 at the time of writing this. Please contact me
        // at jacobsmith256@gmail.com if you cannot access this presentation or have questions.
	Double_t mc_trans_4He[3] = {0.797, 0.755, 0.702}; // GENIE transparencies for He-4
	Double_t data_trans_4He[3] = {0.825, 0.784, 0.725}; // CLAS transparencies for He-4
	Double_t mc_trans_4He_err[3] = {0.00756, 0.01666, 0.04229}; // error on mc_trans_4He
	Double_t data_trans_4He_err[3] = {0.01667, 0.03469, 0.09516}; // error on data_trans_4He
        Double_t mc_trans_C12[3] = {0.716, 0.653, 0.559}; // GENIE transparencies for C-12
        Double_t data_trans_C12[3] = {0.761, 0.711, 0.719}; // CLAS transparencies for C-12
        Double_t mc_trans_C12_err[3] = {0.00411, 0.01207, 0.02844}; // error on mc_trans_C12
        Double_t data_trans_C12_err[3] = {0.01069, 0.03378, 0.08786}; // error on data_trans_C12
        Double_t mc_trans_56Fe[3] = {0.506, 0.385, 0.316}; // GENIE transparencies for Fe-56
        Double_t data_trans_56Fe[3] = {0.512, 0.397, 0.330}; // CLAS transparencies for Fe-56
        Double_t mc_trans_56Fe_err[3] = {0.00501, 0.00931, 0.01898}; // error on mc_trans_56Fe
        Double_t data_trans_56Fe_err[3] = {0.01384, 0.02760, 0.06348}; // error on data_trans_56Fe
	Double_t range_num[3] = {1, 2, 3};
	Double_t range_num_err[3] = {0, 0, 0};

	Double_t mc_trans[3] = {0.78834, 0.70816, 0.48154}; // GENIE transparencies
	Double_t data_trans[3] = {0.81585, 0.75635, 0.49075}; // CLAS data transparencies
	Double_t mc_err[3] = {0.00679, 0.00386, 0.00436}; // error on mc_trans
	Double_t data_err[3] = {0.01485, 0.01013, 0.01225}; // error on data_trans
	Double_t nucleon_num[3] = {4, 12, 56}; // nucleon numbers for He-4, C-12, and Fe-56, respectively
	Double_t nucleon_err[3] = {0, 0, 0}; // there is no error on nucleon number; we know that helium is helium, etc.

	// put the data into TGraphErrors objects
	auto mc_graph_4He = new TGraphErrors(3, range_num, mc_trans_4He, range_num_err, mc_trans_4He_err);
	auto data_graph_4He = new TGraphErrors(3, range_num, data_trans_4He, range_num_err, data_trans_4He_err); 
        auto mc_graph_C12 = new TGraphErrors(3, range_num, mc_trans_C12, range_num_err, mc_trans_C12_err);
        auto data_graph_C12 = new TGraphErrors(3, range_num, data_trans_C12, range_num_err, data_trans_C12_err);
        auto mc_graph_56Fe = new TGraphErrors(3, range_num, mc_trans_56Fe, range_num_err, mc_trans_56Fe_err);
        auto data_graph_56Fe = new TGraphErrors(3, range_num, data_trans_56Fe, range_num_err, data_trans_56Fe_err);

	auto mc_graph = new TGraphErrors(3, nucleon_num, mc_trans, nucleon_err, mc_err);
        auto data_graph = new TGraphErrors(3, nucleon_num, data_trans, nucleon_err, data_err);

        TGraphErrors *graphs[8] = {mc_graph_4He, data_graph_4He, mc_graph_C12, data_graph_C12, mc_graph_56Fe, data_graph_56Fe, mc_graph, data_graph}; // put all of the targets' graphs into an array so they're easier to process in a for loop

        // create a canvas on which to draw the histograms
        TCanvas* c;
        c = new TCanvas( TString::Format("c"), TString::Format("c"), 800, 600);
	c->SetGridy();

	// format the histograms
	mc_graph_4He->SetTitle("Proton Transparency vs. #theta_{electron} Range (He-4, GENIE & CLAS)");
        mc_graph_C12->SetTitle("Proton Transparency vs. #theta_{electron} Range (C-12, GENIE & CLAS)");
        mc_graph_56Fe->SetTitle("Proton Transparency vs. #theta_{electron} Range (Fe-56, GENIE & CLAS)");
        data_graph_4He->SetTitle("Proton Transparency vs. #theta_{electron} Range (He-4, GENIE & CLAS)");
        data_graph_C12->SetTitle("Proton Transparency vs. #theta_{electron} Range (C-12, GENIE & CLAS)");
        data_graph_56Fe->SetTitle("Proton Transparency vs. #theta_{electron} Range (Fe-56, GENIE & CLAS)");
        mc_graph->SetTitle("Proton Transparency vs. Nucleon Number (GENIE & CLAS)");
        data_graph->SetTitle("Proton Transparency vs. Nucleon Number (GENIE & CLAS)");

	mc_graph_4He->GetXaxis()->SetTitle("#theta_{electron} Range");
        mc_graph_C12->GetXaxis()->SetTitle("#theta_{electron} Range");
        mc_graph_56Fe->GetXaxis()->SetTitle("#theta_{electron} Range");
        data_graph_4He->GetXaxis()->SetTitle("#theta_{electron} Range");
        data_graph_C12->GetXaxis()->SetTitle("#theta_{electron} Range");
        data_graph_56Fe->GetXaxis()->SetTitle("#theta_{electron} Range");
        mc_graph->GetXaxis()->SetTitle("Nucleon Number");
        data_graph->GetXaxis()->SetTitle("Nucleon Number");

	mc_graph_4He->GetXaxis()->SetRangeUser( 0.75, 3.25);
        mc_graph_C12->GetXaxis()->SetRangeUser( 0.75, 3.25);
        mc_graph_56Fe->GetXaxis()->SetRangeUser( 0.75, 3.25);
        data_graph_4He->GetXaxis()->SetRangeUser( 0.75, 3.25);
        data_graph_C12->GetXaxis()->SetRangeUser( 0.75, 3.25);
        data_graph_56Fe->GetXaxis()->SetRangeUser( 0.75, 3.25);
        mc_graph->GetXaxis()->SetRangeUser( 1, 100);
        data_graph->GetXaxis()->SetRangeUser( 1, 100);

	mc_graph_4He->GetXaxis()->SetNdivisions( 4);
        mc_graph_C12->GetXaxis()->SetNdivisions( 4);
        mc_graph_56Fe->GetXaxis()->SetNdivisions( 4);
        data_graph_4He->GetXaxis()->SetNdivisions( 4);
        data_graph_C12->GetXaxis()->SetNdivisions( 4);
        data_graph_56Fe->GetXaxis()->SetNdivisions( 4);

	mc_graph_4He->GetYaxis()->SetRangeUser( 0.65, 0.85);
	data_graph_4He->GetYaxis()->SetRangeUser( 0.65, 0.85);
        mc_graph_C12->GetYaxis()->SetRangeUser( 0.55, 0.78);
        data_graph_C12->GetYaxis()->SetRangeUser( 0.55, 0.78);
        mc_graph_56Fe->GetYaxis()->SetRangeUser( 0.29, 0.53);
        data_graph_56Fe->GetYaxis()->SetRangeUser( 0.29, 0.53);
	mc_graph->GetYaxis()->SetRangeUser( 0.4, 1.1);
        data_graph->GetYaxis()->SetRangeUser( 0.4, 1.1);

	mc_graph_4He->SetMarkerStyle( kOpenSquare);
        mc_graph_C12->SetMarkerStyle( kOpenSquare);
        mc_graph_56Fe->SetMarkerStyle( kOpenSquare);
        data_graph_4He->SetMarkerStyle( kOpenCircle);
        data_graph_C12->SetMarkerStyle( kOpenCircle);
        data_graph_56Fe->SetMarkerStyle( kOpenCircle);
	mc_graph->SetMarkerStyle( kFullSquare);
	data_graph->SetMarkerStyle( kFullCircle);

	// this stuff in the for loop is all the same for all of the graphs array
	for (int i = 0; i < 8; i++) {
		graphs[i]->SetLineColor( kBlack);
		graphs[i]->GetYaxis()->SetTitle("Proton Transparency (True QE / Inclusive)");
		graphs[i]->GetXaxis()->CenterTitle( true);
                graphs[i]->GetYaxis()->CenterTitle( true);

		// format GENIE and CLAS slightly different
		if (i%2 == 0) { graphs[i]->SetMarkerColor( kBlue); }
		else { graphs[i]->SetMarkerColor( kBlack); }
	}

	gStyle->SetErrorX( 0); // we manually set it so that there is no x error, but this is just a failsafe

        // draw the He-4 histograms
        mc_graph_4He->Draw("AP"); // draw the axes with this first one
        data_graph_4He->Draw("P"); // just plot the points

        // create and draw a legend for our non-target-specific histograms
        TLegend *legend_4He = new TLegend( 0.15, 0.15, 0.275, 0.275);
        legend_4He->AddEntry( mc_graph_4He, "GENIE MC");
        legend_4He->AddEntry( data_graph_4He, "CLAS Data");
        legend_4He->Draw();

        // save the histogram
        c->Update();
        c->SaveAs(TString::Format("../output/final_trans_plots/trans_vs_4HeRange_beamEnergy2261_exclPlus_FINAL_PLOT.pdf"));

        // draw the C-12 histograms
        mc_graph_C12->Draw("AP"); // draw the axes with this first one
        data_graph_C12->Draw("P"); // just plot the points

        // create and draw a legend for our non-target-specific histograms
        TLegend *legend_C12 = new TLegend( 0.15, 0.15, 0.275, 0.275);
        legend_C12->AddEntry( mc_graph_C12, "GENIE MC");
        legend_C12->AddEntry( data_graph_C12, "CLAS Data");
        legend_C12->Draw();

        // save the histogram
        c->Update();
        c->SaveAs(TString::Format("../output/final_trans_plots/trans_vs_C12Range_beamEnergy2261_exclPlus_FINAL_PLOT.pdf"));

        // draw the Fe-56 histograms
        mc_graph_56Fe->Draw("AP"); // draw the axes with this first one
        data_graph_56Fe->Draw("P"); // just plot the points

        // create and draw a legend for our non-target-specific histograms
        TLegend *legend_56Fe = new TLegend( 0.15, 0.15, 0.275, 0.275);
        legend_56Fe->AddEntry( mc_graph_56Fe, "GENIE MC");
        legend_56Fe->AddEntry( data_graph_56Fe, "CLAS Data");
        legend_56Fe->Draw();

        // save the histogram
	c->Update();
        c->SaveAs(TString::Format("../output/final_trans_plots/trans_vs_56FeRange_beamEnergy2261_exclPlus_FINAL_PLOT.pdf"));

        // draw the non-target-specific histograms
	mc_graph->Draw("AP"); // draw the axes with this first one
	data_graph->Draw("P"); // just plot the points

	// create and draw a legend for our non-target-specific histograms
    	TLegend *legend = new TLegend( 0.75, 0.75, 0.875, 0.875);
    	legend->AddEntry( mc_graph, "GENIE MC");
    	legend->AddEntry( data_graph, "CLAS Data");
    	legend->Draw();
	
	// save the histogram
        c->SetLogy();
        c->SetLogx();
        c->SetGridx();
	c->Update();
	c->SaveAs(TString::Format("../output/final_trans_plots/trans_vs_target_beamEnergy2261_exclPlus_FINAL_PLOT.pdf"));

        // print out some useful information about the sigma-level significance of disagreement between the
        // target ranges and average target transparencies
        Double_t sigma_diff_mcAndData_4He[3] = {-1, -1, -1};
        Double_t sigma_diff_mcAndData_C12[3] = {-1, -1, -1};
        Double_t sigma_diff_mcAndData_56Fe[3] = {-1, -1, -1};
        Double_t sigma_diff_data_4He[3] = {-1, -1, -1};
        Double_t sigma_diff_data_C12[3] = {-1, -1, -1};
        Double_t sigma_diff_data_56Fe[3] = {-1, -1, -1};
        for (int i = 0; i < 3; i++) {
                sigma_diff_mcAndData_4He[i] = std::abs(mc_trans_4He[i] - data_trans_4He[i]) / (mc_trans_4He_err[i] + data_trans_4He_err[i]);
                sigma_diff_data_4He[i] = std::abs(mc_trans_4He[i] - data_trans_4He[i]) / data_trans_4He_err[i];
                sigma_diff_mcAndData_C12[i] = std::abs(mc_trans_C12[i] - data_trans_C12[i]) / (mc_trans_C12_err[i] + data_trans_C12_err[i]);
                sigma_diff_data_C12[i] = std::abs(mc_trans_C12[i] - data_trans_C12[i]) / data_trans_C12_err[i];
                sigma_diff_mcAndData_56Fe[i] = std::abs(mc_trans_56Fe[i] - data_trans_56Fe[i]) / (mc_trans_56Fe_err[i] + data_trans_56Fe_err[i]);
                sigma_diff_data_56Fe[i] = std::abs(mc_trans_56Fe[i] - data_trans_56Fe[i]) / data_trans_56Fe_err[i];
        }

	std::cout << "PROGRAM OUTPUT ------------------------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "Sigma-level difference in He-4 (MC and data errors): " << sigma_diff_mcAndData_4He[0] << ", " << sigma_diff_mcAndData_4He[1] << ", " << sigma_diff_mcAndData_4He[2] << std::endl;
        std::cout << "Sigma-level difference in He-4 (just data errors): " << sigma_diff_data_4He[0] << ", " << sigma_diff_data_4He[1] << ", " << sigma_diff_data_4He[2] << std::endl;
        std::cout << "Sigma-level difference in C-12 (MC and data errors): " << sigma_diff_mcAndData_C12[0] << ", " << sigma_diff_mcAndData_C12[1] << ", " << sigma_diff_mcAndData_C12[2] << std::endl;
        std::cout << "Sigma-level difference in C-12 (just data errors): " << sigma_diff_data_C12[0] << ", " << sigma_diff_data_C12[1] << ", " << sigma_diff_data_C12[2] << std::endl;
        std::cout << "Sigma-level difference in Fe-56 (MC and data errors): " << sigma_diff_mcAndData_56Fe[0] << ", " << sigma_diff_mcAndData_56Fe[1] << ", " << sigma_diff_mcAndData_56Fe[2] << std::endl;
        std::cout << "Sigma-level difference in Fe-56 (just data errors): " << sigma_diff_data_56Fe[0] << ", " << sigma_diff_data_56Fe[1] << ", " << sigma_diff_data_56Fe[2] << std::endl;
}
