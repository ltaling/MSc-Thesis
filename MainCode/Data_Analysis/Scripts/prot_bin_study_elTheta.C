// The purpose of this script is to extract the h2_prot_mom_el_theta
// plot from a given run of genie_analysis and investigate how best
// to refine the inclusive cuts made on the electron. This script
// produces the TH2 described above as well as a partially cropped 
// version of that TH2 (using the proton momentum cut made in the 
// exclusive+ analysis). Additionally, we project this un-cut and
// partially cut TH2 onto the non-proton-momentum axis.
//
// Author: Jacob Smith (smithja)
// Creation Date: 11/15/2021

// These include statements were in the files that I received Graham.
// Delete/Modify them at your own peril.
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
#include <stdio.h> // needed to sort array of bin values
#include <stdlib.h> // needed to sort array of bin values

#include "parse_file.C"

// -----------------------------------------------------------------------------------------------------------------------
// This function compares two elements to use in conjunction with
// qsort(), which will sort the array of bin values that the user
// specifies during run time.
int elTheta_comp (const void * elem1, const void * elem2)
{
    Double_t f = *((Double_t*)elem1);
    Double_t s = *((Double_t*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

// One of the two main functions for this script. This function should be used when calling this script individually from the command line.
void prot_bin_study_elTheta( std::string file_name){
 
   // Give a little preamble so the user can choose the bins they want to use.
    std::cout << std::endl << std::endl; // separate the preamble message
    Int_t n_bins;

    std::cout << "How many bins would you like to make? ";
    std::cin >> n_bins;
    const Int_t N_BINS = n_bins; // set as a constant so program doesn't output warnings about using a non-constant for array declaration
    
    Double_t low_bin_bounds[N_BINS];
    Double_t high_bin_bounds[N_BINS];
    
    std::cout << "\nPlease specify " << N_BINS << " edges for your bins." << std::endl;
    std::cout << "For example, entering 1, 2, and 3, 4 will give you two bins, one from 1 to 2 and one from 3 to 4." << std::endl;
    for (int i = 0; i < N_BINS; i++) {
        std::cout << "Low Bin edge " << i << ": ";
        std::cin >> low_bin_bounds[i];
        std::cout << "High Bin edge " << i << ": ";
        std::cin >> high_bin_bounds[i];
    }
 
    std::cout << "\n\nSorting your bins..." << std::endl;
    qsort( low_bin_bounds, sizeof(low_bin_bounds)/sizeof(*low_bin_bounds), sizeof(*low_bin_bounds), elTheta_comp);
    qsort( high_bin_bounds, sizeof(high_bin_bounds)/sizeof(*high_bin_bounds), sizeof(*high_bin_bounds), elTheta_comp);
   
    std::cout << "You entered " << N_BINS << " bins:" << std::endl;
    for (int i = 0; i < N_BINS; i++) {
        std::cout << "{" << low_bin_bounds[i] << ", " << high_bin_bounds[i] << "}" << "\n";
    }

    std::cout << std::endl << std::endl;
    std::cout << "END OF USER INPUT SECTION" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl << std::endl;
    // END OF USER INPUT ---------------------------------------------------------------------------------------------------------------------------------------------------

    // CAUTION: These rebin_factor variables need to be large enough so that this
    //          program is not overwhelmed by the amount of statistics originally
    //          in the TH2 extracted from the ROOT file passed to this script.
    //          As a courtesey to others, I have documented value pairs for which 
    //          this program executes.
    // Form (rebin_factor_x, rebin_factor_y):
    Int_t rebin_factor_x = 1;
    Int_t rebin_factor_y = 50;

    std::string info ("");
    std::string cuts ("");

    // import the ROOT data file from the output of doing ./genie_analysis in the e4nu directory
    TFile *input1 = TFile::Open( TString( file_name));

    // Grab the TList in the ROOT File for info on the run
    TList *Run_Info = (TList*)input1->Get("Run_Info");

    std::vector<std::string> parse = parse_file(Run_Info);

    info = parse[0];
    cuts = parse[1];

    // Electron theta cuts
    double theta_el_lb = ((TParameter<double>*)Run_Info->FindObject("El theta lb"))->GetVal();
    double theta_el_ub = ((TParameter<double>*)Run_Info->FindObject("El theta ub"))->GetVal();
    double protMom_lb = ((TParameter<double>*)Run_Info->FindObject("Prot mom lb"))->GetVal();
    std::string run_type = std::string(((TNamed*)Run_Info->FindObject("Run type"))->GetTitle());
    if(theta_el_lb < 20.0) theta_el_lb = 20.0;
    if(theta_el_ub > 60.0) theta_el_ub = 60.0;
    
    // histogram initialization
    TH2F* h2_prot_mom_el_theta = (TH2F*)input1->Get(TString::Format("h2_prot_mom_el_theta"));
    
    // Rebin to reduce amount of data
    h2_prot_mom_el_theta->RebinX( rebin_factor_x); 
    h2_prot_mom_el_theta->RebinY( rebin_factor_y); 

    TH2F* bin_selected_h2_prot_mom_el_theta[N_BINS];
    //Changing this to look at bins of proton momentum
    TH1D* bin_selected_h1_prot_mom_proj[N_BINS];
    Int_t n_xbins, bin_for_lb, bin_for_ub;

    double max_1D = 0.0;
    
    // Actually we're going to project out a certain range of protom momenta corresponding to a certain range of electron theta (from bin_for_lb to bin_for_ub)
    for (int i = 0; i < N_BINS; i++) {
        bin_selected_h2_prot_mom_el_theta[i] = (TH2F*)h2_prot_mom_el_theta->Clone(); // histogram will have proton momentum bin cuts applied
        
    	n_xbins = bin_selected_h2_prot_mom_el_theta[i]->GetNbinsX();
    	bin_for_lb = (Int_t)( (low_bin_bounds[i] / 180) * n_xbins ); // factor of 6 because proton momentum ranges from 0 to 6 
    	bin_for_ub = (Int_t)( (high_bin_bounds[i] / 180) * n_xbins ); // factor of 6 because proton momentum ranges from 0 to 6

        bin_selected_h1_prot_mom_proj[i] = (TH1D*)h2_prot_mom_el_theta->ProjectionY( TString::Format("h2_prot_mom_el_theta_bin%i_py", i+1), bin_for_lb, bin_for_ub)->Clone(); // project TH2 onto x-axis, only take values in given proton momentum band
        
        if(bin_selected_h1_prot_mom_proj[i]->GetMaximum() > max_1D) max_1D = bin_selected_h1_prot_mom_proj[i]->GetMaximum();

    }
    
    // format the histograms
    h2_prot_mom_el_theta->SetStats( 0);
    h2_prot_mom_el_theta->GetXaxis()->SetRangeUser(theta_el_lb, theta_el_ub);
    //h2_prot_mom_el_theta->GetYaxis()->SetRangeUser(protMom_lb, bin_bounds[N_BINS]);
    h2_prot_mom_el_theta->GetZaxis()->SetRangeUser( 1, h2_prot_mom_el_theta->GetMaximum()); 
    h2_prot_mom_el_theta->GetXaxis()->SetTitle("#theta_{electron} [deg]");
    h2_prot_mom_el_theta->GetYaxis()->SetTitle("Proton Momentum [GeV]");
    h2_prot_mom_el_theta->GetXaxis()->CenterTitle( true);
    h2_prot_mom_el_theta->GetYaxis()->CenterTitle( true);
    h2_prot_mom_el_theta->SetTitle(TString(info)+TString(cuts));


    for (int i = 0; i < N_BINS; i++) {
    	if (N_BINS != 1) { 
            bin_selected_h2_prot_mom_el_theta[i]->SetStats(0);
        	//bin_selected_h2_prot_mom_el_theta[i]->GetXaxis()->SetRangeUser( theta_el_lb, theta_el_ub);
            bin_selected_h2_prot_mom_el_theta[i]->GetXaxis()->SetRangeUser( low_bin_bounds[i], high_bin_bounds[i]);
        	//bin_selected_h2_prot_mom_el_theta[i]->GetYaxis()->SetRangeUser( protMom_lb, high_bin_bounds[i]);
    	    bin_selected_h2_prot_mom_el_theta[i]->GetZaxis()->SetRangeUser( 1,bin_selected_h2_prot_mom_el_theta[i]->GetMaximum()); 
        	bin_selected_h2_prot_mom_el_theta[i]->GetXaxis()->SetTitle("#theta_{electron} [deg]");
        	bin_selected_h2_prot_mom_el_theta[i]->GetYaxis()->SetTitle("Proton Momentum [GeV]");
        	bin_selected_h2_prot_mom_el_theta[i]->GetXaxis()->CenterTitle( true);
        	bin_selected_h2_prot_mom_el_theta[i]->GetYaxis()->CenterTitle( true);
        	bin_selected_h2_prot_mom_el_theta[i]->SetTitle("Electron theta isolated "+TString(info)+TString(cuts)+" ("+TString( std::to_string( low_bin_bounds[i]))+" deg < #theta_{e} < "+TString( std::to_string( high_bin_bounds[i]))+" deg)");
    	}

        bin_selected_h1_prot_mom_proj[i]->SetStats( 0);
    	bin_selected_h1_prot_mom_proj[i]->Sumw2();
        //bin_selected_h1_prot_mom_proj[i]->GetXaxis()->SetRangeUser( protMom_lb, theta_el_ub);
        bin_selected_h1_prot_mom_proj[i]->GetYaxis()->SetRangeUser(.01,max_1D);
        bin_selected_h1_prot_mom_proj[i]->GetXaxis()->SetTitle("p_{p} [GeV]");
       	bin_selected_h1_prot_mom_proj[i]->GetYaxis()->SetTitle("Events");
        bin_selected_h1_prot_mom_proj[i]->GetXaxis()->CenterTitle( true);
        bin_selected_h1_prot_mom_proj[i]->GetYaxis()->CenterTitle( true);
    
        if (run_type.find("Data") !=  string::npos)  {
            bin_selected_h1_prot_mom_proj[i]->SetMarkerColor( kBlack);
            bin_selected_h1_prot_mom_proj[i]->SetMarkerStyle( kFullCircle);
            bin_selected_h1_prot_mom_proj[i]->SetLineColor( kBlack);
        }
	   bin_selected_h1_prot_mom_proj[i]->SetTitle(TString(info)+TString(cuts)+" ("+TString( std::to_string( low_bin_bounds[i]))+" deg < #theta_{e} < "+TString( std::to_string( high_bin_bounds[i]))+" deg)");
    }

    gStyle->SetErrorX(0);

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // draw the histogram that does not separate the data into proton bands and save it as a PDF file
    h2_prot_mom_el_theta->Draw("colz");
    c->SetLeftMargin( 0.14);
    c->SetTopMargin( 0.16);
    c->SetRightMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->SetLogz();
    c->Update();    

    std::string save_file_path ("invalid");
    if (run_type.find("Genie") !=  string::npos){ 
	   save_file_path = "../output/prot_bin_study_elTheta/"+TString(info)+TString(cuts)+"protmom_vs_eltheta_MottXSecEq1.pdf"; 
    }
    else if (run_type.find("Data") !=  string::npos) { 
	   save_file_path = "../output/prot_bin_study_elTheta/"+TString(info)+TString(cuts)+".pdf"; 
    }
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));
    
    // reset the top margin because the bin_selected plots won't have as big titles
    c->SetTopMargin( 0.1);
    c->Update();
    
    TLegend* lg2 = new TLegend( 0.6, 0.6, 0.9, 0.9 );
    // this is the loop to name all of the proton band separated files for GENIE MC
    if (run_type.find("Genie") !=  string::npos) {
        for (int i = 0; i < N_BINS; i++) {
	       //if (N_BINS != 1) { // if N_BINS == 1, then this plot will be the same as the TH2 saved above; we don't want duplicate plots
           //bin_selected_h2_prot_mom_el_theta[i]->Draw("colz");
	       //save_file_path = "../output/prot_bin_study_elTheta/"+TString(info)+TString(cuts)+"_bin"+std::to_string( bin_bounds[i])+"to"+std::to_string( bin_bounds[i+1])+"_MottXSecEq1.pdf";
	       //c->SetLogy( 0); // set y scale to be linear
	       //c->SaveAs( TString::Format("%s", save_file_path.c_str()));
	       //}

           lg2->AddEntry(bin_selected_h1_prot_mom_proj[i], (std::to_string( low_bin_bounds[i])+" to "+std::to_string( high_bin_bounds[i])).c_str(), "l");

	       if (i == 0) bin_selected_h1_prot_mom_proj[i]->Draw("PLC HIST");
           else bin_selected_h1_prot_mom_proj[i]->Draw("PLC HIST SAME");
       }
        save_file_path = "../output/prot_bin_study_elTheta/"+TString(info)+TString(cuts)+"prot_mom_bin_studies_protmom_vs_eltheta_MottXSecEq1.pdf";
        lg2->Draw( "same" );
        c->SetLogy();
        c->SaveAs( TString::Format("%s", save_file_path.c_str()));
    }
    // this is the loop to name all of the proton band separated files for CLAS data
    else if (run_type.find("Data") !=  string::npos) {
        for (int i = 0; i < N_BINS; i++) {
            if (N_BINS != 1) { // if N_BINS == 1, then this plot will be the same as the TH2 saved above; we don't want duplicate plots
		bin_selected_h2_prot_mom_el_theta[i]->Draw("colz");
            	save_file_path = "../output/prot_bin_study_elTheta/"+TString(info)+TString(cuts)+"_bin"+std::to_string( low_bin_bounds[i])+"to"+std::to_string( high_bin_bounds[i])+".pdf";
     	    	c->SetLogy( 0); // set y scale to be linear
            	c->SaveAs( TString::Format("%s", save_file_path.c_str()));
	    }

            bin_selected_h1_prot_mom_proj[i]->Draw("E1");
            save_file_path = "../output/prot_bin_study_elTheta/"+TString(info)+TString(cuts)+"_bin"+std::to_string( low_bin_bounds[i])+"to"+std::to_string( high_bin_bounds[i])+".pdf";
	    c->SetLogy();
	    c->SaveAs( TString::Format("%s", save_file_path.c_str()));
        }
    }
}
/*

// One of the two main functions for this script. This function should be used when calling this script from the master script, prot_bin_study.C.
void prot_bin_study_elTheta( std::string file_name, Double_t bin_bounds[], const Int_t N_BINS){
    // CAUTION: These rebin_factor variables need to be large enough so that this
    //          in the TH2 extracted from the ROOT file passed to this script.
    //          As a courtesey to others, I have documented value pairs for which 
    //          this program executes.
    // Form (rebin_factor_x, rebin_factor_y):
    Int_t rebin_factor_x = 1;
    Int_t rebin_factor_y = 10;
    
    // read the file name to determine if we are plotting GENIE MC or CLAS data
    auto file_type_start_pos = file_name.rfind("/"); // find the last instance of '/'
    auto file_type_end_pos = file_name.find("_", file_type_start_pos + 1);
    std::string file_type = file_name.substr( file_type_start_pos + 1, file_type_end_pos - (file_type_start_pos + 1));
    std::string file_type_plot ("invalid");

    if (file_type.compare("genie") == 0) { file_type_plot = "GENIE MC"; }
    else if (file_type.compare("data") == 0) { file_type_plot = "data"; }
    else { // if we do not have GENIE MC or CLAS data
            std::cout << "/----------------------------------------------------------------------------------------\\" << std::endl;
            std::cout << "| It seems you are trying to process a file with an invalid file name. Please make sure  |" << std::endl;
            std::cout << "| you have either the keyword 'genie' or the keyword 'data' at the beginning of the file |" << std::endl;
            std::cout << "| name, e.g. <PATH_TO_FILE>/genie_<OMITTED>.root or <PATH_TO_FILE>/data_<OMITTED>.root.  |" << std::endl;
            std::cout << "\\----------------------------------------------------------------------------------------/" << std::endl;
            return;
    }

    // declare variables that will be used in the plotting/formatting of histograms and file names
    std::string target ("invalid");
    std::string target_plot ("invalid");
    std::string beam_energy ("invalid");
    std::string beam_energy_plot ("invalid");
    std::string range ("invalid");
    std::string range_plot ("invalid");

    std::string elSectors ("invalid");
    std::string elSectors_plot ("invalid");
    std::string deltaPhiEl ("invalid");
    std::string deltaPhiEl_plot ("invalid");
    std::string thetaEl_lb ("invalid");
    std::string thetaEl_ub ("invalid");
    std::string thetaEl_lb_plot ("invalid");
    std::string thetaEl_ub_plot ("invalid");
    std::string elMom_lb ("invalid");
    std::string elMom_lb_plot ("invalid");
    std::string protSectors ("invalid");
    std::string protSectors_plot ("invalid");
    std::string deltaPhiProt ("invalid");
    std::string deltaPhiProt_plot ("invalid");
    std::string thetaProt_lb ("invalid");
    std::string thetaProt_ub ("invalid");
    std::string thetaProt_lb_plot ("invalid");
    std::string thetaProt_ub_plot ("invalid");
    std::string protMom_lb ("invalid");
    std::string protMom_lb_plot ("invalid");
    std::string protMom_ub ("invalid");
    std::string protMom_ub_plot ("invalid");

    // import the ROOT data file from the output of doing ./genie_analysis in the e4nu directory
    TFile *input1 = TFile::Open( TString( file_name));
    
    // Parse the file_name to extract the necessary information that routinely changes for plotting, formatting, etc.
    // Note that the hard-coded numbers for each grouping represent the length of the keyword in question.
    auto target_kwd_pos = file_name.find("target");
    if (target_kwd_pos == std::string::npos) { std::cout << "'target' keyword not found in file_name. Correct file_name!" << std::endl; }
    auto target_end_pos = file_name.find("_", target_kwd_pos);
    target = file_name.substr( target_kwd_pos + 6, target_end_pos - (target_kwd_pos + 6)); // shift by 6 because length of "target" is 6
    if (target.compare("4He") == 0) { target_plot = "He-4"; }
    else if (target.compare("C12") == 0) { target_plot = "C-12"; }
    else if (target.compare("56Fe") == 0) { target_plot = "Fe-56"; }
    else { std::cout << "You have entered an invalid/unknown target. Consider adding code to accomodate this target in prot_bin_study_elTheta.C." << std::endl; }
    //std::cout << "target: " << target << "; target_plot: " << target_plot << std::endl;

    auto beam_kwd_pos = file_name.find("beamEnergy");
    if (beam_kwd_pos == std::string::npos) { std::cout << "'beamEnergy' keyword not found in file_name. Correct file_name!" << std::endl; }
    auto beam_end_pos = file_name.find("_", beam_kwd_pos);
    beam_energy = file_name.substr( beam_kwd_pos + 10, beam_end_pos - (beam_kwd_pos + 10));
    if (beam_energy.compare("1161") == 0) { beam_energy_plot = "1.161 GeV"; }
    else if (beam_energy.compare("2261") == 0) { beam_energy_plot = "2.261 GeV"; }
    else if (beam_energy.compare("4461") == 0) { beam_energy_plot = "4.461 GeV"; }
    else { std::cout << "You have entered an invalid/unknown beam energy. Consider adding code to accomodate this beam energy in prot_bin_study_elTheta.C." << std::endl; }
    //std::cout << "beam_energy: " << beam_energy << "; beam_energy_plot: " << beam_energy_plot << std::endl;
    
    auto elSectors_kwd_pos = file_name.find("elSectors");
    if (elSectors_kwd_pos != std::string::npos){
        auto elSectors_end_pos = file_name.find("_", elSectors_kwd_pos);
        elSectors = file_name.substr( elSectors_kwd_pos + 9, elSectors_end_pos - (elSectors_kwd_pos + 9));
        unsigned int i;
        elSectors_plot = ""; // elSectors is valid so first do away with 'invalid'
        for (i = 0; i < elSectors.length() - 1; i++){ elSectors_plot = elSectors_plot + elSectors.substr( i, 1) + " "; }
        elSectors_plot = elSectors_plot + "& " + elSectors.substr( i, 1);
    }
    //std::cout << "elSectors: " << elSectors << "; elSectors_plot: " << elSectors_plot << std::endl;
    
    auto deltaPhiElEq_kwd_pos = file_name.find("deltaPhiElEq");
    if (deltaPhiElEq_kwd_pos != std::string::npos){
        auto deltaPhiElEq_end_pos = file_name.find("_", deltaPhiElEq_kwd_pos);
        deltaPhiEl = file_name.substr( deltaPhiElEq_kwd_pos + 12, deltaPhiElEq_end_pos - (deltaPhiElEq_kwd_pos + 12));
	deltaPhiEl_plot = elTheta_pt_plot_converter( deltaPhiEl);
    }
    //std::cout << "deltaPhiEl: " << deltaPhiEl << ", deltaPhiEl_plot: " << deltaPhiEl_plot << std::endl;
    
    auto thetaEl_kwd_pos = file_name.find("thetaEl");
    if (thetaEl_kwd_pos != std::string::npos){
        auto thetaEl_lb_pos = thetaEl_kwd_pos + 7;
        auto thetaEl_div_pos = file_name.find("to", thetaEl_kwd_pos);
        auto thetaEl_ub_pos = thetaEl_div_pos + 2;
        auto thetaEl_end_pos = file_name.find("_", thetaEl_kwd_pos);
        thetaEl_lb = file_name.substr( thetaEl_lb_pos, thetaEl_div_pos - thetaEl_lb_pos);
        thetaEl_ub = file_name.substr( thetaEl_ub_pos, thetaEl_end_pos - thetaEl_ub_pos);
	thetaEl_lb_plot = elTheta_pt_plot_converter( thetaEl_lb);
        thetaEl_ub_plot = elTheta_pt_plot_converter( thetaEl_ub);
        if (elTheta_pt_converter( thetaEl_lb) == 20 && elTheta_pt_converter( thetaEl_ub) == 25) { 
	    range = "range1"; 
	    range_plot = "Range 1";
	}
        if (elTheta_pt_converter( thetaEl_lb) == 30 && elTheta_pt_converter( thetaEl_ub) == 35) { 
            range = "range2";  
            range_plot = "Range 2";
        }
        if (elTheta_pt_converter( thetaEl_lb) == 40 && elTheta_pt_converter( thetaEl_ub) == 45) { 
            range = "range3";  
            range_plot = "Range 3";
        }
    }
    //std::cout << "thetaEl_lb: " << thetaEl_lb << ", thetaEl_lb_plot: " << thetaEl_lb_plot << "; thetaEl_ub: " << thetaEl_ub << ", thetaEl_ub_plot: " << thetaEl_ub_plot << std::endl;
    
    auto elMom_kwd_pos = file_name.find("elMomLBEq");
    if (elMom_kwd_pos != std::string::npos){
        auto elMom_whl_pos = elMom_kwd_pos + 9;
        auto elMom_end_pos = file_name.find("_", elMom_kwd_pos);
        elMom_lb = file_name.substr( elMom_kwd_pos + 9, elMom_end_pos - (elMom_kwd_pos + 9));
        elMom_lb_plot = elTheta_pt_plot_converter( elMom_lb) + " GeV";
    }
    //std::cout << "elMom_lb: " << elMom_lb << "; elMom_lb_plot: " << elMom_lb_plot << std::endl;

    auto protSectors_kwd_pos = file_name.find("protSectors");
    if (protSectors_kwd_pos != std::string::npos){
        auto protSectors_end_pos = file_name.find("_", protSectors_kwd_pos);
        protSectors = file_name.substr( protSectors_kwd_pos + 11, protSectors_end_pos - (protSectors_kwd_pos + 11));
        unsigned int i;
        protSectors_plot = ""; // protSectors is valid so first do away with 'invalid'
        for (i = 0; i < protSectors.length() - 1; i++){ protSectors_plot = protSectors_plot + protSectors.substr( i, 1) + " "; }
        protSectors_plot = protSectors_plot + "& " + protSectors.substr( i, 1);
    }
    //std::cout << "protSectors: " << protSectors << "; protSectors_plot: " << protSectors_plot << std::endl;

    auto deltaPhiProtEq_kwd_pos = file_name.find("deltaPhiProtEq");
    if (deltaPhiProtEq_kwd_pos != std::string::npos){
        auto deltaPhiProtEq_end_pos = file_name.find("_", deltaPhiProtEq_kwd_pos);
        deltaPhiProt = file_name.substr( deltaPhiProtEq_kwd_pos + 14, deltaPhiProtEq_end_pos - (deltaPhiProtEq_kwd_pos + 14));
        deltaPhiProt_plot = elTheta_pt_plot_converter( deltaPhiProt);
    }
    //std::cout << "deltaPhiProt: " << deltaPhiProt << ", deltaPhiProt_plot: " << deltaPhiProt_plot << std::endl;
    
    auto thetaProt_kwd_pos = file_name.find("thetaProt");
    if (thetaProt_kwd_pos != std::string::npos){
        auto thetaProt_lb_pos = thetaProt_kwd_pos + 9;
        auto thetaProt_div_pos = file_name.find("to", thetaProt_kwd_pos);
        auto thetaProt_ub_pos = thetaProt_div_pos + 2;
        auto thetaProt_end_pos = file_name.find("_", thetaProt_kwd_pos);
        thetaProt_lb = file_name.substr( thetaProt_lb_pos, thetaProt_div_pos - thetaProt_lb_pos);
        thetaProt_ub = file_name.substr( thetaProt_ub_pos, thetaProt_end_pos - thetaProt_ub_pos);
        thetaProt_lb_plot = elTheta_pt_plot_converter( thetaProt_lb);
        thetaProt_ub_plot = elTheta_pt_plot_converter( thetaProt_ub);
    }
    //std::cout << "thetaProt_lb: " << thetaProt_lb << ", thetaProt_lb_plot: " << thetaProt_lb_plot << "; thetaProt_ub: " << thetaProt_ub << ", thetaProt_ub_plot: " << thetaProt_ub_plot << std::endl;

    auto protMom_lb_kwd_pos = file_name.find("protMomLBEq");
    auto protMom_ub_kwd_pos = file_name.find("protMomUBEq");
    auto protMom_kwd_pos = file_name.find("protMom");
    if (protMom_lb_kwd_pos != std::string::npos){
        auto protMom_whl_pos = protMom_lb_kwd_pos + 11;
        auto protMom_div_pos = file_name.find("pt", protMom_lb_kwd_pos);
        auto protMom_dec_pos = protMom_div_pos + 2;
        auto protMom_end_pos = file_name.find("_", protMom_lb_kwd_pos);
        protMom_lb = file_name.substr( protMom_lb_kwd_pos + 11, protMom_end_pos - (protMom_lb_kwd_pos + 11));
        protMom_lb_plot = file_name.substr( protMom_whl_pos, protMom_div_pos - protMom_whl_pos) + "." + file_name.substr( protMom_dec_pos, protMom_end_pos - protMom_dec_pos) + " GeV";
    }
     //std::cout << "protMom_lb: " << protMom_lb << "; protMom_lb_plot: " << protMom_lb_plot << std::endl;
    else if (protMom_ub_kwd_pos != std::string::npos){
        auto protMom_whl_pos = protMom_ub_kwd_pos + 11;
        auto protMom_div_pos = file_name.find("pt", protMom_ub_kwd_pos);
        auto protMom_dec_pos = protMom_div_pos + 2;
        auto protMom_end_pos = file_name.find("_", protMom_ub_kwd_pos);
        protMom_ub = file_name.substr( protMom_ub_kwd_pos + 11, protMom_end_pos - (protMom_ub_kwd_pos + 11));
        protMom_ub_plot = file_name.substr( protMom_whl_pos, protMom_div_pos - protMom_whl_pos) + "." + file_name.substr( protMom_dec_pos, protMom_end_pos - protMom_dec_pos) + " GeV";
    }
    //std::cout << "protMom_ub: " << protMom_ub << "; protMom_ub_plot: " << protMom_ub_plot << std::endl;
    else if (protMom_kwd_pos != std::string::npos){
        auto protMom_lb_whl_pos = protMom_kwd_pos + 7;
        auto protMom_lb_div_pos = file_name.find("pt", protMom_lb_whl_pos);
        auto protMom_lb_dec_pos = protMom_lb_div_pos + 2;
        auto protMom_div_pos = file_name.find("to", protMom_lb_dec_pos);
        auto protMom_ub_whl_pos = protMom_div_pos + 2;
        auto protMom_ub_div_pos = file_name.find("pt", protMom_ub_whl_pos);
        auto protMom_ub_dec_pos = protMom_ub_div_pos + 2;
        auto protMom_end_pos = file_name.find("_", protMom_ub_dec_pos);
        protMom_lb = file_name.substr( protMom_lb_whl_pos, protMom_div_pos - protMom_lb_whl_pos);
        protMom_ub = file_name.substr( protMom_ub_whl_pos, protMom_end_pos - protMom_ub_whl_pos);
        protMom_lb_plot = elTheta_pt_plot_converter( protMom_lb);
        protMom_ub_plot = elTheta_pt_plot_converter( protMom_ub);
    }
    //std::cout << "protMom_lb: " << protMom_lb << ", protMom_lb_plot: " << protMom_lb_plot << "; protMom_ub: " << protMom_ub << ", protMom_ub_plot: " << protMom_ub_plot << std::endl;
    else {
        std::cout << "ERROR: You need an exclusive+ cut simulation to run this file, which is" << std::endl;
        std::cout << "       determined by having a proton momentum cut in the file name with" << std::endl;
        std::cout << "       the keyword 'protMomLBEq'." << std::endl;
        std::cout << "PROGRAM ENDED; NO PLOTS GENERATED" << std::endl;
        return;
    }
    
    // histogram initialization
    TH2F* h2_prot_mom_el_theta = (TH2F*)input1->Get(TString::Format("h2_prot_mom_el_theta"));
    
    // CAUTION: These rebinning statements are required for this program to execute. There are originally 6000 * 6000 bins for this TH2, which is cloned N_BINS many times below--that's a LOT of data.
    h2_prot_mom_el_theta->RebinX( rebin_factor_x); // speeds up code by getting rid of statistics we won't necessarily look at
    h2_prot_mom_el_theta->RebinY( rebin_factor_y); // speeds up code by getting rid of statistics we won't necessarily look at

    TH2F* bin_selected_h2_prot_mom_el_theta[N_BINS];
    TH1D* bin_selected_h1_el_theta_proj[N_BINS];
    Int_t n_ybins;
    Int_t bin_for_lb;
    Int_t bin_for_ub;
    for (int i = 0; i < N_BINS; i++) {
        bin_selected_h2_prot_mom_el_theta[i] = (TH2F*)h2_prot_mom_el_theta->Clone(); // histogram will have proton momentum bin cuts applied
        
	n_ybins = bin_selected_h2_prot_mom_el_theta[i]->GetNbinsY();
	bin_for_lb = (Int_t)( (bin_bounds[i] / 6) * n_ybins ); // factor of 6 because proton momentum ranges from 0 to 6 
	bin_for_ub = (Int_t)( (bin_bounds[i+1] / 6) * n_ybins ); // factor of 6 because proton momentum ranges from 0 to 6

        bin_selected_h1_el_theta_proj[i] = (TH1D*)h2_prot_mom_el_theta->ProjectionX( TString::Format("h2_prot_mom_el_theta_bin%i_px", i+1), bin_for_lb, bin_for_ub)->Clone(); // project TH2 onto x-axis, only take values in given proton momentum band
    }
    
    // format the histograms
    h2_prot_mom_el_theta->SetStats( 0);
    h2_prot_mom_el_theta->GetXaxis()->SetRangeUser( elTheta_pt_converter( thetaEl_lb), elTheta_pt_converter( thetaEl_ub));
    h2_prot_mom_el_theta->GetYaxis()->SetRangeUser( elTheta_pt_converter( protMom_lb), bin_bounds[N_BINS]);
    h2_prot_mom_el_theta->GetZaxis()->SetRangeUser( 1, h2_prot_mom_el_theta->GetMaximum()); // do not assign a color to bins with <= 0 events in them
    h2_prot_mom_el_theta->GetXaxis()->SetTitle("#theta_{electron} [deg]");
    h2_prot_mom_el_theta->GetYaxis()->SetTitle("Proton Momentum [GeV]");
    h2_prot_mom_el_theta->GetXaxis()->CenterTitle( true);
    h2_prot_mom_el_theta->GetYaxis()->CenterTitle( true);
    if (protMom_lb_kwd_pos != std::string::npos) {
        h2_prot_mom_el_theta->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, p_{p^{+}} > "+TString(protMom_lb_plot)+"}");
    }
    else if (protMom_ub_kwd_pos != std::string::npos) {
        h2_prot_mom_el_theta->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, p_{p^{+}} < "+TString(protMom_ub_plot)+"}");
    }
    else if (protMom_kwd_pos != std::string::npos) {
        h2_prot_mom_el_theta->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, "+TString(protMom_lb_plot)+" < p_{p^{+}} < "+TString(protMom_ub_plot)+"}");
    }

    for (int i = 0; i < N_BINS; i++) {
        if (N_BINS != 1) { // if N_BINS == 1, this plot will be the same as TH2 histogram above bit with a different file name; we want to avoid that
            bin_selected_h2_prot_mom_el_theta[i]->SetStats( 0);
            bin_selected_h2_prot_mom_el_theta[i]->GetXaxis()->SetRangeUser( elTheta_pt_converter( thetaEl_lb), elTheta_pt_converter( thetaEl_ub));
            bin_selected_h2_prot_mom_el_theta[i]->GetYaxis()->SetRangeUser( bin_bounds[i], bin_bounds[i+1]);
            bin_selected_h2_prot_mom_el_theta[i]->GetZaxis()->SetRangeUser( 1, bin_selected_h2_prot_mom_el_theta[i]->GetMaximum()); // do not assign a color to bins with <=0 events in them
            bin_selected_h2_prot_mom_el_theta[i]->GetXaxis()->SetTitle("#theta_{electron} [deg]");
            bin_selected_h2_prot_mom_el_theta[i]->GetYaxis()->SetTitle("Proton Momentum [GeV]");
            bin_selected_h2_prot_mom_el_theta[i]->GetXaxis()->CenterTitle( true);
            bin_selected_h2_prot_mom_el_theta[i]->GetYaxis()->CenterTitle( true);
            bin_selected_h2_prot_mom_el_theta[i]->SetTitle("Proton Momentum Isolated "+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam, "+TString(file_type_plot)+" "+TString(range_plot)+" ("+TString( elTheta_pt_plot_converter( std::to_string( bin_bounds[i])))+" GeV < p_{p^{+}} < "+TString( elTheta_pt_plot_converter( std::to_string( bin_bounds[i+1])))+" GeV)");
        }

        bin_selected_h1_el_theta_proj[i]->SetStats( 0);
        bin_selected_h1_el_theta_proj[i]->Sumw2();
        bin_selected_h1_el_theta_proj[i]->GetXaxis()->SetRangeUser( elTheta_pt_converter( thetaEl_lb), elTheta_pt_converter( thetaEl_ub));
        bin_selected_h1_el_theta_proj[i]->GetXaxis()->SetTitle("#theta_{electron} [deg]");
        bin_selected_h1_el_theta_proj[i]->GetYaxis()->SetTitle("Events");
        bin_selected_h1_el_theta_proj[i]->GetXaxis()->CenterTitle( true);
        bin_selected_h1_el_theta_proj[i]->GetYaxis()->CenterTitle( true);
        if (file_type.compare("data") == 0) {
            bin_selected_h1_el_theta_proj[i]->SetMarkerColor( kBlack);
            bin_selected_h1_el_theta_proj[i]->SetMarkerStyle( kFullCircle);
            bin_selected_h1_el_theta_proj[i]->SetLineColor( kBlack);
        }
	bin_selected_h1_el_theta_proj[i]->SetTitle(TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam, "+TString(file_type_plot)+" "+TString(range_plot)+" ("+TString( elTheta_pt_plot_converter( std::to_string( bin_bounds[i])))+" GeV < p_{p^{+}} < "+TString( elTheta_pt_plot_converter( std::to_string( bin_bounds[i+1])))+" GeV)");
    }
    gStyle->SetErrorX( 0);

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // draw the histogram that does not separate the data into proton bands and save it as a PDF file
    h2_prot_mom_el_theta->Draw("colz");
    c->SetLeftMargin( 0.14);
    c->SetTopMargin( 0.16);
    c->SetRightMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->SetLogz();
    c->Update();    
    std::string save_file_path ("invalid"); 
    if (file_type.compare("genie") == 0){ // include in file name that we have set the Mott cross section equal to 1 with the flag "MottXSecEq1"
        save_file_path = "../output/prot_bin_study_elTheta/"+TString(file_type)+"_h2_prot_mom_el_theta_target"+TString(target)+"_beamEnergy"+TString(beam_energy)+"_MottXSecEq1_"+TString(range)+"_elMomLBEq"+TString(elMom_lb)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_protMom"+TString(protMom_lb)+"to"+TString(protMom_ub)+".pdf";
    }
    else if (file_type.compare("data") == 0) { // we can't set data to have a Mott cross section equal to 1, so we remove that flag here
        save_file_path = "../output/prot_bin_study_elTheta/"+TString(file_type)+"_h2_prot_mom_el_theta_target"+TString(target)+"_beamEnergy"+TString(beam_energy)+"_"+TString(range)+"_elMomLBEq"+TString(elMom_lb)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_protMom"+TString(protMom_lb)+"to"+TString(protMom_ub)+".pdf";
    }
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));
 
    // reset the top margin because the bin_selected plots won't have as big titles
    c->SetTopMargin( 0.1);
    c->Update();
    
    // this is the loop to name all of the proton band separated files for GENIE MC
    if (file_type.compare("genie") == 0) {
        for (int i = 0; i < N_BINS; i++) {
            if (N_BINS != 1) { // if N_BINS == 1, then this plot will be the same as the TH2 saved above; we don't want duplicate plots
                bin_selected_h2_prot_mom_el_theta[i]->Draw("colz");
                save_file_path = "../output/prot_bin_study_elTheta/"+TString(file_type)+"_h2_prot_mom_el_theta_target"+target+"_beamEnergy"+beam_energy+"_"+range+"_elMomLBEq"+TString(elMom_lb)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_bin"+elTheta_dot_converter( bin_bounds[i])+"to"+elTheta_dot_converter( bin_bounds[i+1])+"_MottXSecEq1.pdf";
                c->SetLogy( 0); // set y scale to be linear
                c->SaveAs( TString::Format("%s", save_file_path.c_str()));
            }

            bin_selected_h1_el_theta_proj[i]->Draw("HIST");
            save_file_path = "../output/prot_bin_study_elTheta/"+TString(file_type)+"_h1_el_theta_proj_target"+target+"_beamEnergy"+beam_energy+"_"+range+"_elMomLBEq"+TString(elMom_lb)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_bin"+elTheta_dot_converter( bin_bounds[i])+"to"+elTheta_dot_converter( bin_bounds[i+1])+"_MottXSecEq1.pdf";
            c->SetLogy();
            c->SaveAs( TString::Format("%s", save_file_path.c_str()));
        }
    }
    // this is the loop to name all of the proton band separated files for CLAS data
    else if (file_type.compare("data") == 0) {
        for (int i = 0; i < N_BINS; i++) {
            if (N_BINS != 1) { // if N_BINS == 1, then this plot will be the same as the TH2 saved above; we don't want duplicate plots
                bin_selected_h2_prot_mom_el_theta[i]->Draw("colz");
                save_file_path = "../output/prot_bin_study_elTheta/"+TString(file_type)+"_h2_prot_mom_el_theta_target"+target+"_beamEnergy"+beam_energy+"_"+range+"_elMomLBEq"+TString(elMom_lb)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_bin"+elTheta_dot_converter( bin_bounds[i])+"to"+elTheta_dot_converter( bin_bounds[i+1])+".pdf";
                c->SetLogy( 0); // set y scale to be linear
                c->SaveAs( TString::Format("%s", save_file_path.c_str()));
            }

            bin_selected_h1_el_theta_proj[i]->Draw("E1");
            save_file_path = "../output/prot_bin_study_elTheta/"+TString(file_type)+"_h1_el_theta_proj_target"+target+"_beamEnergy"+beam_energy+"_"+range+"_elMomLBEq"+TString(elMom_lb)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_bin"+elTheta_dot_converter( bin_bounds[i])+"to"+elTheta_dot_converter( bin_bounds[i+1])+".pdf";
            c->SetLogy();
            c->SaveAs( TString::Format("%s", save_file_path.c_str()));
        }
    }
}
*/