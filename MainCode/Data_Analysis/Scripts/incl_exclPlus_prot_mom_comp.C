// File to plot the inclusive and exclusive+ case histograms of the 
// proton momentum for a pair of *either* GENIE Monte Carlo (MC) or 
// CLAS data runs. This script puts both the proton momentum plots 
// onto the same canvas.
//
// Author: Jacob Smith (smithja)
// Creation Date: 11/15/2021

// These include statements came with Graham's scripts which I have used
// as the starting point for many of my plotting scripts. Remove them at
// your own peril.
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
//#include "myFunctions.cpp"

#include <algorithm> // smithja

void incl_exclPlus_prot_mom_comp( std::string incl_file_name, std::string exclPlus_file_name){
    std::cout << std::endl; // separate this output message
    std::cout << "/-------------------------------------------------------------\\" << std::endl;
    std::cout << "| CAUTION: MAKE SURE YOU PASS THE INCLUSIVE FILE AS THE FIRST |" << std::endl;
    std::cout << "| ARGUMENT AND EXCLUSIVE+ FILE AS THE SECOND ARGUMENT!!!      |" << std::endl;
    std::cout << "| This message is printed by default as a reminder.           |" << std::endl;
    std::cout << "| -- Jacob Smith (smithja) on November 15th, 2021             |" << std::endl;
    std::cout << "\\-------------------------------------------------------------/" << std::endl;    
    std::cout << std::endl << std::endl; // separate this output message

    // read the file names to make sure we are dealing with either GENIE MC or CLAS data
    auto incl_file_type_start_pos = incl_file_name.rfind("/"); // find the last instance of '/'
    auto incl_file_type_end_pos = incl_file_name.find("_", incl_file_type_start_pos + 1);
    std::string incl_file_type = incl_file_name.substr( incl_file_type_start_pos + 1, incl_file_type_end_pos - (incl_file_type_start_pos + 1));
    
    if (incl_file_type.compare("genie") != 0 && incl_file_type.compare("data") != 0) { // if we do not have GENIE MC or CLAS data as the first file
            std::cout << "/-----------------------------------------------------------------------------------------\\" << std::endl;
            std::cout << "| The first argument passed refers to a file with an invalid name. Please make sure you   |" << std::endl;
            std::cout << "| have one of the keywords 'genie' or 'data' at the beginning of the file and that the    |" << std::endl;
	    std::cout << "| first file argument you have specified only has inclusive cuts applied.                 |" << std::endl;
            std::cout << "| File name example: <PATH_TO_FILE>/<FILE_NAME_WITH_INCLUSIVE_CUTS>.root                  |" << std::endl;
            std::cout << "\\-----------------------------------------------------------------------------------------/" << std::endl;
            return;
    }

    auto exclPlus_file_type_start_pos = exclPlus_file_name.rfind("/"); // find the last instance of '/'
    auto exclPlus_file_type_end_pos = exclPlus_file_name.find("_", exclPlus_file_type_start_pos + 1);
    std::string exclPlus_file_type = exclPlus_file_name.substr( exclPlus_file_type_start_pos + 1, exclPlus_file_type_end_pos - (exclPlus_file_type_start_pos + 1));
    if (exclPlus_file_type.compare( incl_file_type) != 0) { // if the file types do not match
            std::cout << "/-----------------------------------------------------------------------------------------\\" << std::endl;
            std::cout << "| The second argument passed refers to a file with an invalid name. Please make sure you  |" << std::endl;
            std::cout << "| have the keyword at the beginning of the second file match that of the first file.      |" << std::endl;
            std::cout << "\\-----------------------------------------------------------------------------------------/" << std::endl;
            return;
    }

    // declare variables that will be used in the plotting/formatting of histograms and file names
    std::string incl_target ("invalid");
    std::string exclPlus_target ("invalid");
    std::string target ("invalid");
    std::string target_plot ("invalid");
    std::string incl_beam_energy ("invalid");
    std::string exclPlus_beam_energy ("invalid");
    std::string beam_energy ("invalid");
    std::string beam_energy_plot ("invalid");

    std::string elSectors ("invalid");
    std::string elSectors_plot ("invalid");
    std::string deltaPhiEl ("invalid");
    std::string thetaEl_lb ("invalid");
    std::string thetaEl_ub ("invalid");
    std::string elMom_lb ("invalid");
    std::string elMom_lb_plot ("invalid");
    std::string protSectors ("invalid");
    std::string protSectors_plot ("invalid");
    std::string deltaPhiProt ("invalid");
    std::string thetaProt_lb ("invalid");
    std::string thetaProt_ub ("invalid");
    std::string protMom_lb ("invalid");
    std::string protMom_lb_plot ("invalid");

    // import the ROOT data files from the output of doing ./genie_analysis in the e4nu directory
    TFile *input1 = TFile::Open( TString( incl_file_name));
    TFile *input2 = TFile::Open( TString( exclPlus_file_name));

    // Parse the file_name to extract the necessary information that routinely changes for plotting, formatting, etc.
    // Note that the hard-coded numbers for each grouping represent the length of the keyword in question.
    auto incl_target_kwd_pos = incl_file_name.find("target");
    if (incl_target_kwd_pos == std::string::npos) { std::cout << "'target' keyword not found in first file name. Correct the file name!" << std::endl; }
    auto incl_target_end_pos = incl_file_name.find("_", incl_target_kwd_pos);
    incl_target = incl_file_name.substr( incl_target_kwd_pos + 6, incl_target_end_pos - (incl_target_kwd_pos + 6)); // shift by 6 because length of "target" is 6
    auto exclPlus_target_kwd_pos = exclPlus_file_name.find("target");
    if (exclPlus_target_kwd_pos == std::string::npos) { std::cout << "'target' keyword not found in second file name. Correct the file name!" << std::endl; }
    auto exclPlus_target_end_pos = exclPlus_file_name.find("_", exclPlus_target_kwd_pos);
    exclPlus_target = exclPlus_file_name.substr( exclPlus_target_kwd_pos + 6, exclPlus_target_end_pos - (exclPlus_target_kwd_pos + 6)); // shift by 6 because length of "target" is 6
    if (incl_target.compare( exclPlus_target) != 0) {
        std::cout << "ERROR: 'target' keywords in file names do not match!!!" << std::endl;
	std::cout << "Please correct 'target' keywords. No plots created." << std::endl;
        std::cout << "incl_target: " << incl_target << "; exclPlus_target: " << exclPlus_target << std::endl;
        return; // stop the script
    }
    if (incl_target.compare("4He") == 0) { target_plot = "He-4"; }
    else if (incl_target.compare("C12") == 0) { target_plot = "C-12"; }
    else if (incl_target.compare("56Fe") == 0) { target_plot = "Fe-56"; }
    else { std::cout << "You have entered an invalid/unknown target. Consider adding code to accomodate this target in incl_exclPlus_prot_mom_comp.C." << std::endl; }
    //std::cout << incl_target: " << incl_target << ", exclPlus_target: " << exclPlus_target << "; target_plot: " << target_plot << std::endl;

    auto beam1_kwd_pos = incl_file_name.find("beamEnergy");
    if (beam1_kwd_pos == std::string::npos) { std::cout << "'beamEnergy' keyword not found in the first file name. Correct the file name!" << std::endl; }
    auto beam1_end_pos = incl_file_name.find("_", beam1_kwd_pos);
    incl_beam_energy = incl_file_name.substr( beam1_kwd_pos + 10, beam1_end_pos - (beam1_kwd_pos + 10));
    auto beam2_kwd_pos = exclPlus_file_name.find("beamEnergy");
    if (beam2_kwd_pos == std::string::npos) { std::cout << "'beamEnergy' keyword not found in the second file name. Correct the file name!" << std::endl; }
    auto beam2_end_pos = exclPlus_file_name.find("_", beam2_kwd_pos);
    exclPlus_beam_energy = exclPlus_file_name.substr( beam2_kwd_pos + 10, beam2_end_pos - (beam2_kwd_pos + 10));
    if (incl_beam_energy.compare( exclPlus_beam_energy) != 0) {
        std::cout << "ERROR: 'beamEnergy' keywords in file names do not match!!!" << std::endl;
        std::cout << "Please correct 'beamEnergy' keywords. No plots created." << std::endl;
        std::cout << "incl_beam_energy: " << incl_beam_energy << "; exclPlus_beam_energy: " << exclPlus_beam_energy << std::endl;
        return; // stop the script
    }
    if (incl_beam_energy.compare("1161") == 0) { beam_energy_plot = "1.161 GeV"; }
    else if (incl_beam_energy.compare("2261") == 0) { beam_energy_plot = "2.261 GeV"; }
    else if (incl_beam_energy.compare("4461") == 0) { beam_energy_plot = "4.461 GeV"; }
    else { std::cout << "You have entered an invalid/unknown beam energy. Consider adding code to accomodate this beam energy in incl_exclPlus_prot_mom_comp.C" << std::endl; }
    //std::cout << "incl_beam_energy: " << incl_beam_energy << ", exclPlus_beam_energy: " << exclPlus_beam_energy << "; beam_energy_plot: " << beam_energy_plot << std::endl;

    // Future users may want to require that individual (or all) cuts match between
    // the two files used in this script. Rather than propagate the changes to ensure
    // this, we stop at the target and beam energy. Going forward, we alias
    // exclPlus_file_name to just be file_name and exclPlus_file_type to just be file_type. The only 
    // parts of the output this change affects is how the data must be processed and
    // how the final PDF file is saved (i.e. we use the exclusive+ file's cut information
    // since the inclusive cuts should match in the file passed as the first argument). 
    std::string file_name = exclPlus_file_name;
    std::string file_type = exclPlus_file_type;   

    // Since we require that the target and beamEnergy keywords of each file must
    // match, we alias exclPlus_target and exclPlus_beam_energy to just be target and beam_energy
    target = exclPlus_target;
    beam_energy = exclPlus_beam_energy;

    // extract information from the exclusive+ file name to label the output PDF file
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
    }
    //std::cout << "deltaPhiEl: " << deltaPhiEl << std::endl;

    auto thetaEl_kwd_pos = file_name.find("thetaEl");
    if (thetaEl_kwd_pos != std::string::npos){
        auto thetaEl_lb_pos = thetaEl_kwd_pos + 7;
        auto thetaEl_div_pos = file_name.find("to", thetaEl_kwd_pos);
        auto thetaEl_ub_pos = thetaEl_div_pos + 2;
        auto thetaEl_end_pos = file_name.find("_", thetaEl_kwd_pos);
        thetaEl_lb = file_name.substr( thetaEl_lb_pos, thetaEl_div_pos - thetaEl_lb_pos);
        thetaEl_ub = file_name.substr( thetaEl_ub_pos, thetaEl_end_pos - thetaEl_ub_pos);
    }
    //std::cout << "thetaEl_lb: " << thetaEl_lb << "; thetaEl_ub: " << thetaEl_ub << std::endl;

    auto elMom_kwd_pos = file_name.find("elMomLBEq");
    if (elMom_kwd_pos != std::string::npos){
        auto elMom_whl_pos = elMom_kwd_pos + 9;
        auto elMom_div_pos = file_name.find("pt", elMom_kwd_pos);
        auto elMom_dec_pos = elMom_div_pos + 2;
        auto elMom_end_pos = file_name.find("_", elMom_kwd_pos);
        elMom_lb = file_name.substr( elMom_kwd_pos + 9, elMom_end_pos - (elMom_kwd_pos + 9));
        elMom_lb_plot = file_name.substr( elMom_whl_pos, elMom_div_pos - elMom_whl_pos) + "." + file_name.substr( elMom_dec_pos, elMom_end_pos - elMom_dec_pos) + " GeV";
    }
    //std::cout << "elMom_lb: " << elMom_lb << "; elMom_lb_plot: " << elMom_lb_plot << std::endl;

    auto protSectors_kwd_pos = file_name.find("protSectors");
    if (protSectors_kwd_pos != std::string::npos){
        auto protSectors_end_pos = file_name.find("_", protSectors_kwd_pos);
        protSectors = file_name.substr( protSectors_kwd_pos + 11, protSectors_end_pos - (protSectors_kwd_pos + 11));
        unsigned int i;
        for (i = 0; i < protSectors.length() - 1; i++){ protSectors_plot = protSectors_plot + protSectors.substr( i, 1) + " "; }
        protSectors_plot = protSectors_plot + "& " + protSectors.substr( i, 1);
    }
    //std::cout << "protSectors: " << protSectors << "; protSectors_plot: " << protSectors_plot << std::endl;

    auto deltaPhiProtEq_kwd_pos = file_name.find("deltaPhiProtEq");
    if (deltaPhiProtEq_kwd_pos != std::string::npos){
        auto deltaPhiProtEq_end_pos = file_name.find("_", deltaPhiProtEq_kwd_pos);
        deltaPhiProt = file_name.substr( deltaPhiProtEq_kwd_pos + 14, deltaPhiProtEq_end_pos - (deltaPhiProtEq_kwd_pos + 14));
    }
    //std::cout << "deltaPhiProt: " << deltaPhiProt << std::endl;

    auto thetaProt_kwd_pos = file_name.find("thetaProt");
    if (thetaProt_kwd_pos != std::string::npos){
        auto thetaProt_lb_pos = thetaProt_kwd_pos + 9;
        auto thetaProt_div_pos = file_name.find("to", thetaProt_kwd_pos);
        auto thetaProt_ub_pos = thetaProt_div_pos + 2;
        auto thetaProt_end_pos = file_name.find("_", thetaProt_kwd_pos);
        thetaProt_lb = file_name.substr( thetaProt_lb_pos, thetaProt_div_pos - thetaProt_lb_pos);
        thetaProt_ub = file_name.substr( thetaProt_ub_pos, thetaProt_end_pos - thetaProt_ub_pos);
    }
    //std::cout << "thetaProt_lb: " << thetaProt_lb << "; thetaProt_ub: " << thetaProt_ub << std::endl;

    auto protMom_kwd_pos = file_name.find("protMomLBEq");
    if (protMom_kwd_pos != std::string::npos){
        auto protMom_whl_pos = protMom_kwd_pos + 11;
        auto protMom_div_pos = file_name.find("pt", protMom_kwd_pos);
        auto protMom_dec_pos = protMom_div_pos + 2;
        auto protMom_end_pos = file_name.find("_", protMom_kwd_pos);
        protMom_lb = file_name.substr( protMom_kwd_pos + 11, protMom_end_pos - (protMom_kwd_pos + 11));
        protMom_lb_plot = file_name.substr( protMom_whl_pos, protMom_div_pos - protMom_whl_pos) + "." + file_name.substr( protMom_dec_pos, protMom_end_pos - protMom_dec_pos) + " GeV";
    }
    //std::cout << "protMom_lb: " << protMom_lb << "; protMom_lb_plot: " << protMom_lb_plot << std::endl;

    // histogram initialization
    TH1F* h1_prot_momentum_incl_file[6][4];
    TH1F* h1_prot_momentum_exclPlus_file[6][4];

    // extract and compile the GENIE MC histograms, which are separated into sector and interaction
    if (file_type.compare("genie") == 0) {
        for (int i = 0; i < 6; i++) { // for all the sectors
    	    for (int j = 0; j < 4; j++) { // for all the interactions
	        // extract the histrograms
	        h1_prot_momentum_incl_file[i][j] = (TH1F*)input1->Get(TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", j+1, i));
		h1_prot_momentum_exclPlus_file[i][j] = (TH1F*)input2->Get(TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", j+1, i));
	    
	        // compile the histograms by interaction, but not by sector; this is done later because it must be done for both GENIE MC and CLAS data
	        if (j == 0) { continue; } // do not add the first interaction to itself; we are adding everything to the first interaction histogram
	        h1_prot_momentum_incl_file[i][0]->Add( h1_prot_momentum_incl_file[i][j]);
                h1_prot_momentum_exclPlus_file[i][0]->Add( h1_prot_momentum_exclPlus_file[i][j]);
 	    }
	    // compile the interaction-compiled histograms by sector
	    if (i == 0) { continue; } // do not add first sector to itself; we are adding everything to the first sector histogram
	    h1_prot_momentum_incl_file[0][0]->Add( h1_prot_momentum_incl_file[i][0]);
            h1_prot_momentum_exclPlus_file[0][0]->Add( h1_prot_momentum_exclPlus_file[i][0]);
        }
    }
    
    // extract and compile the CLAS data histograms, which are only separated into sectors
    else if (file_type.compare("data") == 0) {
        for (int i = 0; i < 6; i++) { // for all the sectors
	    h1_prot_momentum_incl_file[i][0] = (TH1F*)input1->Get( TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", i));
            h1_prot_momentum_exclPlus_file[i][0] = (TH1F*)input2->Get( TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_prot_mom_QE__%i", i));
        
	    // compile the histograms by sector
	    if (i == 0) { continue; } // do not add first sector to itself; we are adding everything to the first sector histogram
            h1_prot_momentum_incl_file[0][0]->Add( h1_prot_momentum_incl_file[i][0]);
            h1_prot_momentum_exclPlus_file[0][0]->Add( h1_prot_momentum_exclPlus_file[i][0]);
        }
    }

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // format the histograms
    h1_prot_momentum_incl_file[0][0]->Sumw2();
    h1_prot_momentum_incl_file[0][0]->Rebin( 50);
    h1_prot_momentum_incl_file[0][0]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
    h1_prot_momentum_incl_file[0][0]->GetXaxis()->SetTitle("Proton Momentum [GeV/c]");
    h1_prot_momentum_incl_file[0][0]->GetYaxis()->SetTitle("Events");
    h1_prot_momentum_incl_file[0][0]->GetXaxis()->SetRangeUser( 0, 2.0);
    h1_prot_momentum_incl_file[0][0]->GetXaxis()->CenterTitle( true);
    h1_prot_momentum_incl_file[0][0]->GetYaxis()->CenterTitle( true);
//    h1_prot_momentum_incl_file[0][0]->GetXaxis()->SetTitleSize( 0.06);
//    h1_prot_momentum_incl_file[0][0]->GetYaxis()->SetTitleSize( 0.06);
    Double_t incl_file_y_max = h1_prot_momentum_incl_file[0][0]->GetMaximum(); // get y-axis max for plotting purposes later

    h1_prot_momentum_exclPlus_file[0][0]->Sumw2();
    h1_prot_momentum_exclPlus_file[0][0]->Rebin( 50);
    h1_prot_momentum_exclPlus_file[0][0]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
    h1_prot_momentum_exclPlus_file[0][0]->GetXaxis()->SetTitle("Proton Momentum [GeV/c]");
    h1_prot_momentum_exclPlus_file[0][0]->GetYaxis()->SetTitle("Events");
    h1_prot_momentum_exclPlus_file[0][0]->GetXaxis()->SetRangeUser( 0, 2.0);
    h1_prot_momentum_exclPlus_file[0][0]->GetXaxis()->CenterTitle( true);
    h1_prot_momentum_exclPlus_file[0][0]->GetYaxis()->CenterTitle( true);
//    h1_prot_momentum_exclPlus_file[0][0]->GetXaxis()->SetTitleSize( 0.06);
//    h1_prot_momentum_exclPlus_file[0][0]->GetYaxis()->SetTitleSize( 0.06);
    Double_t exclPlus_file_y_max = h1_prot_momentum_exclPlus_file[0][0]->GetMaximum(); // get y-axis max for plotting purposes directly below

    Double_t y_max = std::max( incl_file_y_max, exclPlus_file_y_max);
    h1_prot_momentum_incl_file[0][0]->GetYaxis()->SetRangeUser( 0, 1.1*y_max); // bound plot so we can see all the mc with some headroom
    h1_prot_momentum_exclPlus_file[0][0]->GetYaxis()->SetRangeUser( 0, 1.1*y_max); // bound plot so we can see all the data with some headroom

    gStyle->SetErrorX(0); // don't plot x error bars for any of the distributions

    if (file_type.compare("genie") == 0) {
	h1_prot_momentum_incl_file[0][0]->SetLineColor( kBlack);
        h1_prot_momentum_exclPlus_file[0][0]->SetLineColor( kRed);
    }
    else if (file_type.compare("data") == 0){
        h1_prot_momentum_incl_file[0][0]->SetMarkerStyle( kFullCircle);
        h1_prot_momentum_incl_file[0][0]->SetMarkerColor( kBlack);
        h1_prot_momentum_incl_file[0][0]->SetLineColor( kBlack);
    
        h1_prot_momentum_exclPlus_file[0][0]->SetMarkerStyle( kFullCircle);
        h1_prot_momentum_exclPlus_file[0][0]->SetMarkerColor( kRed);
        h1_prot_momentum_exclPlus_file[0][0]->SetLineColor( kRed);
    }

    // title the plots
    h1_prot_momentum_incl_file[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt)+", "+TString(thetaProt_lb)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub)+"#circ, p_{p^{+}} > "+TString(protMom_lb)+"}");
    h1_prot_momentum_exclPlus_file[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt)+", "+TString(thetaProt_lb)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub)+"#circ, p_{p^{+}} > "+TString(protMom_lb)+"}");

    // draw the plots
    if (file_type.compare("genie") == 0) {
        h1_prot_momentum_incl_file[0][0]->Draw("HIST");
	h1_prot_momentum_exclPlus_file[0][0]->Draw("HIST SAME");
    }
    else if (file_type.compare("data") == 0) {
	h1_prot_momentum_incl_file[0][0]->Draw("E1");
        h1_prot_momentum_exclPlus_file[0][0]->Draw("E1 SAME");
    }
    
    // draw a legend for our plot
    std::string file_type_leg = file_type; // alias file_type so we can include information about if we are running over GENIE MC or CLAS data in the legend
    std::transform( file_type_leg.begin(), file_type_leg.end(), file_type_leg.begin(), ::toupper);
    TLegend *legend = new TLegend( 0.675, 0.725, 0.875, 0.825);
    legend->AddEntry( h1_prot_momentum_incl_file[0][0], "Inclusive, "+TString( file_type_leg));
    legend->AddEntry( h1_prot_momentum_exclPlus_file[0][0], "Exclusive+, "+TString( file_type_leg));
    legend->SetBorderSize( 0);
    legend->SetFillStyle( 0);
    legend->Draw();
    
    // crop the margins of the canvas
    c->SetLeftMargin( 0.14);
    c->SetBottomMargin( 0.14);
    c->SetTopMargin( 0.14);
    c->Update();
    
    // save this histogram as a PDF file 
    std::string save_file_path = "../output/incl_exclPlus_prot_mom_comp/"+file_type+"_incl_exclPlus_prot_mom_target"+target+"_beamEnergy"+beam_energy+"_MottXSecEq1_elSectors"+elSectors+"_deltaPhiEq"+deltaPhiEl+"_thetaEl"+thetaEl_lb+"to"+thetaEl_ub+"_elMomLBEq"+elMom_lb+"_protSectors"+protSectors+"_deltaPhiProtEq"+deltaPhiProt+"_thetaProt"+thetaProt_lb+"to"+thetaProt_ub+"_protMomLBEq"+protMom_lb+".pdf"; 
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));
    
    // smithja: this section of code integrates the distribution and prints it to the screen
    int num_events_incl_file = h1_prot_momentum_incl_file[0][0]->Integral();
    int num_events_exclPlus_file = h1_prot_momentum_exclPlus_file[0][0]->Integral();
    int num_events = num_events_incl_file + num_events_exclPlus_file;

    cout << "// --- PROGRAM OUTPUT -----------------------------------------------------------" << endl;
    cout << "There are " << num_events_incl_file << " many events in the inclusive file's h1_prot_momentum histogram." << endl;
    cout << "There are " << num_events_exclPlus_file << " many events in the exclusive+ file's h1_prot_momentum histogram." << endl;
    cout << "There are " << num_events << " total events in both of the files' h1_prot_momentum histograms." << endl;
    cout << "// ------------------------------------------------------------------------------" << endl << endl; // separate program output
}
