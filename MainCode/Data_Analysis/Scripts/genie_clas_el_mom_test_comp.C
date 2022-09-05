// File to plot histograms of the electron momentum for a pair of GENIE Monte
// Carlo (MC) and CLAS data runs.
//
// NOTE: At this time (11/10/2021), I am requiring that the target
//       and beam energy be equal between the GENIE MC file and
//       CLAS data file. Future users (this includes my future self)
//       may want to remove this functionality or explicitly require
//       that all the kinematic cuts be the same. To do away with
//       target matching and/or beam energy matching, I (smithja) 
//       suggest adding an argument to the function call signature
//       to override this as necessary. To explicitly match some/all
//       of the kinematic cuts, you will have to restructure the
//       code quite a bit. 
//
// NOTE: In instances where we extract target, beam energy, and 
//       kinematic cut information from the input file names to use
//       in plotting or saving the histograms as PDF files, we
//       always use the GENIE MC file name. 
//
// Author: Jacob Smith (smithja)
// Creation Date: 11/10/2021

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


// This auxillary function takes in a string (which should be the
// written out version of a double) with "pt" in its name and
// converts the string to the applicable double value.
// Ex.) (std::string) 1pt25 --> (double) 1.25
double genie_clas_el_mom_test_comp_pt_converter( std::string str) {
    if (str.find("pt") != std::string::npos) { str = std::regex_replace( str, std::regex("pt"), "."); } // replace 'pt' -> '.'
    double val = stod( str);
    return val;
}

// -----------------------------------------------------------------------------------------------------------------------
// This auxillary function utilizes genie_clas_el_mom_test_comp_pt_converter to take the resulting
// value from genie_clas_el_mom_test_comp_pt_converter, remove any trailing zeros, and then recast
// the double from genie_clas_el_mom_test_comp_pt_converter to be a std::string.
std::string genie_clas_el_mom_test_comp_pt_plot_converter( std::string str) {
    str = std::to_string( genie_clas_el_mom_test_comp_pt_converter( str)); // replace 'pt' -> '.' and recast as a std::string
    while (str.substr( str.length() - 1).compare("0") == 0){ // while the last character in str is a 0
        str = str.substr( 0, str.length() - 1); // take off the trailing 0
    }
    if (str.substr( str.length() - 1).compare(".") == 0) { // if the last character is a decimal point with no numbers after it
        str = str.substr( 0, str.length() - 1); // take off the decimal point
    }
    return str.c_str();
}

// MAIN FUNCTION ---------------------------------------------------------------------------------------------------------
void genie_clas_el_mom_test_comp( std::string mc_file_name, std::string data_file_name){
    // This COMMONLY EDITED SECTION contains variables that a user may want to routinely modify.
    Int_t rebin_factor = 10; // how many bins to merge together into one for the final histogram; h1_el_mom_test plots created by genie_analysis range from 0 to 6 and have 6000 bins

    // End of COMMONLY EDITED SECTION --------------

    // read the file names to make sure we are dealing with either GENIE MC or CLAS data
    // and ensure that the files are passed in the right order when running this script

    auto mc_file_type_start_pos = mc_file_name.rfind("/"); // find the last instance of '/'
    auto mc_file_type_end_pos = mc_file_name.find("_", mc_file_type_start_pos + 1);
    std::string mc_file_type = mc_file_name.substr( mc_file_type_start_pos + 1, mc_file_type_end_pos - (mc_file_type_start_pos + 1));
    if (mc_file_type.compare("genie") != 0) { // if we do not have GENIE MC as the first file
            std::cout << "/-----------------------------------------------------------------------------------------\\" << std::endl;
            std::cout << "| The first argument passed refers to a file with an invalid name. Please make sure you   |" << std::endl;
            std::cout << "| have the keyword 'genie' at the beginning of the file and that this file is GENIE data. |" << std::endl;
            std::cout << "| File name example: <PATH_TO_FILE>/genie_<OMITTED>.root                                  |" << std::endl;
            std::cout << "\\-----------------------------------------------------------------------------------------/" << std::endl;
            return;
    }

    auto data_file_type_start_pos = data_file_name.rfind("/"); // find the last instance of '/'
    auto data_file_type_end_pos = data_file_name.find("_", data_file_type_start_pos + 1);
    std::string data_file_type = data_file_name.substr( data_file_type_start_pos + 1, data_file_type_end_pos - (data_file_type_start_pos + 1));
    if (data_file_type.compare("data") != 0) { // if we do not have CLAS data as the second file
            std::cout << "/-----------------------------------------------------------------------------------------\\" << std::endl;
            std::cout << "| The second argument passed refers to a file with an invalid name. Please make sure you  |" << std::endl;
            std::cout << "| have the keyword 'data' at the beginning of the file and that this file is CLAS data.   |" << std::endl;
            std::cout << "| File name example: <PATH_TO_FILE>/data_<OMITTED>.root.                                  |" << std::endl;
            std::cout << "\\-----------------------------------------------------------------------------------------/" << std::endl;
            return;
    }

    // declare variables that will be used in the plotting/formatting of histograms and file names
    std::string mc_target ("invalid");
    std::string data_target ("invalid");
    std::string target_plot ("invalid");
    std::string mc_beam_energy ("invalid");
    std::string data_beam_energy ("invalid");
    std::string beam_energy_plot ("invalid");

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

    // import the ROOT data files from the output of doing ./genie_analysis in the e4nu directory
    TFile *mc_input = TFile::Open( TString( mc_file_name));
    TFile *data_input = TFile::Open( TString( data_file_name));

    // Parse the file_name to extract the necessary information that routinely changes for plotting, formatting, etc.
    // Note that the hard-coded numbers for each grouping represent the length of the keyword in question.
    auto mc_target_kwd_pos = mc_file_name.find("target");
    if (mc_target_kwd_pos == std::string::npos) { std::cout << "'target' keyword not found in MC file_name. Correct file_name!" << std::endl; }
    auto mc_target_end_pos = mc_file_name.find("_", mc_target_kwd_pos);
    mc_target = mc_file_name.substr( mc_target_kwd_pos + 6, mc_target_end_pos - (mc_target_kwd_pos + 6)); // shift by 6 because length of "target" is 6
    auto data_target_kwd_pos = data_file_name.find("target");
    if (data_target_kwd_pos == std::string::npos) { std::cout << "'target' keyword not found in CLAS data file_name. Correct file_name!" << std::endl; }
    auto data_target_end_pos = data_file_name.find("_", data_target_kwd_pos);
    data_target = data_file_name.substr( data_target_kwd_pos + 6, data_target_end_pos - (data_target_kwd_pos + 6)); // shift by 6 because length of "target" is 6
    if (mc_target.compare( data_target) != 0) {
        std::cout << "ERROR: 'target' keywords in file names do not match!!!" << std::endl;
	std::cout << "Please correct 'target' keywords. No plots created." << std::endl;
        std::cout << "mc_target: " << mc_target << "; data_target: " << data_target << std::endl;
        return; // stop the script
    }
    if (mc_target.compare("4He") == 0) { target_plot = "He-4"; }
    else if (mc_target.compare("C12") == 0) { target_plot = "C-12"; }
    else if (mc_target.compare("56Fe") == 0) { target_plot = "Fe-56"; }
    else { std::cout << "You have entered an invalid/unknown target. Consider adding code to accomodate this target in genie_clas_el_mom_test_comp.C." << std::endl; }
    //std::cout << "mc_target: " << mc_target << "; data_target: " << data_target << "; target_plot: " << target_plot << std::endl;

    auto mc_beam_kwd_pos = mc_file_name.find("beamEnergy");
    if (mc_beam_kwd_pos == std::string::npos) { std::cout << "'beamEnergy' keyword not found in the MC file_name. Correct the MC file_name!" << std::endl; }
    auto mc_beam_end_pos = mc_file_name.find("_", mc_beam_kwd_pos);
    mc_beam_energy = mc_file_name.substr( mc_beam_kwd_pos + 10, mc_beam_end_pos - (mc_beam_kwd_pos + 10));
    auto data_beam_kwd_pos = data_file_name.find("beamEnergy");
    if (data_beam_kwd_pos == std::string::npos) { std::cout << "'beamEnergy' keyword not found in the CLAS data file_name. Correct the CLAS file_name!" << std::endl; }
    auto data_beam_end_pos = data_file_name.find("_", data_beam_kwd_pos);
    data_beam_energy = data_file_name.substr( data_beam_kwd_pos + 10, data_beam_end_pos - (data_beam_kwd_pos + 10));
    if (mc_beam_energy.compare( data_beam_energy) != 0) {
        std::cout << "ERROR: 'beamEnergy' keywords in file names do not match!!!" << std::endl;
        std::cout << "Please correct 'beamEnergy' keywords. No plots created." << std::endl;
        std::cout << "mc_beam_energy: " << mc_beam_energy << "; data_beam_energy: " << data_beam_energy << std::endl;
        return; // stop the script
    }
    if (mc_beam_energy.compare("1161") == 0) { beam_energy_plot = "1.161 GeV"; }
    else if (mc_beam_energy.compare("2261") == 0) { beam_energy_plot = "2.261 GeV"; }
    else if (mc_beam_energy.compare("4461") == 0) { beam_energy_plot = "4.461 GeV"; }
    else { std::cout << "You have entered an invalid/unknown beam energy. Consider adding code to accomodate this beam energy in genie_clas_el_mom_test_comp.C" << std::endl; }
    //std::cout << "mc_beam_energy: " << mc_beam_energy << "data_beam_energy: " << data_beam_energy << "; beam_energy_plot: " << beam_energy_plot << std::endl;

    auto elSectors_kwd_pos = mc_file_name.find("elSectors");
    if (elSectors_kwd_pos != std::string::npos){
        auto elSectors_end_pos = mc_file_name.find("_", elSectors_kwd_pos);
        elSectors = mc_file_name.substr( elSectors_kwd_pos + 9, elSectors_end_pos - (elSectors_kwd_pos + 9));
        unsigned int i;
        elSectors_plot = ""; // elSectors is valid so first do away with 'invalid'
        for (i = 0; i < elSectors.length() - 1; i++){ elSectors_plot = elSectors_plot + elSectors.substr( i, 1) + " "; }
        elSectors_plot = elSectors_plot + "& " + elSectors.substr( i, 1);
    }
    //std::cout << "elSectors: " << elSectors << "; elSectors_plot: " << elSectors_plot << std::endl;

    auto deltaPhiElEq_kwd_pos = mc_file_name.find("deltaPhiElEq");
    if (deltaPhiElEq_kwd_pos != std::string::npos){
        auto deltaPhiElEq_end_pos = mc_file_name.find("_", deltaPhiElEq_kwd_pos);
        deltaPhiEl = mc_file_name.substr( deltaPhiElEq_kwd_pos + 12, deltaPhiElEq_end_pos - (deltaPhiElEq_kwd_pos + 12));
        deltaPhiEl_plot = genie_clas_el_mom_test_comp_pt_plot_converter( deltaPhiEl);
    }
    //std::cout << "deltaPhiEl: " << deltaPhiEl << std::endl;

    auto thetaEl_kwd_pos = mc_file_name.find("thetaEl");
    if (thetaEl_kwd_pos != std::string::npos){
        auto thetaEl_lb_whl_pos = thetaEl_kwd_pos + 7;
        auto thetaEl_lb_div_pos = mc_file_name.find("pt", thetaEl_lb_whl_pos);
        auto thetaEl_lb_dec_pos = thetaEl_lb_div_pos + 2;
        auto thetaEl_div_pos = mc_file_name.find("to", thetaEl_lb_dec_pos);
        auto thetaEl_ub_whl_pos = thetaEl_div_pos + 2;
        auto thetaEl_ub_div_pos = mc_file_name.find("pt", thetaEl_ub_whl_pos);
        auto thetaEl_ub_dec_pos = thetaEl_ub_div_pos + 2;
        auto thetaEl_end_pos = mc_file_name.find("_", thetaEl_ub_dec_pos);
        thetaEl_lb = mc_file_name.substr( thetaEl_lb_whl_pos, thetaEl_div_pos - thetaEl_lb_whl_pos);
        thetaEl_ub = mc_file_name.substr( thetaEl_ub_whl_pos, thetaEl_end_pos - thetaEl_ub_whl_pos);
        thetaEl_lb_plot = genie_clas_el_mom_test_comp_pt_plot_converter( thetaEl_lb);
        thetaEl_ub_plot = genie_clas_el_mom_test_comp_pt_plot_converter( thetaEl_ub);
    }
    //std::cout << "thetaEl_lb: " << thetaEl_lb << ", thetaEl_lb_plot: " << thetaEl_lb_plot << "; thetaEl_ub: " << thetaEl_ub << ", thetaEl_ub_plot: " << thetaEl_ub_plot << std::endl;

    auto elMom_kwd_pos = mc_file_name.find("elMomLBEq");
    if (elMom_kwd_pos != std::string::npos){
        auto elMom_whl_pos = elMom_kwd_pos + 9;
        auto elMom_div_pos = mc_file_name.find("pt", elMom_kwd_pos);
        auto elMom_dec_pos = elMom_div_pos + 2;
        auto elMom_end_pos = mc_file_name.find("_", elMom_kwd_pos);
        elMom_lb = mc_file_name.substr( elMom_kwd_pos + 9, elMom_end_pos - (elMom_kwd_pos + 9));
        elMom_lb_plot = mc_file_name.substr( elMom_whl_pos, elMom_div_pos - elMom_whl_pos) + "." + mc_file_name.substr( elMom_dec_pos, elMom_end_pos - elMom_dec_pos) + " GeV";
    }
    //std::cout << "elMom_lb: " << elMom_lb << "; elMom_lb_plot: " << elMom_lb_plot << std::endl;

    auto protSectors_kwd_pos = mc_file_name.find("protSectors");
    if (protSectors_kwd_pos != std::string::npos){
        auto protSectors_end_pos = mc_file_name.find("_", protSectors_kwd_pos);
        protSectors = mc_file_name.substr( protSectors_kwd_pos + 11, protSectors_end_pos - (protSectors_kwd_pos + 11));
        unsigned int i;
        protSectors_plot = ""; // protSectors is valid so first do away with 'invalid'
        for (i = 0; i < protSectors.length() - 1; i++){ protSectors_plot = protSectors_plot + protSectors.substr( i, 1) + " "; }
        protSectors_plot = protSectors_plot + "& " + protSectors.substr( i, 1);
    }
    //std::cout << "protSectors: " << protSectors << "; protSectors_plot: " << protSectors_plot << std::endl;

    auto deltaPhiProtEq_kwd_pos = mc_file_name.find("deltaPhiProtEq");
    if (deltaPhiProtEq_kwd_pos != std::string::npos){
        auto deltaPhiProtEq_end_pos = mc_file_name.find("_", deltaPhiProtEq_kwd_pos);
        deltaPhiProt = mc_file_name.substr( deltaPhiProtEq_kwd_pos + 14, deltaPhiProtEq_end_pos - (deltaPhiProtEq_kwd_pos + 14));
        deltaPhiProt_plot = genie_clas_el_mom_test_comp_pt_plot_converter( deltaPhiProt);
    }
    //std::cout << "deltaPhiProt: " << deltaPhiProt << std::endl;

    auto thetaProt_kwd_pos = mc_file_name.find("thetaProt");
    if (thetaProt_kwd_pos != std::string::npos){
        auto thetaProt_lb_whl_pos = thetaProt_kwd_pos + 9;
        auto thetaProt_lb_div_pos = mc_file_name.find("pt", thetaProt_lb_whl_pos);
        auto thetaProt_lb_dec_pos = thetaProt_lb_div_pos + 2;
        auto thetaProt_div_pos = mc_file_name.find("to", thetaProt_lb_dec_pos);
        auto thetaProt_ub_whl_pos = thetaProt_div_pos + 2;
        auto thetaProt_ub_div_pos = mc_file_name.find("pt", thetaProt_ub_whl_pos);
        auto thetaProt_ub_dec_pos = thetaProt_ub_div_pos + 2;
        auto thetaProt_end_pos = mc_file_name.find("_", thetaProt_ub_dec_pos);
        thetaProt_lb = mc_file_name.substr( thetaProt_lb_whl_pos, thetaProt_div_pos - thetaProt_lb_whl_pos);
        thetaProt_ub = mc_file_name.substr( thetaProt_ub_whl_pos, thetaProt_end_pos - thetaProt_ub_whl_pos);
        thetaProt_lb_plot = genie_clas_el_mom_test_comp_pt_plot_converter( thetaProt_lb);
        thetaProt_ub_plot = genie_clas_el_mom_test_comp_pt_plot_converter( thetaProt_ub);
    }
    //std::cout << "thetaProt_lb: " << thetaProt_lb << ", thetaProt_lb_plot: " << thetaProt_lb_plot << "; thetaProt_ub: " << thetaProt_ub << ", thetaProt_ub_plot: " << thetaProt_ub_plot << std::endl;

    auto protMom_lb_kwd_pos = mc_file_name.find("protMomLBEq");
    auto protMom_ub_kwd_pos = mc_file_name.find("protMomUBEq");
    auto protMom_kwd_pos = mc_file_name.find("protMom");
    if (protMom_lb_kwd_pos != std::string::npos){
        auto protMom_whl_pos = protMom_lb_kwd_pos + 11;
        auto protMom_div_pos = mc_file_name.find("pt", protMom_lb_kwd_pos);
        auto protMom_dec_pos = protMom_div_pos + 2;
        auto protMom_end_pos = mc_file_name.find("_", protMom_lb_kwd_pos);
        protMom_lb = mc_file_name.substr( protMom_lb_kwd_pos + 11, protMom_end_pos - (protMom_lb_kwd_pos + 11));
        protMom_lb_plot = mc_file_name.substr( protMom_whl_pos, protMom_div_pos - protMom_whl_pos) + "." + mc_file_name.substr( protMom_dec_pos, protMom_end_pos - protMom_dec_pos) + " GeV";
    }
     //std::cout << "protMom_lb: " << protMom_lb << "; protMom_lb_plot: " << protMom_lb_plot << std::endl;
    else if (protMom_ub_kwd_pos != std::string::npos){
        auto protMom_whl_pos = protMom_ub_kwd_pos + 11;
        auto protMom_div_pos = mc_file_name.find("pt", protMom_ub_kwd_pos);
        auto protMom_dec_pos = protMom_div_pos + 2;
        auto protMom_end_pos = mc_file_name.find("_", protMom_ub_kwd_pos);
        protMom_ub = mc_file_name.substr( protMom_ub_kwd_pos + 11, protMom_end_pos - (protMom_ub_kwd_pos + 11));
        protMom_ub_plot = mc_file_name.substr( protMom_whl_pos, protMom_div_pos - protMom_whl_pos) + "." + mc_file_name.substr( protMom_dec_pos, protMom_end_pos - protMom_dec_pos) + " GeV";
    }
    //std::cout << "protMom_ub: " << protMom_ub << "; protMom_ub_plot: " << protMom_ub_plot << std::endl;
    else if (protMom_kwd_pos != std::string::npos){
        auto protMom_lb_whl_pos = protMom_kwd_pos + 7;
        auto protMom_lb_div_pos = mc_file_name.find("pt", protMom_lb_whl_pos);
        auto protMom_lb_dec_pos = protMom_lb_div_pos + 2;
        auto protMom_div_pos = mc_file_name.find("to", protMom_lb_dec_pos);
        auto protMom_ub_whl_pos = protMom_div_pos + 2;
        auto protMom_ub_div_pos = mc_file_name.find("pt", protMom_ub_whl_pos);
        auto protMom_ub_dec_pos = protMom_ub_div_pos + 2;
        auto protMom_end_pos = mc_file_name.find("_", protMom_ub_dec_pos);
        protMom_lb = mc_file_name.substr( protMom_lb_whl_pos, protMom_div_pos - protMom_lb_whl_pos);
        protMom_ub = mc_file_name.substr( protMom_ub_whl_pos, protMom_end_pos - protMom_ub_whl_pos);
        protMom_lb_plot = genie_clas_el_mom_test_comp_pt_plot_converter( protMom_lb);
        protMom_ub_plot = genie_clas_el_mom_test_comp_pt_plot_converter( protMom_ub);
    }
    //std::cout << "protMom_lb: " << protMom_lb << ", protMom_lb_plot: " << protMom_lb_plot << "; protMom_ub: " << protMom_ub << ", protMom_ub_plot: " << protMom_ub_plot << std::endl;

    // Determine what case we are in (i.e. no cuts, inclusive w/o electron momentum cut,
    // inclusive, exclusive, or exclusive with protMom cut) based on the null and non-null
    // values of the above pointers.
    std::string cut_case ("invalid");
    if (mc_file_name.find("noCuts") != std::string::npos) { cut_case = "noCuts"; }
    if ((elSectors_kwd_pos != std::string::npos) && (deltaPhiElEq_kwd_pos != std::string::npos) && (thetaEl_kwd_pos != std::string::npos)) { cut_case = "incl_noElMomCut"; }
    if ((cut_case.compare("incl_noElMomCut") == 0) && (elMom_kwd_pos != std::string::npos)) { cut_case = "incl"; }
    if ((cut_case.compare("incl") == 0) && (protSectors_kwd_pos != std::string::npos) && (deltaPhiProtEq_kwd_pos != std::string::npos) && (thetaProt_kwd_pos != std::string::npos)) { cut_case = "excl"; }
    if (cut_case.compare("excl") == 0) {
        if (protMom_lb_kwd_pos != std::string::npos) { cut_case = "excl_withProtMomLB"; }
        if (protMom_ub_kwd_pos != std::string::npos) { cut_case = "excl_withProtMomUB"; }
        if (protMom_kwd_pos != std::string::npos) { cut_case = "excl_withProtMomSlice"; }
    }

    // histogram initialization
    TH1F* mc_h1_el_mom_test_momentum[6][4];
    TH1F* data_h1_el_mom_test_momentum[6];

    // extract and compile the GENIE MC histograms, which are separated into sector and interaction
    for (int i = 0; i < 6; i++) { // for all the sectors
	for (int j = 0; j < 4; j++) { // for all the interactions
	    // extract the histrograms
	    mc_h1_el_mom_test_momentum[i][j] = (TH1F*)mc_input->Get(TString::Format("h1_%i_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__%i", j+1, i));
	    
	    // compile the histograms by interaction, but not by sector; this is done later because it must be done for both GENIE MC and CLAS data
	    if (j == 0) { continue; } // do not add the first interaction to itself; we are adding everything to the first interaction histogram
	    mc_h1_el_mom_test_momentum[i][0]->Add( mc_h1_el_mom_test_momentum[i][j]);
	}
	
    }
    
    // extract and compile the CLAS data histograms, which are only separated into sect
    for (int i = 0; i < 6; i++) { // for all the sectors
	data_h1_el_mom_test_momentum[i] = (TH1F*)data_input->Get( TString::Format("h1_0_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_el_mom__%i", i));
    }

    // compile all the sectors into one histogram
    for( int i = 1; i < 6; i++) { 
	mc_h1_el_mom_test_momentum[0][0]->Add( mc_h1_el_mom_test_momentum[i][0]);
	data_h1_el_mom_test_momentum[0]->Add( data_h1_el_mom_test_momentum[i]); 
    }

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // format the histograms
    mc_h1_el_mom_test_momentum[0][0]->Sumw2();
    mc_h1_el_mom_test_momentum[0][0]->Rebin( rebin_factor);
    mc_h1_el_mom_test_momentum[0][0]->Scale( data_h1_el_mom_test_momentum[0]->Integral()/mc_h1_el_mom_test_momentum[0][0]->Integral()); // normalize MC to CLAS data
    mc_h1_el_mom_test_momentum[0][0]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
    mc_h1_el_mom_test_momentum[0][0]->GetXaxis()->SetTitle("Electron Momentum [GeV/c]");
    mc_h1_el_mom_test_momentum[0][0]->GetYaxis()->SetTitle("Events");
    mc_h1_el_mom_test_momentum[0][0]->GetXaxis()->SetRangeUser( 1.99, 2.2);
    mc_h1_el_mom_test_momentum[0][0]->GetXaxis()->CenterTitle( true);
    mc_h1_el_mom_test_momentum[0][0]->GetYaxis()->CenterTitle( true);
    Double_t mc_y_max = mc_h1_el_mom_test_momentum[0][0]->GetMaximum(); // get y-axis max for plotting purposes later

    data_h1_el_mom_test_momentum[0]->Sumw2();
    data_h1_el_mom_test_momentum[0]->Rebin( rebin_factor);
    data_h1_el_mom_test_momentum[0]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
    data_h1_el_mom_test_momentum[0]->GetXaxis()->SetTitle("Electron Momentum [GeV/c]");
    data_h1_el_mom_test_momentum[0]->GetYaxis()->SetTitle("Events");
    data_h1_el_mom_test_momentum[0]->GetXaxis()->SetRangeUser( 1.99, 2.2);  
    data_h1_el_mom_test_momentum[0]->GetXaxis()->CenterTitle( true);
    data_h1_el_mom_test_momentum[0]->GetYaxis()->CenterTitle( true);
    data_h1_el_mom_test_momentum[0]->SetMarkerStyle( kFullCircle);
    data_h1_el_mom_test_momentum[0]->SetMarkerColor( kBlack);
    data_h1_el_mom_test_momentum[0]->SetLineColor( kBlack);
    Double_t data_y_max = data_h1_el_mom_test_momentum[0]->GetMaximum(); // get y-axis max for plotting purposes directly below

    Double_t mc_data_y_max = std::max( mc_y_max, data_y_max);
    mc_h1_el_mom_test_momentum[0][0]->GetYaxis()->SetRangeUser( 0, 1.1*mc_data_y_max); // bound plot so we can see all the mc with some headroom
    data_h1_el_mom_test_momentum[0]->GetYaxis()->SetRangeUser( 0, 1.1*mc_data_y_max); // bound plot so we can see all the data with some headroom

    gStyle->SetErrorX(0); // don't plot x error bars for any of the distributions

    if (cut_case.compare("noCuts") == 0) { 
	mc_h1_el_mom_test_momentum[0][0]->SetTitle(TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam (No Cuts)"); 
        data_h1_el_mom_test_momentum[0]->SetTitle(TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam (No Cuts)"); 
    }
    else if ( cut_case.compare("incl_noElMomCut") == 0) { 
	mc_h1_el_mom_test_momentum[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ}"); 
        data_h1_el_mom_test_momentum[0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ}"); 
    }
    else if ( cut_case.compare("incl") == 0) { 
	mc_h1_el_mom_test_momentum[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+"}"); 
        data_h1_el_mom_test_momentum[0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+"}");
    }
    else if ( cut_case.compare("excl") == 0) { 
	mc_h1_el_mom_test_momentum[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ}");
        data_h1_el_mom_test_momentum[0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ}");
    }
    else if ( cut_case.compare("excl_withProtMomLB") == 0) { 
	mc_h1_el_mom_test_momentum[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, p_{p^{+}} > "+TString(protMom_lb_plot)+"}"); 
        data_h1_el_mom_test_momentum[0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, p_{p^{+}} > "+TString(protMom_lb_plot)+"}");
    }
    else if ( cut_case.compare("excl_withProtMomUB") == 0) {
        mc_h1_el_mom_test_momentum[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, p_{p^{+}} < "+TString(protMom_ub_plot)+"}");
        data_h1_el_mom_test_momentum[0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, p_{p^{+}} < "+TString(protMom_ub_plot)+"}");
    }
    else if ( cut_case.compare("excl_withProtMomSlice") == 0) {
        mc_h1_el_mom_test_momentum[0][0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, "+TString(protMom_lb_plot)+" < p_{p^{+}} < "+TString(protMom_ub_plot)+"}");
        data_h1_el_mom_test_momentum[0]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl_plot)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb_plot)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub_plot)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt_plot)+", "+TString(thetaProt_lb_plot)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub_plot)+"#circ, "+TString(protMom_lb_plot)+" < p_{p^{+}} < "+TString(protMom_ub_plot)+"}");
    }

    mc_h1_el_mom_test_momentum[0][0]->Draw("HIST");
    data_h1_el_mom_test_momentum[0]->Draw("E1 SAME");

    // draw a legend for our plot
    TLegend *legend = new TLegend( 0.2, 0.725, 0.35, 0.825);
    legend->AddEntry( mc_h1_el_mom_test_momentum[0][0], "GENIE MC");
    legend->AddEntry( data_h1_el_mom_test_momentum[0], "CLAS data");
    legend->SetBorderSize( 0);
    legend->SetFillStyle( 0);
    legend->Draw();

    // crop the margins of the canvas
    c->SetLeftMargin( 0.14);
    c->SetBottomMargin( 0.14);
    if ( cut_case.compare("noCuts") != 0) { c->SetTopMargin( 0.16); }
    c->Update();

    // save this histogram as a PDF file
    std::string save_file_path ("invalid");
    if ( cut_case.compare("noCuts") == 0) { save_file_path = "../output/genie_clas_el_mom_test_comp/genie_clas_comparison_h1_el_mom_test_target"+mc_target+"_beamEnergy"+mc_beam_energy+"_MottXSecEq1_noCuts.pdf"; }
    else if ( cut_case.compare("incl_noElMomCut") == 0) { save_file_path = "../output/genie_clas_el_mom_test_comp/genie_clas_comparison_h1_el_mom_test_target"+TString(mc_target)+"_beamEnergy"+TString(mc_beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+".pdf"; }
    else if ( cut_case.compare("incl") == 0) { save_file_path = "../output/genie_clas_el_mom_test_comp/genie_clas_comparison_h1_el_mom_test_target"+TString(mc_target)+"_beamEnergy"+TString(mc_beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+".pdf"; }
    else if ( cut_case.compare("excl") == 0) { save_file_path = "../output/genie_clas_el_mom_test_comp/genie_clas_comparison_h1_el_mom_test_target"+TString(mc_target)+"_beamEnergy"+TString(mc_beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+"_protSectors"+TString(protSectors)+"_deltaPhiProtEq"+TString(deltaPhiProt)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+".pdf"; }
    else if ( cut_case.compare("excl_withProtMomLB") == 0) { save_file_path = "../output/genie_clas_el_mom_test_comp/genie_clas_comparison_h1_el_mom_test_target"+TString(mc_target)+"_beamEnergy"+TString(mc_beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+"_protSectors"+TString(protSectors)+"_deltaPhiProtEq"+TString(deltaPhiProt)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_protMomLBEq"+TString(protMom_lb)+".pdf"; }
    else if ( cut_case.compare("excl_withProtMomUB") == 0) { save_file_path = "../output/genie_clas_el_mom_test_comp/genie_clas_comparison_h1_el_mom_test_target"+TString(mc_target)+"_beamEnergy"+TString(mc_beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+"_protSectors"+TString(protSectors)+"_deltaPhiProtEq"+TString(deltaPhiProt)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_protMomUBEq"+TString(protMom_ub)+".pdf"; }
    else if ( cut_case.compare("excl_withProtMomSlice") == 0) { save_file_path = "../output/genie_clas_el_mom_test_comp/genie_clas_comparison_h1_el_mom_test_target"+TString(mc_target)+"_beamEnergy"+TString(mc_beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+"_protSectors"+TString(protSectors)+"_deltaPhiProtEq"+TString(deltaPhiProt)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_protMom"+TString(protMom_lb)+"to"+TString(protMom_ub)+".pdf"; }
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));


    // smithja: this section of code integrates the distribution and prints it to the screen
    int num_events_mc = mc_h1_el_mom_test_momentum[0][0]->Integral();
    int num_events_data = data_h1_el_mom_test_momentum[0]->Integral();
    int num_events = num_events_mc + num_events_data;

    cout << "// --- PROGRAM OUTPUT -----------------------------------------------------------" << endl;
    cout << "There are " << num_events_mc << " many events in the GENIE MC h1_el_mom_test_momentum histogram." << endl;
    cout << "There are " << num_events_data << " many events in the CLAS data h1_el_mom_test_momentum histogram." << endl;
    cout << "There are " << num_events << " total events in both of the GENIE and CLAS h1_el_mom_test_momentum histograms." << endl;
    cout << "// ------------------------------------------------------------------------------" << endl << endl; // separate program output
}
