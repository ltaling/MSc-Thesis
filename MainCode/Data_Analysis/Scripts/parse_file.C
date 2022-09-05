#include <TFile.h>
#include <TList.h>
#include <TParameter.h>
#include <TNamed.h>

#include <iostream>

// -----------------------------------------------------------------------------------------------------------------------

std::string double_converter( double x) {

    std::string str = std::to_string(x);

    while (str.substr( str.length() - 1).compare("0") == 0){ // while the last character in str is a 0
        str = str.substr( 0, str.length() - 1); // take off the trailing 0
    }
    if (str.substr( str.length() - 1).compare(".") == 0) { // if the last character is a decimal point with no numbers after it
        str = str.substr( 0, str.length() - 1); // take off the decimal point
    }
    return str.c_str();
}

// -----------------------------------------------------------------------------------------------------------------------

std::vector<std::string> parse_file(TList *Run_Info) {

	std::vector<std::string> out;

	// declare variables that will be used in the plotting/formatting of histograms and file names
    std::string info ("");
    std::string cut_case ("");
    std::string elSectors ("");
    std::string deltaPhiEl ("");
    std::string thetaEl ("");
    std::string elMom ("");
    std::string protSectors ("");
    std::string deltaPhiProt ("");
    std::string thetaProt ("");
    std::string protMom ("");
    std::string cuts ("");

    // Parse TList
    // First Grab Run Info
    info = std::string(Run_Info->FindObject("Run type")->GetTitle()) + 
        "_" + std::string(Run_Info->FindObject("Target")->GetTitle()) +
        "_" + double_converter(((TParameter<double>*)Run_Info->FindObject("Beam Energy"))->GetVal()) + "_";

    std::cout << "Info = " << info << std::endl;
    out.push_back(info);

    // Check if there are any cuts
    TParameter<bool> *anycuts = (TParameter<bool>*)Run_Info->FindObject("AnyCuts");
    if(anycuts->GetVal() == true) cut_case="NoCuts_";

    // Electron and Proton Sectors
    elSectors = std::string(Run_Info->FindObject("electron sectors")->GetTitle()) + "_";
    protSectors = std::string(Run_Info->FindObject("proton sectors")->GetTitle()) + "_";

    // Electron and Proton DeltaPhi Cuts
    double El_deltaPhi = ((TParameter<double>*)Run_Info->FindObject("El delta phi"))->GetVal();
    double Prot_deltaPhi = ((TParameter<double>*)Run_Info->FindObject("Prot delta phi"))->GetVal();

    if(El_deltaPhi != 9999) { deltaPhiEl = double_converter(El_deltaPhi) + "_"; }
    if(Prot_deltaPhi != 9999) { deltaPhiProt = double_converter(Prot_deltaPhi) + "_"; }

    // Electron theta cuts
    double theta_el_lb = ((TParameter<double>*)Run_Info->FindObject("El theta lb"))->GetVal();
    double theta_el_ub = ((TParameter<double>*)Run_Info->FindObject("El theta ub"))->GetVal();

    if(theta_el_lb != 0 || theta_el_ub != 180) {thetaEl = "thetaEl_" + double_converter(theta_el_lb) 
        + "_" + double_converter(theta_el_ub) + "_";}

    // Proton theta cuts
    double theta_prot_lb = ((TParameter<double>*)Run_Info->FindObject("Prot theta lb"))->GetVal();
    double theta_prot_ub = ((TParameter<double>*)Run_Info->FindObject("Prot theta ub"))->GetVal();

    if(theta_prot_lb != 0 || theta_prot_ub != 180) {thetaProt = "thetaProt_" + double_converter(theta_prot_lb) 
        + "_" + double_converter(theta_prot_ub) + "_";}

    // Electron momentum cuts
    double el_mom_lb = ((TParameter<double>*)Run_Info->FindObject("El mom lb"))->GetVal();

    if(el_mom_lb != 0) { elMom = "elMomlb_" + double_converter(el_mom_lb) + "_"; }

    // Proton momentum cuts
    double Prot_mom_lb = ((TParameter<double>*)Run_Info->FindObject("Prot mom lb"))->GetVal();
    double Prot_mom_ub = ((TParameter<double>*)Run_Info->FindObject("Prot mom ub"))->GetVal();

    if(el_mom_lb != 0) { protMom = "ProtMom_" + double_converter(Prot_mom_lb) 
        + "_" + double_converter(Prot_mom_ub); }

    // Add all of the cuts together
    if ( cut_case.compare("noCuts") == 0) cuts = cut_case;
    else cuts = elSectors + deltaPhiEl + thetaEl + elMom + protSectors + deltaPhiProt + thetaProt + protMom;

    std::cout << "Cuts = " << cuts << std::endl;
    out.push_back(cuts);
    return out;

}