// File to plot histograms of the proton theta for a simulation. The histograms generated
// with this macro are separated into the different interaction components: quasi-elastic (QE),
// mesonic current (MEC), resonance (RES), and deep inelastic scattering (DIS). This functionality
// is indicated by "intSeparated" being in the macro's file name.
//
// Author: Jacob Smith (smithja)
// Creation Date: 10/19/2021

// These include statements were here when I received Graham's code. Delete/Modify them at
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

void h1_el_theta_intSeparated( std::string file_name){
    // read the file name to determine if we are plotting GENIE MC or CLAS data
    auto file_type_start_pos = file_name.rfind("/"); // find the last instance of '/'
    auto file_type_end_pos = file_name.find("_", file_type_start_pos + 1);
    std::string file_type = file_name.substr( file_type_start_pos + 1, file_type_end_pos - (file_type_start_pos + 1));

    if (file_type.compare("genie") != 0){ // if we do not have GENIE MC
            std::cout << "/----------------------------------------------------------------------------------------\\" << std::endl;
            std::cout << "| It seems you are trying to process a file with an invalid file name. Please make sure  |" << std::endl;
            std::cout << "| you have the keyword 'genie' at the beginning of the file name since interaction and   |" << std::endl;
            std::cout << "| FSI separation can only be done on GENIE MC files.                                     |" << std::endl;
            std::cout << "\\----------------------------------------------------------------------------------------/" << std::endl;
            return;
    }

    // declare variables that will be used in the plotting/formatting of histograms and file names
    std::string target ("invalid");
    std::string target_plot ("invalid");
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
    else { std::cout << "You have entered an invalid/unknown target. Consider adding code to accomodate this target in h1_el_theta_intSeparated.C." << std::endl; }
    //std::cout << "target: " << target << "; target_plot: " << target_plot << std::endl;

    auto beam_kwd_pos = file_name.find("beamEnergy");
    if (beam_kwd_pos == std::string::npos) { std::cout << "'beamEnergy' keyword not found in file_name. Correct file_name!" << std::endl; }
    auto beam_end_pos = file_name.find("_", beam_kwd_pos);
    beam_energy = file_name.substr( beam_kwd_pos + 10, beam_end_pos - (beam_kwd_pos + 10));
    if (beam_energy.compare("1161") == 0) { beam_energy_plot = "1.161 GeV"; }
    else if (beam_energy.compare("2261") == 0) { beam_energy_plot = "2.261 GeV"; }
    else if (beam_energy.compare("4461") == 0) { beam_energy_plot = "4.461 GeV"; }
    else { std::cout << "You have entered an invalid/unknown beam energy. Consider adding code to accomodate this beam energy in h1_el_theta_intSeparated.C" << std::endl; }
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
        protSectors_plot = ""; // protSectors is valid so first do away with 'invalid'
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

    // Determine what case we are in (i.e. no cuts, inclusive w/o electron momentum cut,
    // inclusive, exclusive, or exclusive with protMom cut) based on the null and non-null 
    // values of the above pointers.
    std::string cut_case ("invalid");
    if (file_name.find("noCuts") != std::string::npos) { cut_case = "noCuts"; }
    if ((elSectors_kwd_pos != std::string::npos) && (deltaPhiElEq_kwd_pos != std::string::npos) && (thetaEl_kwd_pos != std::string::npos)) { cut_case = "incl_noElMomCut"; }
    if ((cut_case.compare("incl_noElMomCut") == 0) && (elMom_kwd_pos != std::string::npos)) { cut_case = "incl"; }
    if ((cut_case.compare("incl") == 0) && (protSectors_kwd_pos != std::string::npos) && (deltaPhiProtEq_kwd_pos != std::string::npos) && (thetaProt_kwd_pos != std::string::npos)) { cut_case = "excl"; }
    if ((cut_case.compare("excl") == 0) && (protMom_kwd_pos != std::string::npos)) { cut_case = "excl_withProtMomCut"; }

    // histogram initialization
    TH1F* h1_el_theta[24];

    h1_el_theta[0] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector5_el_theta"));
    h1_el_theta[1] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector5_el_theta"));
    h1_el_theta[2] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector5_el_theta"));
    h1_el_theta[3] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector5_el_theta"));

    h1_el_theta[4] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector4_el_theta"));
    h1_el_theta[5] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector4_el_theta"));
    h1_el_theta[6] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector4_el_theta"));
    h1_el_theta[7] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector4_el_theta"));

    h1_el_theta[8] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector3_el_theta"));
    h1_el_theta[9] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector3_el_theta"));
    h1_el_theta[10] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector3_el_theta"));
    h1_el_theta[11] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector3_el_theta"));

    h1_el_theta[12] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector2_el_theta"));
    h1_el_theta[13] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector2_el_theta"));
    h1_el_theta[14] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector2_el_theta"));
    h1_el_theta[15] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector2_el_theta"));

    h1_el_theta[16] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector1_el_theta"));
    h1_el_theta[17] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector1_el_theta"));
    h1_el_theta[18] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector1_el_theta"));
    h1_el_theta[19] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector1_el_theta"));

    h1_el_theta[20] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq1_NoQ4Weight_InSector0_el_theta"));
    h1_el_theta[21] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq2_NoQ4Weight_InSector0_el_theta"));
    h1_el_theta[22] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq3_NoQ4Weight_InSector0_el_theta"));
    h1_el_theta[23] = (TH1F*)input1->Get(TString::Format("h1_InteractionEq4_NoQ4Weight_InSector0_el_theta"));

    // compile all the histograms into the four interaction types (QE, MEC, RES, and DIS for quasi-elastic,
    // mesonic current, resonance, and deep inelastic scattering, respectively.
    for( int i = 1; i < 6; i++) {
        h1_el_theta[20]->Add(h1_el_theta[20-(4*i)]);
	h1_el_theta[21]->Add(h1_el_theta[21-(4*i)]); 
	h1_el_theta[22]->Add(h1_el_theta[22-(4*i)]); 
	h1_el_theta[23]->Add(h1_el_theta[23-(4*i)]); 
    }

    // create a canvas on which to draw the histograms
    TCanvas* c;
    c = new TCanvas(TString::Format("c"), TString::Format("c"), 800, 600);

    // initialize things such as color options that will be used in plotting the histograms below
    const int color_options[5] = {kBlue, kAzure+10, kGreen-3, kViolet, kBlue+2};
    const char * draw_options[5] = {"HIST","HIST SAME","HIST SAME", "HIST SAME","HIST SAME"};
    const char * draw_options2[5] = {"E SAME", "E SAME","E SAME", "E SAME","E SAME"};

    double hist_max = 0.;
    // format the histograms and then make them
    for( int i = 20; i < 24; i++){
	// formatting statements
        h1_el_theta[i]->Rebin( 4);
        h1_el_theta[i]->SetStats( 0); // get rid of the stats box that usually appears at the top right of plots
        h1_el_theta[i]->SetLineColor( color_options[i-20]);
        h1_el_theta[i]->GetXaxis()->SetTitle("#theta_{electron} [deg]");
        h1_el_theta[i]->GetYaxis()->SetTitle("Events");
        if (thetaProt_kwd_pos != std::string::npos) { h1_el_theta[i]->GetXaxis()->SetRangeUser( std::stod( thetaProt_lb), std::stod( thetaProt_ub)); }
	else { h1_el_theta[i]->GetXaxis()->SetRangeUser( 0, 180); }
        h1_el_theta[i]->GetXaxis()->CenterTitle( true);
        h1_el_theta[i]->GetYaxis()->CenterTitle( true);
//        h1_el_theta[i]->GetXaxis()->SetTitleSize( 0.06);
//        h1_el_theta[i]->GetYaxis()->SetTitleSize( 0.06);

        if (cut_case.compare("noCuts") == 0) { h1_el_theta[i]->SetTitle(TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam (No Cuts)"); }
	else if ( cut_case.compare("incl_noElMomCut") == 0) { h1_el_theta[i]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub)+"#circ}"); }
        else if ( cut_case.compare("incl") == 0) { h1_el_theta[i]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+"}"); }
        else if ( cut_case.compare("excl") == 0) { h1_el_theta[i]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt)+", "+TString(thetaProt_lb)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub)+"#circ}"); }
	else if ( cut_case.compare("excl_withProtMomCut") == 0) { h1_el_theta[i]->SetTitle("#splitline{"+TString(target_plot)+" with "+TString(beam_energy_plot)+" Beam}{#Delta#phi_{e^{-}} = "+TString(deltaPhiEl)+"#circ, e^{-} Sectors "+TString(elSectors_plot)+", "+TString(thetaEl_lb)+"#circ < #theta_{e^{-}} < "+TString(thetaEl_ub)+"#circ, p_{e^{-}} > "+TString(elMom_lb_plot)+", p^{+} Sectors "+TString(protSectors_plot)+", #Delta#phi_{p^{+}} = "+TString(deltaPhiProt)+", "+TString(thetaProt_lb)+"#circ < #theta_{p^{+}} < "+TString(thetaProt_ub)+"#circ, p_{p^{+}} > "+TString(protMom_lb)+"}"); }
    
        // take the maximum of the four interactions to use later in plotting
        if (h1_el_theta[i]->GetMaximum() > hist_max) { hist_max = h1_el_theta[i]->GetMaximum(); }
    }

    // draw the histograms
    for( int i = 20; i < 24; i++){
	h1_el_theta[i]->GetYaxis()->SetRangeUser( 0, hist_max*1.1); // smithja: hist_max mult. by 1.1 to give the graph some head room
        h1_el_theta[i]->Draw( draw_options[i-20]);
        h1_el_theta[i]->Draw( draw_options2[i-20]);
    }

    // crop the margins of the canvas
    c->SetLeftMargin( 0.14);
    c->SetBottomMargin( 0.14);
    if ( cut_case.compare("noCuts") != 0) { c->SetTopMargin( 0.16); }
    c->Update();

    // create and draw a legend for our histograms
    TLegend *legend = new TLegend( 0.7, 0.655, 1.2, 0.8);
    legend->AddEntry( h1_el_theta[20],"QE","l");
    legend->AddEntry( h1_el_theta[21],"MEC","l");
    legend->AddEntry( h1_el_theta[22],"RES","l");
    legend->AddEntry( h1_el_theta[23],"DIS","l");
    legend->SetBorderSize( 0);
    legend->SetFillStyle( 0);
    legend->Draw();

    // save this histogram as a PDF file
    std::string save_file_path ("invalid");
    if ( cut_case.compare("noCuts") == 0) { save_file_path = "../output/h1_el_theta_intSeparated/"+TString(file_type)+"_h1_el_theta_intSeparated_target"+target+"_beamEnergy"+beam_energy+"_MottXSecEq1_noCuts.pdf"; }
    else if ( cut_case.compare("incl_noElMomCut") == 0) { save_file_path = "../output/h1_el_theta_intSeparated/"+TString(file_type)+"_h1_el_theta_intSeparated_target"+TString(target)+"_beamEnergy"+TString(beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+".pdf"; }
    else if ( cut_case.compare("incl") == 0) { save_file_path = "../output/h1_el_theta_intSeparated/"+TString(file_type)+"_h1_el_theta_intSeparated_target"+TString(target)+"_beamEnergy"+TString(beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+".pdf"; }
    else if ( cut_case.compare("excl") == 0) { save_file_path = "../output/h1_el_theta_intSeparated/"+TString(file_type)+"_h1_el_theta_intSeparated_target"+TString(target)+"_beamEnergy"+TString(beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+"_protSectors"+TString(protSectors)+"_deltaPhiProtEq"+TString(deltaPhiProt)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+".pdf"; }
    else if ( cut_case.compare("excl_withProtMomCut") == 0) { save_file_path = "../output/h1_el_theta_intSeparated/"+TString(file_type)+"_h1_el_theta_intSeparated_target"+TString(target)+"_beamEnergy"+TString(beam_energy)+"_MottXSecEq1_elSectors"+TString(elSectors)+"_deltaPhiEq"+TString(deltaPhiEl)+"_thetaEl"+TString(thetaEl_lb)+"to"+TString(thetaEl_ub)+"_elMomLBEq"+TString(elMom_lb)+"_protSectors"+TString(protSectors)+"_deltaPhiProtEq"+TString(deltaPhiProt)+"_thetaProt"+TString(thetaProt_lb)+"to"+TString(thetaProt_ub)+"_protMomLBEq"+TString(protMom_lb)+".pdf"; }
    c->SaveAs( TString::Format("%s", save_file_path.c_str()));

    // smithja: this section of code integrates the distributions and prints
    //          it to the screen
    int num_QE_events = h1_el_theta[20]->Integral();
    int num_MEC_events = h1_el_theta[21]->Integral();
    int num_RES_events = h1_el_theta[22]->Integral();
    int num_DIS_events = h1_el_theta[23]->Integral();
    int num_events = num_QE_events + num_MEC_events + num_RES_events + num_DIS_events;

    cout << "// --- PROGRAM OUTPUT ------------------------------------------------------" << endl;
    cout << "There are " << num_QE_events << " many QE events in this histogram." << endl;
    cout << "There are " << num_MEC_events << " many MEC events in this histogram." << endl;
    cout << "There are " << num_RES_events << " many RES events in this histogram." << endl;
    cout << "There are " << num_DIS_events << " many DIS events in this histogram." << endl;
    cout << "There are " << num_events << " many total events in the h1_el_theta_intSeparated histogram." << endl;
    cout << "// -------------------------------------------------------------------------" << endl << endl; // separate program output
}
