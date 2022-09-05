// This is a file to access each of the primary (i.e. non-comparison such as
// genie_class_prot_mom_comp.C) plotting scripts in this directory and go through them
// all in succession so the user does not have to run each macro.
//
// Note that this macro does not necessarily go over the macros that are not involved 
// in the ongoing analysis for Jacob Smith's (smithja) fall 2021 SULI project.
//
// Author: Jacob Smith (smithja)
// Date: 10/1/2021

#include "h1_el_mom.C"
#include "h1_el_mom_intSeparated.C"
#include "el_phiAndTheta_plot.C"
#include "h2_el_prot_phi.C"
#include "h2_el_prot_theta.C"
#include "h1_prot_mom.C"
#include "h1_prot_mom_intSeparated.C"
#include "h1_prot_mom_intAndFSISeparated.C"
#include "prot_phiAndTheta_plot.C"
#include "h1_el_theta_intSeparated.C"
#include "h1_prot_theta_intSeparated.C"
#include "h1_prot_theta_intAndFSISeparated.C"

void primary_plots( std::string file_name){
	// read the file name to determine if we are plotting GENIE MC or CLAS data
	auto file_type_start_pos = file_name.rfind("/"); // find the last instance of '/'
	auto file_type_end_pos = file_name.find("_", file_type_start_pos + 1);
	std::string file_type = file_name.substr( file_type_start_pos + 1, file_type_end_pos - (file_type_start_pos + 1));

        // only	plot these distributions if we are working with	GENIE MC data–the
        // CLAS	experiment cannot separate things into interaction components or
        // FSI vs. no FSI contributions, so these plotting scripts will always 
        // generate empty histograms
        if (file_type.compare("genie") == 0) { // if we have GENIE MC
                h1_el_mom_intSeparated( file_name); // this plotting script uses a plot from the inclusive section of the code
                h1_el_theta_intSeparated( file_name);

                h1_prot_mom_intSeparated( file_name);
                h1_prot_mom_intAndFSISeparated( file_name);
                h1_prot_theta_intSeparated( file_name);
                h1_prot_theta_intAndFSISeparated( file_name);
        }
	else if (file_type.compare("data") == 0) { // if we have CLAS data
		h1_el_mom( file_name); //this plotting script uses a plot from the inclusive section of the code
	}
        else { // if we do not have either GENIE MC or CLAS data
                std::cout << "/----------------------------------------------------------------------------------------\\" << std::endl;
	        std::cout << "| It seems you are trying to process a file with an invalid file name. Please make sure  |" << std::endl;
                std::cout << "| you have either the keyword 'genie' or the keyword 'data' at the beginning of the file |" << std::endl;
                std::cout << "| name, e.g. <PATH_TO_FILE>/genie_<OMITTED>.root or <PATH_TO_FILE>/data_<OMITTED>.root.  |" << std::endl;
                std::cout << "\\----------------------------------------------------------------------------------------/" << std::endl;
                return;
        }

	// run these plotting scripts regardless of if we have GENIE MC or CLAS data
	el_phiAndTheta_plot( file_name);	
	prot_phiAndTheta_plot( file_name);
        h1_prot_mom( file_name);


// -----------------------------------------------------------
// These files have yet to be worked into run_all_scripts.C.
//	h2_el_prot_phi( file_name);
//	h2_el_prot_theta( file_name);
}

