// This is a file to access each of the comparison plotting scripts in 
// this directory and go through them all in succession so the user does 
// not have to run each macro.
//
// Note that this macro does not necessarily go over the macros that are not involved 
// in the ongoing analysis for Jacob Smith's (smithja) fall 2021 SULI project.
//
// Author: Jacob Smith (smithja)
// Date: 11/12/2021

#include "genie_clas_el_mom_comp.C"
#include "genie_clas_el_phiAndTheta_comp.C"
#include "genie_clas_prot_mom_comp.C"
#include "genie_clas_prot_phiAndTheta_comp.C"

void genie_clas_comparison_plots( std::string mc_file_name, std::string data_file_name){
	// call each of the plotting scripts
	genie_clas_el_mom_comp( mc_file_name, data_file_name);
	genie_clas_el_phiAndTheta_comp( mc_file_name, data_file_name);
	genie_clas_prot_mom_comp( mc_file_name, data_file_name);
	genie_clas_prot_phiAndTheta_comp( mc_file_name, data_file_name);
}

