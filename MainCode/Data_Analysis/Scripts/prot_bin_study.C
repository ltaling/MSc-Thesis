// This is a file to access each of the proton bin study plotting scripts in
// this directory and go through all of them in succession so the user does not
// have to run each file individually.
//
// Note that this macro does not necessarily go over the macros that are not involved 
// in the ongoing analysis for Jacob Smith's (smithja) fall 2021 SULI project.
//
// Author: Jacob Smith (smithja)
// Date: 11/18/2021

#include <iostream>
#include <vector>
#include <regex>
#include <stdio.h> // needed to sort array of bin values
#include <stdlib.h> // needed to sort array of bin values


#include "prot_bin_study_elPhi.C"
#include "prot_bin_study_elTheta.C"
#include "prot_bin_study_elMom.C"
#include "prot_bin_study_protPhi.C"
#include "prot_bin_study_protTheta.C"


// -----------------------------------------------------------------------------------------------------------------------
// This function compares two elements to use in conjunction with
// qsort(), which will sort the array of bin values that the user
// specifies during run time.
int comp (const void * elem1, const void * elem2)
{
    Double_t f = *((Double_t*)elem1);
    Double_t s = *((Double_t*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

void prot_bin_study( std::string file_name){
	// read the file name to determine if we are plotting GENIE MC or CLAS data
	auto file_type_start_pos = file_name.rfind("/"); // find the last instance of '/'
	auto file_type_end_pos = file_name.find("_", file_type_start_pos + 1);
	std::string file_type = file_name.substr( file_type_start_pos + 1, file_type_end_pos - (file_type_start_pos + 1));

        if (file_type.compare("genie") == 0 || file_type.compare("data") == 0) { // if we have GENIE MC or CLAS data
	    	// Give a little preamble so the user can choose the bins they want to use.
    	    	std::cout << std::endl << std::endl; // separate the preamble message
    	    	Int_t n_bins;

    	    	std::cout << "How many bins would you like to make? ";
    	    	std::cin >> n_bins;
	    	const Int_t N_BINS = n_bins; // set as a constant so program doesn't output warnings about using a non-constant for array declaration

    	    	Double_t bin_bounds[N_BINS+1];

    	    	std::cout << "\nPlease specify " << N_BINS+1 << " edges for your bins." << std::endl;
    	    	std::cout << "For example, entering 1, 2, and 3 will give you two bins, one from 1 to 2 and one from 2 to 3." << std::endl;
    	    	for (int i = 0; i < N_BINS+1; i++) {
                	std::cout << "Bin edge " << i << ": ";
                	std::cin >> bin_bounds[i];
    	    	}

    	    	std::cout << "\n\nSorting your bins..." << std::endl;
    	    	qsort( bin_bounds, sizeof(bin_bounds)/sizeof(*bin_bounds), sizeof(*bin_bounds), comp);

    	    	std::cout << "You entered " << N_BINS << " bins:" << std::endl;
    	    	for (int i = 0; i < N_BINS; i++) {
                	std::cout << "{" << bin_bounds[i] << ", " << bin_bounds[i+1] << "}";
                	if (i != N_BINS-1) { std::cout << ", ";}
    	    	}

    	    	std::cout << std::endl << std::endl;
    	    	std::cout << "END OF USER INPUT SECTION" << std::endl;
    	    	std::cout << "-----------------------------------------------------------" << std::endl << std::endl;
    	    	// END OF USER INPUT -------------------------------------------------------------------------------------------------------------
		
		prot_bin_study_elPhi( file_name, bin_bounds, N_BINS);
		std::cout << std::endl; // separate program output
                prot_bin_study_elTheta( file_name, bin_bounds, N_BINS);
                std::cout << std::endl; // separate program output
                prot_bin_study_elMom( file_name, bin_bounds, N_BINS);
                std::cout << std::endl; // separate program output
                prot_bin_study_protPhi( file_name, bin_bounds, N_BINS);
                std::cout << std::endl; // separate program output
                prot_bin_study_protTheta( file_name, bin_bounds, N_BINS);
        }
        else { // if we do not have either GENIE MC or CLAS data
                std::cout << "/----------------------------------------------------------------------------------------\\" << std::endl;
	        std::cout << "| It seems you are trying to process a file with an invalid file name. Please make sure  |" << std::endl;
                std::cout << "| you have either the keyword 'genie' or the keyword 'data' at the beginning of the file |" << std::endl;
                std::cout << "| name, e.g. <PATH_TO_FILE>/genie_<OMITTED>.root or <PATH_TO_FILE>/data_<OMITTED>.root.  |" << std::endl;
                std::cout << "\\----------------------------------------------------------------------------------------/" << std::endl;
                return;
        }
}

