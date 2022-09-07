### These scripts read in the file output from the program Genie_analysis and plots the histograms

## How to Use

Once a file has been created by running the program Genie_analysis, it can be read by the script and a graph plotted. There are two types of scripts; the ones which were there from the beginning only require an input file whereas the ones made by me have additional parameters. To see which input parameters are required, just read the source code and the function definition will tell you. The general template to run the script is:

`root 'ScriptName.C("root file name and path", params)'`

A specific example, for one of my pion scripts (all pion scripts are written by me), would read `root 'h1_pion_mom.C("../Sorted/C12/Genie_1_Acc_-1P_C12_2.261000.root",1,0,0)'`. This plots a 1D histogram of $\pi^-$ momentum for 2.261 GeV scattering off C-12 with any number of protons and with detector acceptance on. 

NB: It is possible to combine pion types so writing 1,0,1 in the parameter region would plot it for both negative and positive pions. The graph title adjusts automatically.

The produced graphs can be seen in the outputs folder (if you're unsure where it is, just look at the source code and copy the file path).
