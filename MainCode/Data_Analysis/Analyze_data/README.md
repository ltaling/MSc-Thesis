## How To Use The Genie Analysis Program

This folder contains the program that sorts out the GENIE file into the histograms and also counts events. Once the folder is downloaded, you need to edit the filepaths to match your directory. Assuming the folder structure is kept the same, only line 295 needs to be edited which is the location of the *output* file.

Note: It has been heavily edited with `/* ... */` commenting out large parts of code to improve the running time. As a result, many features are missing and may need to be "reactivated" depending on what the project requires. All edits made by me are labelled with `lataling` before it so just do `ctrl+f` and search for "lataling" to find all edits made. 

After making any edits and when running for the first time, run `make clean; make` in the command line. There will be a lot of warnings but there should be no errors. You can remove the warnings but there are a lot and I don't think it'll improve the performance much.

Once the program is compiled, run the code via:

./genie_analysis [target] [beam energy in MeV] [0-3 depending on data (0) or simulation (1-3) to analyze] [number of rotations] [1/0 to use/don't use electron sector cuts] [1/0 to use/don't use electron phi cut] [electron theta min cut] [electron theta max cut] [electron momentum min cut] [1/0 to use/don't use proton sector cuts] [1/0 to use/don't use proton phi cut] [proton theta min cut] [proton theta max cut] [proton momentum min cut] [proton momentum max cut] [number of final state protons (0-3) with -1 meaning any number of protons] [1/0 to turn detector acceptance on/off] 

For example, `./genie_analysis C12 2261 1 1 0 0 0 180 0 0 0 0 180 0 9999 -1 1`. This analyzes the 2.261 GeV C12 file with no cuts, any number of protons in the final state and with detector acceptance on.

The momentum cut can be turned on/off by setting `bool Applymomthresh = true/false;` on line 98. 
