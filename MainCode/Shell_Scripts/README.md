# ReadMe

## Description
These scripts can be run within a directory and they'll run all the code for you so you don't need to keep pressing "enter" for 30 years. However, the directories in these folders are specific to the one on my Fermilab account so they should be amended before using. To find out which variables to include, just type `./filename.sh -h` in the command line. 

### File Contents

- File_Maker.sh - runs the Genie analysis with the input parameters for all targets and energies
- Graph_Maker_2.sh - draws all pion related histograms for the target/energy specified in the code (to change the target/energy, you'll need to edit the source code) 

## Once downloaded
If you try to run it, you'll most likely receive an error. To make it executable, do the following steps:

- Type `chmod u+r+x FileName.sh` and press enter
- Then, type `sed -i -e 's/\r$//' ./FileName.sh` and press enter

This needs to be done for each file and it should be displayed with green text when you execute `ls`.  




