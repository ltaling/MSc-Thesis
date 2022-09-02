#! /usr/bin/bash

# Initialise the shell script with variables
helpFunction()
{
   echo ""
   echo "Usage: $0 -a Acceptance -p Number of Protons"
   echo -e "\t-a Would you like to turn detector acceptance on (1) or off (0)?"
   echo -e "\t-p How many protons (0,1,2,3)?"
   exit 1 # Exit script after printing help
}

while getopts "a:p:" opt
do
   case "$opt" in
      a ) Acc="$OPTARG" ;;
      p ) Nump="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$Acc" ] || [ -z "$Nump" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Define all the file/program/output names and paths here
ProgramName="genie_analysis"
ProgramPath="/genie/app/users/lataling/e4nu/Analyze_data"
declare -a TargetEnergy=("4He 2261" "4He 4461" "56Fe 2261" "56Fe 4461" "C12 1161" "C12 4461") #Leave out "C12 2261" as it doesn't work atm

#Go to Program Directory
cd $ProgramPath

#Run the Program with all combinations
for nuc in "${TargetEnergy[@]}"; do
  ./$ProgramName $nuc 1 1 0 0 0 180 0 0 0 0 180 0 9999 $Nump $Acc
done
