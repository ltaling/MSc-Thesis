#! /usr/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -a Acceptance -p Number of Protons"
   echo -e "\t-a Would you like to turn detector acceptance on (1) or off (0)?"
   echo -e "\t-p How many protons (-1,0,1,2,3)? -1 means any number of protons"
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

# Initialising variables
if [ "$Acc" = "1" ]; then
  Acceptance="Acc"
else
  Acceptance="noAcc"
fi

declare -a Programs=("h1_pion_mom.C" "h1_pion_phi.C" "h1_pion_theta.C" "h1_pion_w.C" "h1_pion_q2.C" "h1_pion_omega.C" "h1_pion_electrontheta.C" "h1_pion_electronmomentum.C" "h1_pion_x.C")
ProgramPath="/genie/app/users/lataling/e4nu/Scripts/"
FilePath="../Sorted/C12/"
type="C12_2.261000"

#Go to the script location
cd $ProgramPath

for prog in "${Programs[@]}"; do
  root $prog'("'"$FilePath""Genie_1_""$Acceptance"'_'"$Nump"'P_'"$type"'.root",1 ,0 ,0)' -l -b -q
 # root $prog'("'"$FilePath""Genie_1_""$Acceptance"'_'"$Nump"'P_'"$type"'.root",0 ,1 ,0)' -l -b -q
  root $prog'("'"$FilePath""Genie_1_""$Acceptance"'_'"$Nump"'P_'"$type"'.root",0 ,0 ,1)' -l -b -q
done
