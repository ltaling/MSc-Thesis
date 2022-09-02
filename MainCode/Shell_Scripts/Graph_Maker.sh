#! /usr/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -a Acceptance -p Number of Protons"
   echo -e "\t-a Would you like to turn detector acceptance on (1) or off (0)?"
   echo -e "\t-p How many protons (0,1,2,3,n)?"
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

declare -a Programs=("h1_pion_mom.C" "h1_pion_phi.C" "h1_pion_theta.C" "h2_pion_mom_phi.C" "h2_pion_mom_theta.C")
ProgramPath="/genie/app/users/lataling/e4nu/Scripts/"
FilePathHe="../Sorted/4He/"
FilePathFe="../Sorted/56Fe/"
FilePathC="../Sorted/C12/"
Words="Genie_1_"
declare -a nucen=("4He_2.261000" "4He_4.461000" "56Fe_2.261000" "56Fe_4.461000" "C12_1.161000" "C12_4.461000")

#Go to the script location
cd $ProgramPath

for prog in "${Programs[@]}"; do
  for types in "${nucen[@]}"; do
    for pimin in 0 1; do
      for pizero in 0 1; do
        for piplus in 0 1; do
          if [ $pimin = 0 ] && [ $pizero = 0 ] && [ $piplus = 0 ]; then
            continue
          elif [[ "$types" == *"4He"* ]]; then
              root $prog'("'"$FilePathHe""$Words""$Acceptance"'_'"$Nump"'P_'"$types"'.root",'"$pimin"' ,'"$pizero"' ,'"$piplus"')' -l -b -q
	  elif [[ "$types" == *"56Fe"* ]]; then
              root $prog'("'"$FilePathFe""$Words""$Acceptance"'_'"$Nump"'P_'"$types"'.root",'"$pimin"' ,'"$pizero"' ,'"$piplus"')' -l -b -q
	  elif [[ "$types" == *"C12"* ]]; then
              root $prog'("'"$FilePathC""$Words""$Acceptance"'_'"$Nump"'P_'"$types"'.root",'"$pimin"' ,'"$pizero"' ,'"$piplus"')' -l -b -q
          fi
        done
      done
    done
  done
done
