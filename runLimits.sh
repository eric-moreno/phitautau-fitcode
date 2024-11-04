#!/bin/bash

if [[ -n "$5" ]]; then
  modellist=(${3}Model ${4}Model ${5}Model)
elif [[ -n "$4" ]]; then
  modellist=(${3}Model ${4}Model)
elif [[ -n "$3" ]]; then
  modellist=(${3}Model)
else
  modellist=(fullModel hadelModel hadhadModel hadmuModel)
fi

fullstr=""
#masslist=(10 20 30 40 50 75 100 125 150 200 250 300)
#masslist=(30 40 50 75 100 125 150 200 250 300)
#masslist=(50 125 200)
#masslist=(30 200)
masslist=(50)
#addargs='"--toysFreq"'
#addargs='""'
#nopass=""
mode="AsymptoticLimits"
expsig=0
if [[ "$2" == "full" ]]; then
  #addargs='" --freezeParameters rgx{.*multijet_mcstat.*}"'
  fullstr=" --full"
  #addargs='"-t -1 --freezeParameters allConstrainedNuisances"'
  #addargs='"-t -1 --toysFreq"'
elif [[ "$2" == "full200" ]]; then
  fullstr=" --full"
  masslist=(200)
  #addargs='" --freezeParameters jescale,jeresol,l1prefire,uescale"'
  #addargs='" --freezeParameters allConstrainedNuisances"'
  #addargs='""'
elif [[ "$2" == "fullnopass" ]]; then
  fullstr=" --full"
  masslist=(125)
  #addargs='" --freezeParameters r_dy_hadmu,dy_eff_hadmu,r_dy_hadel,dy_eff_hadel,r_dy_hadhad,dy_eff_hadhad"'
  addargs='""'
  nopass="--nopass"
elif [[ "$2" == "fullhtt" ]]; then
  fullstr=" --full"
  masslist=(125)
  expsig=1
  mode="Significance"
elif [[ "$2" == "htt" ]]; then
  masslist=(125)
  expsig=1
  mode="Significance"
elif [[ "$2" == "fullhttnopass" ]]; then
  masslist=(125)
  expsig=1
  mode="Significance"
  addargs='" --freezeParameters r_dy_hadmu,dy_eff_hadmu,r_dy_hadel,dy_eff_hadel,r_dy_hadhad,dy_eff_hadhad"'
  nopass="--nopass"
fi

#python reloadFromPkl.py --indir /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/cards/ --outdir $1 --pkl ${modellist[@]} ${fullstr} --addargs '" --freezeParameters wlnueffSF_hadhad,wlnueffSF_hadel,wlnueffSF_hadmu,topeffSF_hadhad,topeffSF_hadel,topeffSF_hadmu,wlnuLeffSF_hadhad,wlnuLeffSF_hadel,wlnuLeffSF_hadmu,topLeffSF_hadhad,topLeffSF_hadel,topLeffSF_hadmu,toppt,ztt_eff_hadhad,ztt_eff_hadel,ztt_eff_hadmu,r_ztt_hadhad,r_ztt_hadel,r_ztt_hadmu,wlnunormSF_hadhad,wlnunormSF_hadel,wlnunormSF_hadmu,topnormSF_hadhad,topnormSF_hadel,topnormSF_hadmu"'
#python reloadFromPkl.py --indir /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/cards/ --outdir $1 --pkl ${modellist[@]} ${fullstr} #--addargs '" --freezeParameters wlnueffSF_hadhad,wlnueffSF_hadel,wlnueffSF_hadmu,topeffSF_hadhad,topeffSF_hadel,topeffSF_hadmu,wlnuLeffSF_hadhad,wlnuLeffSF_hadel,wlnuLeffSF_hadmu,topLeffSF_hadhad,topLeffSF_hadel,topLeffSF_hadmu"'
#python reloadFromPkl.py --indir /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/cards/ --outdir $1 --pkl ${modellist[@]} ${fullstr} --expsig 0 --mode AsymptoticLimits --masses ${masslist[@]} --addargs '" --freezeParameters r_dy_hadmu,dy_eff_hadmu,r_dy_hadel,dy_eff_hadel,r_dy_hadhad,dy_eff_hadhad"' --nopass
#python reloadFromPkl.py --indir /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/cards/ --outdir $1 --pkl ${modellist[@]} ${fullstr} --expsig 0 --mode AsymptoticLimits --masses ${masslist[@]} --addargs '"-t -1"'
python reloadFromPkl.py --indir /uscms_data/d3/eamoreno/Analysis/phitautau/CMSSW_10_2_13/src/phitautau-fitcode/cards/ --outdir $1 --pkl ${modellist[@]} ${fullstr} --expsig ${expsig} --mode ${mode} --masses ${masslist[@]} --addargs "${addargs}" ${nopass} 
#python reloadFromPkl_old.py --indir /uscms/home/eamoreno/nobackup/Analysis/phitautau/boostedhiggs/cards/ --outdir $1 --pkl ${modellist[@]} ${fullstr} --expsig ${expsig} --mode ${mode} --masses ${masslist[@]} --addargs "${addargs}" ${nopass}
echo "here"
OUTPUT=($1)

for MODEL in ${modellist[@]}; do
  echo ${MODEL}
  for MASS in ${masslist[@]}; do
    echo ${MASS}
    cd ${OUTPUT[-1]}/${MODEL}_m${MASS}/

    # Remove just the text 'passhadel2017=passhadel2017.txt', 'passhadmu2017=passhadmu2017.txt', or 'passhadhad2017=passhadhad2017.txt' without deleting the whole line
    if [[ -f "build.sh" ]]; then
      sed -i 's/passhadel2017=passhadel2017.txt//g; s/passhadmu2017=passhadmu2017.txt//g; s/passhadhad2017=passhadhad2017.txt//g' build.sh
      echo "Removed assignments related to passhadel, passhadmu, passhadhad from build.sh in ${MODEL}_m${MASS}"
    fi

    if [[ -f "build_gof.sh" ]]; then
      sed -i 's/passhadel2017=passhadel2017.txt//g; s/passhadmu2017=passhadmu2017.txt//g; s/passhadhad2017=passhadhad2017.txt//g' build_gof.sh
      echo "Removed assignments related to passhadel, passhadmu, passhadhad from build_gof.sh in ${MODEL}_m${MASS}"
    fi
    
    source build.sh >& out.log
    source build_gof.sh >& out_gof.log
    if [[ "${mode}" == "Significance" ]]; then
      grep "Significance:" out.log
    elif [[ "${mode}" == "AsymptoticLimits" ]]; then
      grep " -- AsymptoticLimits ( CLs ) --" out.log -A 5
    fi
    cd -
  done
done