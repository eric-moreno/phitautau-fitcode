#!/bin/bash
echo "Running limits"
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
#masslist=(150 200 250 300)
#masslist=(50)
masslist=(125)
#addargs='"--freezeParameters allConstrainedNuisances"'
#addargs='"--skipSBFit --skipBOnlyFit "'
#nopass=""
addargs='"-t -1"'
#addargs='"--freezeParameters allConstrainedNuisances"'
mode="AsymptoticLimits"
expsig=0
if [[ "$2" == "full" ]]; then
  #addargs='"--freezeParameters rgx{.*mcstat.*}"'
  #addargs='"--freezeParameters wlnuLeffSF_hadhad"'
  #addargs='" --freezeParameters jescale,jeresol,l1prefire,uescale"'
  fullstr=" --full"
  #addargs='"--freezeParameters rgx{.*CMS_resonance.*} rgx{.*norm.*}"'
  #addargs='"-t -1 --freezeParameters allConstrainedNuisances"'
  #addargs='" -t -1"'
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
#python reloadFromPkl.py --indir /uscms/home/eamoreno/nobackup/Analysis/phitautau/boostedhiggs/cards/ --outdir $1 --pkl ${modellist[@]} ${fullstr} --expsig ${expsig} --mode ${mode} --masses ${masslist[@]} --addargs "${addargs}" ${nopass}
echo "here"
OUTPUT=($1)

patch_datacards() {
  shopt -s nullglob
  for card in *.txt; do
    # Count how many lines we will convert (so we can bump kmax accordingly)
    cnt_wlnu=$(grep -Ec '^[[:space:]]*wlnuLeffSF_(hadel|hadmu|hadhad)[[:space:]]+extArg\b' "$card" || true)
    cnt_top=$(grep -Ec '^[[:space:]]*topLeffSF_(hadel|hadmu|hadhad)[[:space:]]+extArg\b' "$card" || true)
    #cnt_top=0
    repl_count=$((cnt_wlnu + cnt_top))

    # If you want to ALWAYS add +2 whenever we made ANY change, uncomment next line:
    # [[ $repl_count -gt 0 ]] && repl_count=2

    (( repl_count > 0 )) || continue
    echo "Patching $card (+$repl_count nuisance(s))"

    # Swap extArg -> param 1 1 for both families
    sed -i -E 's/^[[:space:]]*(wlnuLeffSF_(hadel|hadmu|hadhad))[[:space:]]+extArg\b.*$/\1 param 1 1/' "$card"
    sed -i -E 's/^[[:space:]]*(topLeffSF_(hadel|hadmu|hadhad))[[:space:]]+extArg\b.*$/\1 param 1 1/' "$card"

    # Bump kmax by the number of replacements, only if kmax is numeric (not '*')
    if grep -Eq '^[[:space:]]*kmax[[:space:]]+[0-9]+\b' "$card"; then
      # GNU sed: change only the first kmax line and do arithmetic with $repl_count
      sed -i -E "0,/^([[:space:]]*kmax[[:space:]]+)([0-9]+)(.*)$/s//echo \"\1\$((\2+$repl_count))\3\"/e" "$card"
    else
      echo "  (kmax is '*' or non-numeric; left unchanged.)"
    fi
  done
}

for MODEL in ${modellist[@]}; do
  echo ${MODEL}
  for MASS in ${masslist[@]}; do
    echo ${MASS}
    cd ${OUTPUT[-1]}/${MODEL}_m${MASS}/

    # Remove just the text 'passhadel2017=passhadel2017.txt', 'passhadmu2017=passhadmu2017.txt', or 'passhadhad2017=passhadhad2017.txt' without deleting the whole line
    # if [[ -f "build.sh" ]]; then
    #   sed -i 's/passhadel2017=passhadel2017.txt//g; s/passhadmu2017=passhadmu2017.txt//g; s/passhadhad2017=passhadhad2017.txt//g' build.sh
    #   sed -i 's/passhadel2016=passhadel2016.txt//g; s/passhadmu2016=passhadmu2016.txt//g; s/passhadhad2016=passhadhad2016.txt//g' build.sh
    #   sed -i 's/passhadel2016APV=passhadel2016APV.txt//g; s/passhadmu2016APV=passhadmu2016APV.txt//g; s/passhadhad2016APV=passhadhad2016APV.txt//g' build.sh
    #   sed -i 's/passhadel2018=passhadel2018.txt//g; s/passhadmu2018=passhadmu2018.txt//g; s/passhadhad2018=passhadhad2018.txt//g' build.sh
    #   echo "Removed assignments related to passhadel, passhadmu, passhadhad from build.sh in ${MODEL}_m${MASS}"
    # fi

    # if [[ -f "build_gof.sh" ]]; then
    #   sed -i 's/passhadel2017=passhadel2017.txt//g; s/passhadmu2017=passhadmu2017.txt//g; s/passhadhad2017=passhadhad2017.txt//g' build_gof.sh
    #   sed -i 's/passhadel2016=passhadel2016.txt//g; s/passhadmu2016=passhadmu2016.txt//g; s/passhadhad2016=passhadhad2016.txt//g' build_gof.sh
    #   sed -i 's/passhadel2016APV=passhadel2016APV.txt//g; s/passhadmu2016APV=passhadmu2016APV.txt//g; s/passhadhad2016APV=passhadhad2016APV.txt//g' build_gof.sh
    #   sed -i 's/passhadel2018=passhadel2018.txt//g; s/passhadmu2018=passhadmu2018.txt//g; s/passhadhad2018=passhadhad2018.txt//g' build_gof.sh
     
    #   echo "Removed assignments related to passhadel, passhadmu, passhadhad from build_gof.sh in ${MODEL}_m${MASS}"
    # fi

    # Remove just the text 'passhadel2017=passhadel2017.txt', 'passhadmu2017=passhadmu2017.txt', or 'passhadhad2017=passhadhad2017.txt' without deleting the whole line
    if [[ -f "build.sh" ]]; then
      sed -i 's/wlnuCRfailhadel2017=wlnuCRfailhadel2017.txt//g; s/wlnuCRfailhadmu2017=wlnuCRfailhadmu2017.txt//g; s/wlnuCRfailhadhad2017=wlnuCRfailhadhad2017.txt//g' build.sh
      sed -i 's/wlnuCRfailhadel2016=wlnuCRfailhadel2016.txt//g; s/wlnuCRfailhadmu2016=wlnuCRfailhadmu2016.txt//g; s/wlnuCRfailhadhad2016=wlnuCRfailhadhad2016.txt//g' build.sh
      sed -i 's/wlnuCRfailhadel2016APV=wlnuCRfailhadel2016APV.txt//g; s/wlnuCRfailhadmu2016APV=wlnuCRfailhadmu2016APV.txt//g; s/wlnuCRfailhadhad2016APV=wlnuCRfailhadhad2016APV.txt//g' build.sh
      sed -i 's/wlnuCRfailhadel2018=wlnuCRfailhadel2018.txt//g; s/wlnuCRfailhadmu2018=wlnuCRfailhadmu2018.txt//g; s/wlnuCRfailhadhad2018=wlnuCRfailhadhad2018.txt//g' build.sh
      echo "Removed assignments related to wlnuCRfailhadel, wlnuCRfailhadmu, wlnuCRfailhadhad from build.sh in ${MODEL}_m${MASS}"
    fi

    if [[ -f "build_gof.sh" ]]; then
      sed -i 's/wlnuCRfailhadel2017=wlnuCRfailhadel2017.txt//g; s/wlnuCRfailhadmu2017=wlnuCRfailhadmu2017.txt//g; s/wlnuCRfailhadhad2017=wlnuCRfailhadhad2017.txt//g' build_gof.sh
      sed -i 's/wlnuCRfailhadel2016=wlnuCRfailhadel2016.txt//g; s/wlnuCRfailhadmu2016=wlnuCRfailhadmu2016.txt//g; s/wlnuCRfailhadhad2016=wlnuCRfailhadhad2016.txt//g' build_gof.sh
      sed -i 's/wlnuCRfailhadel2016APV=wlnuCRfailhadel2016APV.txt//g; s/wlnuCRfailhadmu2016APV=wlnuCRfailhadmu2016APV.txt//g; s/wlnuCRfailhadhad2016APV=wlnuCRfailhadhad2016APV.txt//g' build_gof.sh
      sed -i 's/wlnuCRfailhadel2018=wlnuCRfailhadel2018.txt//g; s/wlnuCRfailhadmu2018=wlnuCRfailhadmu2018.txt//g; s/wlnuCRfailhadhad2018=wlnuCRfailhadhad2018.txt//g' build_gof.sh
     
      echo "Removed assignments related to wlnuCRfailhadel, wlnuCRfailhadmu, wlnuCRfailhadhad from build_gof.sh in ${MODEL}_m${MASS}"
    fi

    # Patch datacards to convert wlnuLeffSF_(hadel|hadmu|hadhad) from extArg to param
    patch_datacards

    #source build.sh >& out.log
    #source build.sh
    #source build_gof.sh >& out_gof.log
    source build_gof.sh
    if [[ "${mode}" == "Significance" ]]; then
      grep "Significance:" out.log
    elif [[ "${mode}" == "AsymptoticLimits" ]]; then
      grep " -- AsymptoticLimits ( CLs ) --" out.log -A 5
    fi
    cd -
  done
done