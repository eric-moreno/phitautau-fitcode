#!/bin/bash

if [[ -n "$5" ]]; then
  modellist=(${4} ${5})
elif [[ -n "$4" ]]; then
  modellist=(${4})
else
  modellist=(hadhad hadmu hadel)
fi

year=$3

masslist=(125)

for mass in ${masslist[@]}; do
  if [[ -n "$2" ]]; then
    dirname=$1_$2_m${mass}
    mkdir -p ${dirname}/
    mkdir -p $1_$2/


    scp /uscms_data/d3/eamoreno/Analysis/phitautau/CMSSW_10_2_13/src/phitautau-fitcode/${1}/${2}/{full,had*}Model_m${mass}/*{Model.root,pdf} ${dirname}/
    scp /uscms_data/d3/eamoreno/Analysis/phitautau/CMSSW_10_2_13/src/phitautau-fitcode/${1}/${2}/limit*pdf ${1}_${2}/
  else
    dirname=$1_m${mass}
    mkdir -p ${dirname}/
  fi

  cp quick_distplot_phi.py ${dirname}/
  cd ${dirname}/

  for cat in ${modellist[@]}; do
    if [[ ${cat} == "hadhad" ]]; then
      prelist=(SR wlnuCR topCR)
      #prelist=(SR)
    elif [[ ${cat} == "limit" ]]; then
      break
    else
      prelist=(SR wlnuCR topCR)
      #prelist=(SR)
    fi
    for prefix in ${prelist[@]}; do
      for ptype in prefit fit_b fit_s; do
        if [[ ${cat} == "full" ]]; then
          for chan in hadhad hadel hadmu; do
            if [ ${prefix} == "SR" ]; then
              python quick_distplot_phi.py --infile fitDiagnostics.${cat}Model.root --tree shapes_${ptype} --regions fail${chan}${year} loosepass${chan}${year}${year} --label ${cat}_${chan}_${ptype}_${year} --prefix ${prefix} --channel ${chan}
              #python quick_distplot_phi.py --infile fitDiagnostics.${cat}Model.root --tree shapes_${ptype} --regions ${prefix}fail${chan}${year} ${prefix}loosepass${chan}${year}${year} --label ${cat}_${chan}_${ptype}_${year}
            else
              python quick_distplot_phi.py --infile fitDiagnostics.${cat}Model.root --tree shapes_${ptype} --regions ${prefix}fail${chan}${year} ${prefix}loosepass${chan}${year}${year} ${prefix}pass${chan}${year} --label ${cat}_${chan}_${ptype}_${year} --prefix ${prefix} --channel ${chan}
            fi
          done
        else
          if [ ${prefix} == "SR" ]; then
            echo ptype ${ptype}
            python quick_distplot_phi.py --infile fitDiagnostics.${cat}Model.root --tree shapes_${ptype} --regions fail${cat}${year} loosepass${cat}${year} --label ${cat}_${ptype}_${year} --prefix ${prefix} --channel ${cat}
            #python quick_distplot_phi.py --infile fitDiagnostics.${cat}Model.root --tree shapes_${ptype} --regions ${prefix}fail${cat}${year} ${prefix}loosepass${cat}${year} --label ${cat}_${ptype}_${year}
          else
            python quick_distplot_phi.py --infile fitDiagnostics.${cat}Model.root --tree shapes_${ptype} --regions ${prefix}fail${cat}${year} ${prefix}loosepass${cat}${year} ${prefix}pass${cat}${year} --label ${cat}_${ptype}_${year} --prefix ${prefix} --channel ${cat}
          fi
          
        fi
      done
    done
  done
  cd -
done
