#!/bin/bash

ORIGPWD=${PWD}
cd $1
MASSRANGE=(30 40 50 75 100 125 150 200)
#MASSRANGE=(30 50 100 125 150 200)
#MASSRANGE=(10 20 30 40 50 75 100 125 150 200 250 300)
eval "combineTool.py -M CollectLimits hadelModel_m{`echo ${MASSRANGE[@]} | tr ' ' ,`}/*AsymptoticLimits* -o limits_hadel.json"
eval "combineTool.py -M CollectLimits hadmuModel_m{`echo ${MASSRANGE[@]} | tr ' ' ,`}/*AsymptoticLimits* -o limits_hadmu.json"
eval "combineTool.py -M CollectLimits hadhadModel_m{`echo ${MASSRANGE[@]} | tr ' ' ,`}/*AsymptoticLimits* -o limits_hadhad.json"
eval "combineTool.py -M CollectLimits fullModel_m{`echo ${MASSRANGE[@]} | tr ' ' ,`}/*AsymptoticLimits* -o limits_comb.json"
python ${ORIGPWD}/plotLimit.py limits_hadel.json --show exp -o limit_hadel --scenario-label "e\\tau_{h}"
python ${ORIGPWD}/plotLimit.py limits_hadmu.json --show exp -o limit_hadmu --scenario-label "\\mu\\tau_{h}"
python ${ORIGPWD}/plotLimit.py limits_hadhad.json --show exp -o limit_hadhad --scenario-label "\\tau_{h}\\tau_{h}"
python ${ORIGPWD}/plotLimit.py limits_comb.json --show exp -o limit_comb 
python ${ORIGPWD}/plotLimit.py limits_hadel.json --show exp -o limit_hadel_logy --logy --scenario-label "e\\tau_{h}"
python ${ORIGPWD}/plotLimit.py limits_hadmu.json --show exp -o limit_hadmu_logy --logy --scenario-label "\\mu\\tau_{h}"
python ${ORIGPWD}/plotLimit.py limits_hadhad.json --show exp -o limit_hadhad_logy --logy --scenario-label "\\tau_{h}\\tau_{h}"
python ${ORIGPWD}/plotLimit.py limits_comb.json --show exp -o limit_comb_logy --logy
cd -
