<<<<<<< HEAD
# Fitting code for PhiTauTau

## Setup

Follow the setup instructions for combine with CMSSW_10_2_13.
Install `rhalphalib` locally (there may be some modifications necessary, things don't always play nice by default).
Clone this repo into `CMSSW_10_2_13/src/`.

### To install rhalphalib:

1. Install rhalphalib
   git checkout -b phitautau origin/phitautau
2. In CMSSW_10_2_13
   scram b
   pip install flake8 --user
   pip install --user https://github.com/drankincms/rhalphalib/archive/master.zip

## Card Making
First, make your cards!
```
python3 make_cards.py --hist /uscms_data/d3/eamoreno/Analysis/phitautau/boostedhiggs/condor/Jun28_2017_UL/hists_sum_  --year 2017 --cat hadel --cat hadmu --cat hadhad --tag Feb7_2017_UL/1
```

Dylan's samples (messed up massreg hadhad) - `/uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/condor/Nov30_2017_UL/hists_sum_`

Eric's samples (new massreg hadhad retraining) - `/uscms_data/d3/eamoreno/Analysis/phitautau/boostedhiggs/condor/Jun28_2017_UL/hists_sum_`

## Running

Next run limits!

`runLimits.sh` is the top script, which will call the other scripts in this directory as needed. An example invocation is below where you call ALL the limits run:
```
./runLimits.sh Feb7_2017_UL/1 full hadel hadmu hadhad
```

If you want to try a specific channel (e.g. hadel) and only one template(e.g. m200) you can run 
```
./runLimits.sh Feb7_2017_UL/1 full200 hadel 
```

you might need to hit `runLimits.sh` for specifics - for example there is `-t -1 --toysFrequentist` enabled right now

## Local scripts

While the fitting is expected to run remotely, the plotting etc can easily be run locally. But they don''t need to be. This makes pretty mass plots. The scripts to copy down the necessary files are in `local_scripts/`. A standard invocation is below:
```
cd local_scripts/
./plotPhi.sh Feb7_2017_UL 1 2017
```

## Make Brazil-band Plots (Run Phi Limits)

To run this, you should have multiple-mass limits run. If you just run one limit the plot will be of zero width and fail. Here is an example: 

```
./makePhiLimits.sh Feb7_2018_UL/1
```
## Old scripts

```
python test/makeCardsPhi.py --hist condor/Nov30_2017_UL/hists_sum_ --year 2017 --lumi 41.5 --tag Dec12_2017 --label 34 --hPtBinsLep None --hPtCut -1 --hPtBinsHad None --shapeRegionsHad fail loosepass loosepass --metCutLep 75. --lowMetCutHad 75. --unblind --singleBinLepCR --singleBinHadFail --highmassone --unifiedBkgEff

python test/makeCardsPhi.py --hist /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/condor/Nov30_2017_UL/hists_sum_ --year 2017 --lumi 41.5 --tag Dec12_2017 --label 34 --hPtBinsLep None --hPtCut -1 --hPtBinsHad None --shapeRegionsHad fail loosepass loosepass --metCutLep 75. --lowMetCutHad 75. --unblind --singleBinLepCR --singleBinHadFail --highmassone --unifiedBkgEff
```


## New Information



=======
### Installation

```
gh repo clone andrzejnovak/combine_postfits
cd combine_postfits
pip install -e .
```

### Run

Example script to be modified as needed can be found in `make_plots.py. Run as:

```
 python make_plots.py -i hadelModel_m125/fitDiagnosticsTest.root -v --MC --style style_jeff.yml
```
>>>>>>> 3231241b96a7b8474965d04f8941b0e5564ed86f
