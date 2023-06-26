# Datacards

## Make datacards

Use directory $IDIR with latest histograms:
```
2016: TBD
2017: /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/condor/Nov30_2017_UL/hists_sum_
2018: TBD
```

Use the following commmands:
```
# for 2017
python make_cards.py --hist /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/condor/Nov30_2017_UL/hists_sum_ --year 2017 --cat hadel --cat hadmu --cat hadhad --tag Jun19_2017
```