# PMT GAIN #
This repository contains all the necessary to measure the PMT gain and save the results to a .csv file. I report here the relevant command lines to have a quick start with this task.
## Make the histogram ##
For the direct light case:
```
makeHistDirect -i calofilesDirect.txt -o ../data/calo/histogramDirectlight.cc -t calo/chargetree
```
For the indirect light case (WEST side PMTs only):
```
makeHistIndirect -i calofilesIndirect.txt -o ../data/calo/histogramIndirectlight.cc -t calo/indirectlight
```
## Fit direct light ##
Fit the direct light histograms (NB if you want to fit all 360 pmts it takes ~2h). If not, just change the `startch` and `endch` variables inside `fitDirectLight.cc` and recompile. Command line to fit histograms
```
fitDirectLight -i ../data/calo/histogramDirectLight.root -o fitDirectLight.root -d fitDirectLight.csv
```
## Fit indirect light ##
to be do
