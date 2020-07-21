# PMT Darkrate #
In this repository there is the neceissary to select and study the dark-rates pulses.
## makeHistPulses ## 
For now we focus on the darkrate in cold. The script makeHistPulses.cc selects the pulses found on each event and groups them on the appropriate PMT. It applies also some selection. a bash script, `processdarkrates.sh`  allows to create to outputs with and witouth cuts. Just do ``` source processdarkrates.sh```
## studty the Threshold ## 
the macro `studythreshold.cc` produce a graph showing the rate as function of a selection threshold in mV on the pulse amplitude
