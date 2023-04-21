Master's project 2023

To run code:

1. scram b
2. (time) cmsRun runAnalysis.py

This will produce unedited plots that have to be accessed through ROOT:

1. root -l stau_M432_analysis.root
2. .ls
3. [plot name]->Draw()

To have edited plots:

1. root -l stau_M432_analysis.root
2. .x genrecofigures.C
