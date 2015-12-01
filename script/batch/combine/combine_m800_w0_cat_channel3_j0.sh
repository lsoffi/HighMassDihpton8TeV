#!/bin/sh
cd /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/script
eval `scramv1 runtime -sh`
cd 
mkdir workspaces
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/datacardWithAllSyst/HighMass-hgg_8TeV_m800.00_w0.10_channel3_GGH.txt ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/HighMass-hgg_8TeV_m800.00_w0.10.inputsig_GGH.root ./workspaces/
cp  /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/HighMass-hgg.inputbkg_m800.00.root ./workspaces/

combine HighMass-hgg_8TeV_m800.00_w0.10_channel3_GGH.txt -M Asymptotic -m 800 -D data_obs --run=both  -U  -S 1 -n _w0.10_channel3 --verbose 3
mv -v higgsCombine_w0.10_channel3.Asymptotic.mH800.root /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/channel3/
