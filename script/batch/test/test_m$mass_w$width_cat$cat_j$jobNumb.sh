#!/bin/sh
/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/scriptcombine /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/datacardWithAllSyst/HighMass-hgg_8TeV_m150.00_w10_channel0_GGH.txt -M Asymptotic -m 150 -D data_obs --run=both  -U --noFitAsimov --saveToys -S 1 -n _w10
