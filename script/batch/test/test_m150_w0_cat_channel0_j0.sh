#!/bin/sh
/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/script 
eval export PATH="/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/bin/slc5_amd64_gcc472:/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/external/slc5_amd64_gcc472/bin:/afs/cern.ch/cms/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_1/bin/slc5_amd64_gcc472:/afs/cern.ch/cms/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_1/external/slc5_amd64_gcc472/bin:/afs/cern.ch/cms/slc5_amd64_gcc472/external/llvm/3.2-cms2/bin:/afs/cern.ch/cms/slc5_amd64_gcc472/external/gcc/4.7.2/bin:/afs/cern.ch/user/s/soffi/bin:/afs/cern.ch/user/s/soffi/scripts:/usr/sue/bin:/afs/cern.ch/group/zh/bin:/usr/local/bin:/usr/local/bin/X11:/usr/bin:/bin:/usr/bin/X11:/cern/pro/bin:/afs/cern.ch/cms/caf/scripts:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/bin:/usr/lib64/qt-3.3/bin:/usr/kerberos/bin";
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/datacardWithAllSyst/HighMass-hgg_8TeV_m150.00_w0.10_channel0_GGH.txt ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/HighMass-hgg_8TeV_m150.00_w0.10.inputsig_GGH.root ./
cp  /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/HighMass-hgg.inputbkg_m150.00.root ./

combine HighMass-hgg_8TeV_m150.00_w0.10_channel0_GGH.txt -M Asymptotic -m 150 -D data_obs --run=both  -U  -S 0 -n _w0.10
mv -v higgs*root /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/channel0/
