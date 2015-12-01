
#!/bin/bash
#append=local

usage(){
    echo "Usage: `basename $0` joblist_file" > /dev/stderr
}

case $# in 
    0)
	echo -n "Error: " > /dev/stderr
	usage
	exit 1
	;;
    *)
	;;
esac

taskName=mkWS

nJobs=1
queue=$1
let nJobs_=$nJobs-1
model=$2

mkdir -p log/
mkdir -p log/$taskName
mkdir -p batch/
mkdir -p batch/${taskName}
echo "################ Task ${taskName} ##############"    

#DiJet


     for jobNumb in `seq 0 $nJobs_`
	  do
	  for  ((mass=712;mass<=850;mass+=2)) #
	    do
	    for width in  7.00 #10.00 7.00 5.00 2.00 0.10
	      do
	      
		echo "Processing jobNumb: $jobNumb,  mass: $mass, width: $width cat: $cat" > /dev/stdout    
				
		sample='m'$mass'_w'$width'_mod'$model
		echo $sample
		jobname=$taskName"_"$sample
		stdout_file=`pwd`/log/$taskName/`basename ${sample}_ .list`-out.log
		
		    
		    cat <<EOF >  batch/$taskName/${jobname}.sh
#!/bin/sh
cd $PWD
eval \`scramv1 runtime -sh\`
cd \$WORKDIR
mkdir workspaces
mkdir datacardWithAllSyst
mkdir plots
mkdir preliminaryPlots
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/macroWS/mk_workspace_${mass}_${width}_${model}.C ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/HighMass-hgg_models_Bkg_8TeV_test.rs ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/ProduceWorkspaces.C ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/HighMass-HggFitter_mgg.cc ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/histograms_CMS-HGG_08052014_MC.root ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/histograms_CMS-HGG_17042014_DATA.root ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/sigShapeCorrections.root ./
cp /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/effXacc2D*.root ./

root -l -b -q mk_workspace_${mass}_${width}_${model}.C > log_ws_${mass}_${width}_${model}.log
mv -v workspaces/HighMass-hgg_8TeV_m${mass}.00_w${width}.inputsig_${model}.root /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/
mv -v workspaces/HighMass-hgg.inputbkg_m${mass}.00.root /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/
mv -v datacardWithAllSyst/HighMass-hgg_8TeV_m${mass}.00_w${width}_channel*_${model}.txt /afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/datacardWithAllSyst/

EOF
		    echo "bsub -J $jobname -q $queue -o $stdout_file < `pwd`/batch/$taskName/${jobname}.sh"
		 bsub -J $jobname -q $queue -o $stdout_file < `pwd`/batch/$taskName/${jobname}.sh
		
	     
	    done
	 done
    done


#std masses:
#150 200 250 300 350 400 450 500 550 600 650 700 750 800 850
#


#done1
# 
#done2
#


#556 558 560 562 564 566 568 570 572 574 576 578 580 582 584 586 588 590 592 594 596 598 600 602 604 #
