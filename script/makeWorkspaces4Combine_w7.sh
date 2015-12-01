#!/bin/bash

for mass in 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 #200 300 400 500 600 700  
  do
  for width in  7 
    do
    cat>>mk_workspace.C<<EOF

{
gROOT->ProcessLine(".L ProduceWorkspaces.C");
ProduceWorkspaces($mass, $width);
}


EOF
    
    echo $mass $width 
    root -l -b -q mk_workspace.C > log_ws_$mass_$width.log
    rm mk_workspace.C 

  done #width
done #mass


