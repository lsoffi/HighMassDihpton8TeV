#!/bin/bash

for ((mass=150;mass<500;mass++))#mass in   150 160 170 180 190 210 220 230 240 260 270 280 290 300 310 320 330 340 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 640 650 660 670 680 690 700 710 720 730 740 760 770 780 790 810 820 830 840 850
  do
  for width in 5 #0.1 2 5 10
    do
    for model in  BwSyst
      do
      cat>>mk_workspace.C<<EOF

{
gROOT->ProcessLine(".L ProduceWorkspaces.C");
ProduceWorkspaces($mass, $width, "${model}");
}


EOF
      cat mk_workspace.C
      echo $mass $width $model
      root -l -b -q mk_workspace.C > log_ws_$mass_$width_$model.log
      rm mk_workspace.C 
      
    done #width
  done #mass
done #model



#std masses:
#  150 250 300 350 400 450 500 550 600 650 700 750 800 850 200 300 400 500 600 700

#done1
#150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 320 160 170 180 190 210 220 

#done1
#230 240 260 270 280 290 310 330 340 360 370 380 390 410 420 430 440 460 470 480 490 510
