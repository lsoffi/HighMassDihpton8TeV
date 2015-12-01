#!/bin/bash

for mass in  173 174 179 185 196 198 201 216 22 226 228 246 261 265 292 301 319 321 322 323 327 351 353 359 365 375 377 396 401 413 414 418 419 422 423 428 435 436 444 462 476 501 521 539 616 625 626 601 614 639 692 701 709 743 754 763 784 768 787
  do
  for width in 0.1 #0.1 2 5 10
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
