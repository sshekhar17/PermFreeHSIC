#!/usr/bin/bash

python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 300 -e 0.5 -s 1234 
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 600 -e 0.4 -s 5148 
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 900 -e 0.3 -s 9256 

#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 300 -e 0.5 -s 1234 
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 600 -e 0.4 -s 5148 
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 900 -e 0.3 -s 9256 


python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 300 -e 0.5 -s 7324 
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 600 -e 0.4 -s 5487 
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 900 -e 0.3 -s 3392 

#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 300 -e 0.5 -pX 4 -pY 3 -s 7324 
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 600 -e 0.4 -pX 4 -pY 3 -s 5487 
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 900 -e 0.3 -pX 4 -pY 3 -s 3392 
