#!/usr/bin/bash
SEED=2024

function update_seed(){
    local n=$1
    local result=$(((n*92341)%10000))
    echo $result
}

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 400 -e 0.5 -s $SEED 

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 800 -e 0.4 -s $SEED  

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 1200 -e 0.3 -s $SEED 


SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 400 -e 0.5 -s $SEED  

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 800 -e 0.4 -s $SEED  

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 1200 -e 0.3 -s $SEED  



SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 400 -e 0.5 -s $SEED  

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 800 -e 0.4 -s $SEED  

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 1200 -e 0.3 -s $SEED  


SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 40 -e 0.5 -s $SEED -t 

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 80 -e 0.4 -s $SEED  -t 

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 120 -e 0.3 -s $SEED -t 


SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 40 -e 0.5 -s $SEED  -t 

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 80 -e 0.4 -s $SEED  -t 

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 120 -e 0.3 -s $SEED  -t 



SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 400 -e 0.5 -s $SEED  -t 

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 800 -e 0.4 -s $SEED  -t 

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 1200 -e 0.3 -s $SEED  -t 


