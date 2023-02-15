#!/usr/bin/bash
SEED=2024

function update_seed(){
    local n=$1
    local result=$(((n*92341)%10000))
    echo $result
}

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 300 -e 0.5 -s $SEED

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 600 -e 0.4 -s $SEED

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -n 900 -e 0.3 -s $SEED


# SEED=$(update_seed $SEED)
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 300 -e 0.5 -s $SEED

# SEED=$(update_seed $SEED)
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 600 -e 0.4 -s $SEED

# SEED=$(update_seed $SEED)
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 900 -e 0.3 -s $SEED


SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 300 -e 0.5 -s $SEED

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 600 -e 0.4 -s $SEED

SEED=$(update_seed $SEED)
python3 ExperimentAlt.py --save_fig --progress_bar -exp HSIC -k RationalQuadratic -n 900 -e 0.3 -s $SEED


SEED=$(update_seed $SEED)
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 300 -e 0.5 -pX 4 -pY 3 -s $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 600 -e 0.4 -pX 4 -pY 3 -s $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentAlt.py --save_fig --progress_bar -exp dCov -n 900 -e 0.3 -pX 4 -pY 3 -s $SEED
