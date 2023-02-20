#!/usr/bin/bash 
SEED=2023 

function update_seed(){
    local n=$1
    local result=$(((n*92341)%10000))
    echo $result
}
####### Null Experiments 
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --pdf_plot --save_fig -n 200 -r 500 -d 10 --seed $SEED
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --pdf_plot  --save_fig -n 200 -r 500 -d 100 --seed $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --pdf_plot --save_fig -n 200 -r 500 -d 10 -k RationalQuadratic --seed $SEED
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --pdf_plot  --save_fig -n 200 -r 500 -d 100 -k RationalQuadratic --seed $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --pdf_plot --save_fig -n 200 -r 500 -d 10 -exp dCov --seed $SEED
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --pdf_plot  --save_fig -n 200 -r 500 -d 100 -exp dCov --seed $SEED


SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 10 -exp HSIC --seed $SEED
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 100 -exp HSIC --seed $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 10 -k RationalQuadratic --seed $SEED
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --typeI_plot  --save_fig -n 200 -r 500 -d 100 -k RationalQuadratic --seed $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 10 -exp dCov --seed $SEED
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --typeI_plot  --save_fig -n 200 -r 500 -d 100 -exp dCov --seed $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 10 -exp HSIC --seed $SEED
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 10 -exp dCov --seed $SEED

SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 100 -exp HSIC --seed $SEED
SEED=$(update_seed $SEED)
python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 100 -exp dCov --seed $SEED


SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --second_moment -df 1  --save_fig -n 500 -r 500 -d 10 -k Linear --seed $SEED 
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --second_moment -df 2  --save_fig -n 500 -r 500 -d 10 -k Linear --seed $SEED 
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --second_moment -df 3  --save_fig -n 500 -r 500 -d 10 -k Linear --seed $SEED 
SEED=$(update_seed $SEED)

#python3 ExperimentNull.py --second_moment -df 1  --save_fig -n 500 -r 500 -d 10 -exp dCov --seed $SEED 
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --second_moment -df 2  --save_fig -n 500 -r 500 -d 10 -exp dCov --seed $SEED 
SEED=$(update_seed $SEED)
#python3 ExperimentNull.py --second_moment -df 3  --save_fig -n 500 -r 500 -d 10 -exp dCov --seed $SEED 
SEED=$(update_seed $SEED)


