#!/usr/bin/bash 
SEED=2023 

function update_seed(){
    local n=$1
    local result=$(((n*92341)%10000))
    echo $result
}
####### Null Experiments 
#python3 ExperimentNull.py --pdf_plot --save_fig -n 200 -r 500 -d 10 --seed 1234
#python3 ExperimentNull.py --pdf_plot  --save_fig -n 200 -r 500 -d 100 --seed 1324

#python3 ExperimentNull.py --pdf_plot --save_fig -n 200 -r 500 -d 10 -k RationalQuadratic --seed 3421
#python3 ExperimentNull.py --pdf_plot  --save_fig -n 200 -r 500 -d 100 -k RationalQuadratic --seed 3142

#python3 ExperimentNull.py --pdf_plot --save_fig -n 200 -r 500 -d 10 -exp dCov --seed 2341
#python3 ExperimentNull.py --pdf_plot  --save_fig -n 200 -r 500 -d 100 -exp dCov --seed 2134


#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 10 -exp HSIC --seed 5678
#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 100 -exp HSIC --seed 5876

#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 10 -k RationalQuadratic --seed 6875
#python3 ExperimentNull.py --typeI_plot  --save_fig -n 200 -r 500 -d 100 -k RationalQuadratic --seed 6587

#python3 ExperimentNull.py --typeI_plot --save_fig -n 200 -r 500 -d 10 -exp dCov --seed 7856
#python3 ExperimentNull.py --typeI_plot  --save_fig -n 200 -r 500 -d 100 -exp dCov --seed 7568

#python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 10 -exp HSIC --seed 1289
#python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 10 -exp dCov --seed 1982

#python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 100 -exp HSIC --seed 2189
#python3 ExperimentNull.py --normality_test  --save_fig -n 200 -r 200 -nt 250 -d 100 -exp dCov --seed 2918


SEED=$(update_seed $SEED)
python3 ExperimentNull.py --second_moment -df 1  --save_fig -n 500 -r 500 -d 10 -k Linear --seed $SEED 
SEED=$(update_seed $SEED)
python3 ExperimentNull.py --second_moment -df 2  --save_fig -n 500 -r 500 -d 10 -k Linear --seed $SEED 
SEED=$(update_seed $SEED)
python3 ExperimentNull.py --second_moment -df 3  --save_fig -n 500 -r 500 -d 10 -k Linear --seed $SEED 
SEED=$(update_seed $SEED)

python3 ExperimentNull.py --second_moment -df 1  --save_fig -n 500 -r 500 -d 10 -exp dCov --seed $SEED 
SEED=$(update_seed $SEED)
python3 ExperimentNull.py --second_moment -df 2  --save_fig -n 500 -r 500 -d 10 -exp dCov --seed $SEED 
SEED=$(update_seed $SEED)
python3 ExperimentNull.py --second_moment -df 3  --save_fig -n 500 -r 500 -d 10 -exp dCov --seed $SEED 
SEED=$(update_seed $SEED)


