## Permutation-free kernel independence test 

This repository contains the code for reproducing the plots in the paper "[A permutation free kernel independence test](https://arxiv.org/abs/2212.09108)".  We introduce a new unbiased estimate of Hilbert-Schmidt Independence Criterion~(HSIC) (that we call the cross-HSIC statistic) that has a standard normal limiting null distribution in different dimension regimes. This is unlike the usual quadratic-time HSIC statistic that is a degenerate U-statistic, and thus has an intractable limiting null distribution. The simple null distribution allows us to propose a permutation-free independence test, that leads to a significant reduction in computation~(usually 100x) at the price of a small reduction~($\approx \sqrt{2}$ factor) in power.  

### Reproducing the results
* Install the dependencies: `
```
pip install -r requirements.txt`
```
* To reproduce the figures run the scripts 

```
# generate the figures under the null
run_alt_experiments.sh

# generate the figures under the alternative 
run_null_experiments.sh
```

### Structure 
* The files `ExperimentNull.py` and `ExperimentAlt.py` contain the main implementations of the experiments under the null and alternative respectively. 
* The file `independence_tests.py` contains the following tests: 
    * cross-HSIC test 
    * cross-dCov test 
    * HSIC permutation test 
    * HSIC dCov test 
* The file `utils.py` contains some helper functions 
* The file `crossHSIC.py` implements the cross-HSIC statistic 