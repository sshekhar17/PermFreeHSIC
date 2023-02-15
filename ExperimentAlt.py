import multiprocess as mp 
import os
import argparse
from functools import partial 
import numpy as np 
from independence_tests import * 
import matplotlib.pyplot as plt 
from utils import CreateDependentGaussianSource, initialize_kernel
import tikzplotlib as tpl 
from time import time 
from constants import ColorsDict, markersDict


def compare_power(Source=None, 
                        tests_dict = None, 
                        kwargs_dict=None, 
                        N=150, #the maximum sample-size 
                        num_trials=20, # number of repetitions 
                        alpha=0.05, # size of the test
                        initial_val=10, # starting value of sample-size
                        num_steps=15,# number of steps from initial_val to N
                        save_fig=False, 
                        fig_dir = './figures/alt', 
                        figname=None, 
                        parallel=False, 
                        seed=None, 
                        plot_time_comparison=False): 
    """
        Compare the power of different tests of independence 

        #Parameters 
            Source      :function handle    
                            takes in a positive integer n, and returns two numpy 
                            arrays X and Y of sizes (n, d_X) and (n, d_Y)
            tests_dict  :dictionary 
                            keys: strings denoting the names of tests 
                            values: function-handles to run the tests indicated by keys
            kwargs_dict :dictionary 
                            keys: same as tests_dict 
                            values: dicts containing kwargs to be sent to the corresponding 
                                    functions contained in tests_dict 
            N           :int 
                            maximum-sample size 
            initial_val :int 
                            starting value of sample-size grid 
            num_steps   :int 
                            number of points in the grid of sample sizes 
                            between [initial_val, N]
            alpha       :float 
                            significance level of the test 
            num_trials  :int
                            number of trials to estimate the power and type-I errors 
            save_fig    :bool 
                            if true save the figure as a tikzplot figure 
    """

    if Source is None:
        # default is a dependent Gaussian source in 10 dimensions 
        Source = CreateDependentGaussianSource()
    
    if tests_dict is None:
        tests_dict = {'x-HSIC':crossHSIC_test, 'HSIC-perm':HSIC_permutation_test} 
        kwargs_dict = {'x-HSIC':{}, 'HSIC-perm':{'num_perms':50}} 

    NN = np.linspace(initial_val, N, num_steps, dtype=int)
    PowerDict, TimesDict = {}, {}
    for test in tests_dict:
        PowerDict[test] = np.zeros(NN.shape)
        TimesDict[test] = np.zeros(NN.shape)

    if seed is not None:
        np.random.seed(seed) 

    if not parallel:
        for trial in tqdm(range(num_trials)):
            for i, n in enumerate(NN):
                if n%2: #only even values work 
                    n+=1 
                X, Y = Source(n)
                for test in tests_dict:
                    t0 = time()
                    test_func = tests_dict[test]
                    kwargs = kwargs_dict[test]
                    # accumulate the number of rejections 
                    PowerDict[test][i] += test_func(XX=X, YY=Y, alpha=alpha, **kwargs)
                    time_taken = time() - t0 
                    TimesDict[test][i] += time_taken 
    else:
        start_time = time()
        seed0 = seed if seed is not None else 1234
        # define the helper function to be used in parallel
        def helper(trial):
            # set the random seed 
            np.random.seed(((1+trial)*seed0*3811)%10000)
            # initialize the local power dictionary 
            PowerDict_ = {}
            for test in tests_dict:
                PowerDict_[test] = np.zeros(NN.shape)


            for i, n in enumerate(NN):
                if n%2: #only even values work 
                    n+=1 
                X, Y = Source(n)
                for test in tests_dict:
                    test_func = tests_dict[test]
                    kwargs = kwargs_dict[test]
                    # accumulate the number of rejections 
                    PowerDict_[test][i] += test_func(XX=X, YY=Y, alpha=alpha, **kwargs)           
            return PowerDict_
        # run helper in parallel
        n_cores = mp.cpu_count()
        with mp.Pool(n_cores) as PP: 
            result = PP.map(helper, range(num_trials))
        for r in result:
            for test in PowerDict:
                PowerDict[test] += r[test]
        run_time = time() - start_time 
        print(f'Parallel processing took {run_time:.2f} seconds')
    # Power = (# of rejections) / (num_trials)
    for test in tests_dict:
        PowerDict[test] /= num_trials    
        TimesDict[test] /= num_trials

    ## Plot power vs sample-size 
    if not plot_time_comparison:
        plt.figure() 
        for test in tests_dict: 
            plt.plot(NN, PowerDict[test], label=test, color=ColorsDict[test])

        plt.xlabel('Sample-size (n)', fontsize=13)
        plt.ylabel('Power', fontsize=13)
        plt.title('Power vs Sample Size', fontsize=15)
        plt.legend()
        # check if a figure name is provided
        if figname is None: 
            figname = 'temp_fig_power'
        # save the figure if requrired 
        if save_fig:
            figname_ = f'{fig_dir}/{figname}'
            plt.savefig(figname_+'.png', dpi=450)
            tpl.save(figname_+'.tex', axis_width=r'\figwidth', axis_height=r'\figheight')
        else:
            plt.show()
    # plot the time-vs-power comparison  
    if plot_time_comparison:
        sizes = (NN/N)*100
        plt.figure()
        for test in tests_dict: 
            times_, power_ = TimesDict[test], PowerDict[test]
            plt.scatter(times_, power_, s=sizes, label=test, alpha=0.6, 
                            color=ColorsDict[test], marker=markersDict[test])

        plt.xscale('log')
        plt.xlabel('Running Time (Seconds)', fontsize=13)
        plt.ylabel('Power', fontsize=13)
        plt.title('Power vs Computation', fontsize=15)
        plt.legend()

        if save_fig:
            figname_ = f'{fig_dir}/{figname}_Running_Times'
            plt.savefig(figname_+'.png', dpi=450)
            tpl.save(figname_+'.tex', axis_width=r'\figwidth', axis_height=r'\figheight')
        else:
            plt.show()


def get_args():
    parser = argparse.ArgumentParser() 
    parser.add_argument('-d', '--ndims', dest='d', default=10, type=int) 
    parser.add_argument('-n', '--n', dest='n', default=200, type=int)
    parser.add_argument('-a', '--alpha', default=0.05, type=float) 
    parser.add_argument('-ns', '--num_steps', default=20, type=int)
    parser.add_argument('-iv', '--initial_val', default=10, type=int)
    parser.add_argument('-nt', '--num_trials', default=250, type=int)
    parser.add_argument('-np', '--num_perms', default=150, type=int)
    parser.add_argument('-exp', '--expt_type', default='HSIC', choices=['HSIC', 'dCov'])
    parser.add_argument('-k', '--kernel_type', default=None, 
                            choices=[None, 'RBF', 'Linear', 'Polynomial', 'Exponential',
                                        'RationalQuadratic'])
    parser.add_argument('--pX', default=2, type=int)
    parser.add_argument('--pY', default=2, type=int)
    parser.add_argument('-pS', '--pSource', default=2, type=int)
    parser.add_argument('-e', '--epsilon', default=0.5, type=float)
    parser.add_argument('--progress_bar', action='store_true')
    parser.add_argument('-p', '--parallel', action='store_true')
    parser.add_argument('--save_fig', action='store_true')
    parser.add_argument('-s', '--seed', default=None, type=int)
    parser.add_argument('--fig_dir', default='./figures/alt')
    parser.add_argument('-t','--plot_time_comparison', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.fig_dir):
        os.makedirs(args.fig_dir)
    if args.seed is not None:
        np.random.seed(args.seed)
    return args

if __name__=='__main__':
    args = get_args()
    # args = None 
    if args is None:
        d = 10 
        n = 900
        alpha = 0.05
        num_steps = 5 
        initial_val = 10
        num_trials = 20
        num_perms=50
        expt_type='HSIC' 
        kernel_type = None # default RBF + median heuristic 
        pX, pY = 2, 2
        epsilon = 0.3 # the level of dependence
        pSource = 2
        progress_bar = True 
        save_fig = False 
        fig_dir = './figures/alt'
        parallel = False
        plot_time_comparison=True
    else:
        d = args.d 
        n = args.n
        alpha= args.alpha
        num_steps = args.num_steps
        initial_val = args.initial_val
        num_trials = args.num_trials
        num_perms = args.num_perms
        expt_type = args.expt_type # other option is 'dCov' 
        kernel_type= args.kernel_type # other options: Linear, Polynomial
        pX, pY = args.pX, args.pY # default values 
        # parameters of the source 
        epsilon = args.epsilon  # denotes the level of dependence; set to 0 under the null
        pSource = args.pSource  # exponent of func used by 'CreateDependentGaussianSource'
        # the experiments to perform 
        progress_bar = args.progress_bar 
        save_fig = args.save_fig 
        fig_dir = args.fig_dir
        parallel = args.parallel
        plot_time_comparison=args.plot_time_comparison

    #####
    if expt_type=='HSIC':
        temp = f'kernel_{kernel_type}'
    else:
        temp = f'pX_{pX}_pY_{pY}'
    print('='*80 + '\n')
    print(rf'Starting {expt_type} Experiment: $\epsilon$={epsilon}'+temp)
    print('='*80 + '\n')
    ####

    Source = CreateDependentGaussianSource(epsilon=epsilon, ndims=d, p=pSource)
    figname = None 
    #### Set up the statistical tests 
    if expt_type == 'HSIC':
        kernel_X, kernel_Y = initialize_kernel(kernel_type)
        tests_dict = {'HSIC-perm':HSIC_permutation_test, 'x-HSIC':crossHSIC_test} 

        kwargs_dict = {'HSIC-perm':{'kernel_X':kernel_X, 'kernel_Y':kernel_Y, 'num_perms':num_perms},
                        'x-HSIC':{'kernel_X':kernel_X, 'kernel_Y':kernel_Y}} 
        if kernel_type is None:
            kernel_type = 'RBF+Median'
        if save_fig:
            figname = f'Power_HSIC_epsilon_{epsilon}_n_{n}_d_{d}_kernel_{kernel_type}'
            figname = figname.replace('.', '_')
    else:
        dcov_stat_func = partial(distance_covariance_test, return_stat=True, pX=pX, pY=pY)
        tests_dict = {'dcov-perm':permutation_test, 'x-dcov':distance_covaraince_DA_test} 
        kwargs_dict = {'dcov-perm':{'stat_func':dcov_stat_func, 'num_perms':num_perms}, 
                    'x-dcov':{'pX':pX, 'pY':pY}
                    } 
        if save_fig:
            figname = f'Power_dCov_epsilon_{epsilon}_n_{n}_d_{d}_pX_{pX}_pY_{pY}'
            figname = figname.replace('.', '_')
    ######################################################################## 
    ### Run the experiment 
    compare_power(Source,
                    N=n,
                    num_trials=num_trials,
                    tests_dict=tests_dict,
                    kwargs_dict=kwargs_dict,
                    initial_val=initial_val, 
                    num_steps=num_steps, 
                    save_fig=save_fig, 
                    figname=figname, 
                    parallel=parallel, 
                    plot_time_comparison=plot_time_comparison
                )
    ######################################################################## 