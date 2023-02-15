"""
Experiment evaluating the null distribution of the statistic. 
"""
import os
import argparse
from math import sqrt
from functools import partial
from crossHSIC import get_studentized_cross_dcov, get_studentized_cross_hsic
from utils import LinearKernel, PolynomialKernel, RBFkernel, ExponentialKernel, RationalQuadraticKernel 
from utils import CreateDependentGaussianSource, CreateIndependentGaussianSource
from utils import CreateIndependentTSource
from utils import initialize_kernel

from time import time 
from tqdm import tqdm
import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import cdist, pdist 
import matplotlib.pyplot as plt
import tikzplotlib as tpl


def compute_stats(Source, n, stat_func, stat_params={}, repetitions=100, 
                    progress_bar=False): 

    StatVals = np.zeros((repetitions,)) 
    range_ = tqdm(range(repetitions)) if progress_bar else range(repetitions) 

    for trial in range_: 
        # generate the data 
        X, Y = Source(n) 
        # calculate the statistic 
        StatVals[trial] = stat_func(X, Y, **stat_params) 
    return StatVals 

def typeI_error(Source, NN, stat_func, stat_params={}, repetitions=100, 
                    progress_bar=False, alpha=0.05):
    typeI = np.zeros(NN.shape) 
    range_ = tqdm(list(enumerate(NN))) if progress_bar else list(enumerate(NN))
    th = stats.norm.ppf(1-alpha) 
    for i, n in range_: 
        stat_vals = compute_stats(Source, n, stat_func, stat_params, repetitions) 
        typeI[i] = np.sum(stat_vals>th)/repetitions 
    return typeI 


def test_for_normality(Source, NN, stat_func, stat_params={}, 
                        repetitions=100, progress_bar=False, 
                        alpha=None, num_trials=250):
    alpha = 0.05 if alpha is None else alpha 
    Power = np.zeros(NN.shape)
    range_ = range(num_trials)
    range_ = tqdm(range_) if progress_bar else range_
    for _ in range_:
        for i, n in enumerate(NN):
            stat_vals = compute_stats(Source, n, stat_func, stat_params, repetitions) 
            _, pval = stats.normaltest(stat_vals)
            Power[i] += 1.0*(pval<=alpha)
    # normalize the power 
    Power /= num_trials
    return Power


def get_args():
    parser = argparse.ArgumentParser() 
    parser.add_argument('-d', '--ndims', dest='d', default=10, type=int) 
    parser.add_argument('-dh', '--d_high', dest='d_high', default=150, type=int) 
    parser.add_argument('-n', '--n', dest='n', default=200, type=int)
    parser.add_argument('-a', '--alpha', default=0.05, type=float) 
    parser.add_argument('-ns', '--num_steps', default=20, type=int)
    parser.add_argument('-iv', '--initial_val', default=10, type=int)
    parser.add_argument('-r', '--repetitions', default=100, type=int, 
                         help="number of times to recompute the statistic")
    parser.add_argument('-nt', '--num_trials', default=100, type=int,
                         help="number of time the 'outer' experiment is run")
    parser.add_argument('-exp', '--expt_type', default='HSIC', choices=['HSIC', 'dCov'])
    parser.add_argument('-k', '--kernel_type', default=None, 
                            choices=[None, 'RBF', 'Linear', 'Polynomial', 'Exponential',
                                        'RationalQuadratic'])
    parser.add_argument('--pX', default=2, type=int)
    parser.add_argument('--pY', default=2, type=int)
    parser.add_argument('-pS', '--pSource', default=2, type=int)
    parser.add_argument('-e', '--epsilon', default=0.0)
    parser.add_argument('--pdf_plot', action='store_true')
    parser.add_argument('--typeI_plot', action='store_true')
    parser.add_argument('--normality_test', action='store_true')
    parser.add_argument('--second_moment', action='store_true')
    parser.add_argument('-df', '--df', default=3, type=int, help='deg of freedom in T dist')
    parser.add_argument('--progress_bar', action='store_true')
    parser.add_argument('--save_fig', action='store_true')
    parser.add_argument('-s', '--seed', default=None, type=int)
    parser.add_argument('--fig_dir', default='./figures/null')
    args = parser.parse_args()

    # fig_dir = os.path.dirname(args.fig_dir)
    ## Note that this assumes that the directory path 
    ## does not have the "/" at the end 
    if not os.path.exists(args.fig_dir):
        os.makedirs(args.fig_dir)
    if args.seed is not None:
        np.random.seed(args.seed)
    return args


if __name__=='__main__': 
    ######################################
    ### all the experiment arguments here 
    args = get_args()
    # args = None 
    if args is not None:
        d = args.d 
        d_high = args.d_high
        n = args.n
        alpha= args.alpha
        num_steps = args.num_steps
        initial_val = args.initial_val
        repetitions = args.repetitions
        num_trials = args.num_trials
        expt_type = args.expt_type # other option is 'dCov' 
        kernel_type= args.kernel_type # other options: Linear, Polynomial
        pX, pY = args.pX, args.pY # default values 
        # parameters of the source 
        epsilon = args.epsilon  # denotes the level of dependence; set to 0 under the null
        pSource = args.pSource  # exponent of func used by 'CreateDependentGaussianSource'
        # the experiments to perform 
        pdf_plot=args.pdf_plot # pdf plot
        typeI_plot=args.typeI_plot  # type-I error plot 
        normality_test= args.normality_test # p values plot 
        second_moment= args.second_moment
        df = args.df
        progress_bar = args.progress_bar 
        save_fig = args.save_fig 
        fig_dir = args.fig_dir
    else:
        d = 10
        d_high = 150
        n = 150
        alpha= 0.05
        num_steps = 15
        initial_val = 10
        repetitions = 500
        num_trials = 20
        expt_type = 'HSIC'# other option is 'dCov' 
        kernel_type= 'RBF+median'# other options: Linear, Polynomial
        pX, pY = 2, 2 # default values 
        # parameters of the source 
        epsilon = 0  # denotes the level of dependence; set to 0 under the null
        pSource = 2  # exponent of func used by 'CreateDependentGaussianSource'
        # the experiments to perform 
        pdf_plot= True 
        typeI_plot= False
        normality_test= False
        second_moment = False
        df = 3 # degrees of freedom in the T-distribution 
        progress_bar = True
        save_fig = False
        fig_dir = './figures'

    #####
    if args is not None:
        if args.seed is None: 
            seed = int(time()%10000)
        else:
            seed = args.seed
        print('\n')
        print(f'Starting Experiment: seed = {seed}')
        print(f'{expt_type}:\t (PDF:{pdf_plot}), \t (TypeI:{typeI_plot}, \t (Normality:{normality_test})')
        print('\n')

    ###### Set up the functions and parmas for computing the statistics 
    if expt_type == 'HSIC':
        # initialize the kernels  for the tests 
        kernel_X, kernel_Y = initialize_kernel(kernel_type)
        ### initialize the stat_func and params 
        stat_func = get_studentized_cross_hsic 
        stat_params = {'kernel_X':kernel_X, 'kernel_Y':kernel_Y}
    elif expt_type == 'dCov':
        stat_func = get_studentized_cross_dcov 
        pX = 2 if pX is None else pX
        pY = 2 if pY is None else pY
        stat_params = {'pX':pX, 'pY':pY}
    else:
        raise Exception(f'Experiment type {expt_type} not recognized')

    ##### Create the source 
    #### TODO: allow other types of sources to be used 
    #### epsilon=0 corresponds to independent sources 
    Source = CreateDependentGaussianSource(ndims=d,
                            epsilon=epsilon, 
                            p=pSource)
    SourceHigh = CreateDependentGaussianSource(ndims=d_high,
                            epsilon=epsilon, 
                            p=pSource)

    #### Plot the empirical pdf  
    if pdf_plot:
        Stat = compute_stats(Source,n, stat_func, stat_params, repetitions=repetitions,
                                progress_bar=progress_bar)
        Stat2 = compute_stats(SourceHigh,n, stat_func, stat_params, repetitions=repetitions,
                                progress_bar=progress_bar)

        #### 
        if expt_type=='HSIC':
            title = f'{kernel_type} (n={n})' 
        else:
            title = f'{expt_type} (n={n})'
        label = f'x-{expt_type}'
        ylabel = 'Density'
        xx = np.linspace(-10, 10, 1000)
        yy = stats.norm.pdf(xx) 
        # xx = xx + Stat.mean()
        #######################
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(x=[Stat, Stat2], density=True, label=[f'd={d}', f'd={d_high}'], alpha=0.8)
        ax.plot(xx, yy, label='N(0,1)', color='k') 
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title, fontsize=16)
        ax.legend(fontsize=14, fancybox=True, framealpha=0.0)
        if save_fig:
            if epsilon==0:
                figname = f'NullPdf_{expt_type}_n_{n}_d_{d}_d_high_{d_high}_seed_{seed}'
            else:
                figname = f'AltPdf_{expt_type}_n_{n}_d_{d}_d_high_{d_high}_epsilon_{epsilon}_seed_{seed}'.replace('.','_')
            if expt_type=='HSIC' and kernel_type is not None:
                figname = f'{figname}_{kernel_type}'
            figname_tex = f'{fig_dir}/{figname}.tex' 
            figname_png = f'{fig_dir}/{figname}.png' 
            plt.savefig(figname_png, dpi=450)
            tpl.save(figname_tex, axis_width=r'\figwidth', axis_height=r'\figheight')
        else:
            plt.show()
    ##### Type-I error 
    if typeI_plot:
        NN = np.linspace(initial_val, n, num=num_steps, dtype=int)
        NN = np.array([(n_//2)*2 for n_ in NN], dtype=int)
        TypeI = typeI_error(Source, NN, stat_func, stat_params, repetitions, 
                            progress_bar=progress_bar, alpha=alpha)

        plt.figure()
        plt.plot(NN, TypeI)
        plt.plot(NN, np.ones(NN.shape)*alpha, 'k', alpha=0.15) 
        plt.ylim([0, 2.0*max(TypeI)])
        plt.title(f'Type-I error vs Sample Size ({expt_type})', fontsize=15)
        plt.xlabel('Sample size (n)', fontsize=13)
        plt.ylabel('Type-I error', fontsize=13)
        if save_fig:
            if epsilon==0:
                figname = f'TypeI_error_{expt_type}_n_{n}_d_{d}_seed_{seed}'
            else:
                figname = f'Power_{expt_type}_n_{n}_d_{d}_epsilon_{epsilon}_seed_{seed}'.replace('.','_')
            if expt_type=='HSIC' and kernel_type is not None:
                figname = f'{figname}_{kernel_type}'
            else:
                figname = f'{figname}_{expt_type}'

            figname_tex = f'{fig_dir}/{figname}.tex' 
            figname_png = f'{fig_dir}/{figname}.png' 
            plt.savefig(figname_png, dpi=450)
            tpl.save(figname_tex, axis_width=r'\figwidth', axis_height=r'\figheight')
        else:
            plt.show()

    ##### Test for Normality 
    if normality_test:
        NN = np.linspace(initial_val, n, num=num_steps, dtype=int)
        NN = np.array([(n_//2)*2 for n_ in NN], dtype=int)
        Power = test_for_normality(Source, NN, stat_func, stat_params,  
                                        repetitions=repetitions, progress_bar=True, 
                                        num_trials=num_trials)
        plt.figure()
        plt.plot(NN, Power)
        plt.plot(NN, np.ones(NN.shape)*alpha, 'k', alpha=0.15) 
        plt.ylim([0, 1.0])
        plt.title('Power of the test for normality', fontsize=15)
        plt.ylabel('Power', fontsize=13)
        plt.xlabel('Sample Size (n)', fontsize=13)
        if save_fig:
            if epsilon==0:
                figname = f'Null_Normality_test_{expt_type}_n_{n}_d_{d}_seed_{seed}'
            else:
                figname = f'Alt_Normaity_test_{expt_type}_n_{n}_d_{d}_epsilon_{epsilon}_seed_{seed}'.replace('.','_')
            if expt_type=='HSIC' and kernel_type is not None:
                figname = f'{figname}_{kernel_type}'
            else:
                figname = f'{figname}_{expt_type}'

            figname_tex = f'{fig_dir}/{figname}.tex' 
            figname_png = f'{fig_dir}/{figname}.png' 
            plt.savefig(figname_png, dpi=450)
            tpl.save(figname_tex, axis_width=r'\figwidth', axis_height=r'\figheight')
        else:
            plt.show()

    if second_moment:
        # create a source with finite second moment / but infinite third and higher moments
        Source = CreateIndependentTSource(df1=df, df2=df, ndims1=d, ndims2=d)
        Source2 = CreateIndependentTSource(df1=df, df2=df, ndims1=d, ndims2=d_high)
        kernel_type = 'Linear'
        ###### Set up the functions and parmas for computing the statistics 
        if expt_type == 'HSIC':
            # initialize the kernels  for the tests 
            kernel_X, kernel_Y = initialize_kernel(kernel_type)
            ### initialize the stat_func and params 
            stat_func = get_studentized_cross_hsic 
            stat_params = {'kernel_X':kernel_X, 'kernel_Y':kernel_Y}
        elif expt_type == 'dCov':
            stat_func = get_studentized_cross_dcov 
            pX = 2 if pX is None else pX
            pY = 2 if pY is None else pY
            stat_params = {'pX':pX, 'pY':pY}
        else:
            raise Exception(f'Experiment type {expt_type} not recognized')

        print(f'Repeitions is {repetitions}')

        Stat = compute_stats(Source,n, stat_func, stat_params, repetitions=repetitions,
                                progress_bar=progress_bar)
        Stat2 = compute_stats(Source2,n, stat_func, stat_params, repetitions=repetitions,
                                progress_bar=progress_bar)
        #### 
        title = f'{kernel_type} Kernel (n={n}, dof={df})' 
        label = f'x-{expt_type}'
        ylabel = 'Density'
        xx = np.linspace(-10, 10, 1000)
        yy = stats.norm.pdf(xx) 
        xx = xx + Stat.mean()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(x=[Stat, Stat2], density=True, label=[f'd={d}', f'd={d_high}'], alpha=0.8)
        ax.plot(xx, yy, label='N(0,1)', color='k') 
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title, fontsize=16)
        ax.legend(fontsize=14, fancybox=True, framealpha=0.0)
        if save_fig:
            figname = f'Second_Moment_NullPdf_{expt_type}_n_{n}_d_{d}_dof_{df}_seed_{seed}'
            figname_tex = f'{fig_dir}/{figname}.tex' 
            figname_png = f'{fig_dir}/{figname}.png' 
            plt.savefig(figname_png, dpi=450)
            tpl.save(figname_tex, axis_width=r'\figwidth', axis_height=r'\figheight')
        else:
            plt.show()
