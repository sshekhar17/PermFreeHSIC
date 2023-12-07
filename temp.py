from math import sqrt 
from tqdm import tqdm 
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.stats as stats
import tikzplotlib as tpl 

from utils import CreateIndependentGaussianSource, CreateDependentGaussianSource
from independence_tests import distance_covariance_test


if False:

    # Compare the distribution of the distance-covariance statistic under the null 
    n = 200
    num_trials = 5000
    Source = CreateIndependentGaussianSource(ndims=200)
    # Source = CreateDependentGaussianSource(ndims=100, epsilon=1.0)
    alpha = 0.05
    idx_ = int(num_trials*(1-alpha))

    Stats = np.zeros((num_trials, )) 
    for i in tqdm(range(num_trials)): 
        XX, YY = Source(n)
        # XX = (XX>0)*1.0
        # YY = (YY>0)*1.0
        Stats[i] = distance_covariance_test(XX=XX, YY=YY, return_stat=True)

    Stats_ = np.sort(Stats)
    th_ = Stats_[idx_]
    th = (stats.norm.ppf((1-alpha/2)))**2

    plt.figure()
    plt.hist(x=Stats, density=True, bins=20, alpha=0.8)
    plt.axvline(x=th_, color='r')
    plt.title("Distribution of distance-covariance statistic")
    figname = './figures/dcov_example2.tex'
    tpl.save(figname, axis_width=r'\figwidth', axis_height=r'\figheight')

    # plt.show()
    # plt.axvline(x=th, color='k')
    # plt.xlim([0.5, 4])


if True: 
    ###=======================================================================
    #### epsilon = 0.3
    NN = np.array([
        10, 72, 135, 197, 260, 323, 385, 448, 511, 573, 636, 698, 761, 824,
        886, 949, 1012, 1074, 1137, 1200
    ])
    power_da = np.array([
        0.056, 0.076, 0.18, 0.296, 0.42, 0.616, 0.712, 0.756, 0.86, 
        0.92, 0.976, 0.988, 0.984, 1.0, 0.996, 0.996
        ])
    power_da = np.concatenate((power_da, np.ones((20-len(power_da),))))
    power_perm = np.array([
        0.076, 0.224, 0.448, 0.588, 0.736, 0.9, 0.948, 0.984, 0.992
    ])


    ###=======================================================================

    ###=======================================================================
    #### eps = 0.4
    # NN = np.array([
    #     10, 51, 93, 134, 176, 217, 259, 301, 342, 384, 425, 
    #     467, 508, 550, 592, 633, 675, 716, 758, 800
    # ])
    # power_da = np.array([
    #     0.048, 0.116, 0.32, 0.356, 0.564, 0.764, 0.876, 0.928, 
    #     0.968, 0.988, 0.992, 0.996
    #     ])
    # power_da = np.concatenate((power_da, np.ones((20-len(power_da),))))

    # power_perm = np.array([
    #     0.152, 0.464, 0.648, 0.788, 0.944, 0.972

    # ])
    ###=======================================================================
    power_perm = np.concatenate((power_perm, np.ones((20-len(power_perm),))))


    power_pred = predict_power(power_da)

    # for p, pp in zip(power_da, power_pred): 
    #     print(f'{pp/p:.2f}')
    
    n = 15
    plt.plot(NN[:n], power_da[:n], label='cross-HSIC')
    plt.plot(NN[:n], power_perm[:n], label='HSIC-perm')
    plt.plot(NN[:n], power_pred[:n], '--', label='predicted')
    plt.title(r"$\epsilon$ = $0.3$")
    plt.legend()
    plt.savefig("predicted_power_eps_0_3.png")
    # plt.show()