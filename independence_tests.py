from math import sqrt
from functools import partial
from statistics import median
from crossHSIC import get_HSIC, get_cross_HSIC, get_variance, get_K_L_matrices
from utils import RBFkernel, distance_kernel, median_heuristic, permutation_test 
from utils import CreateDependentGaussianSource, CreateIndependentGaussianSource

from tqdm import tqdm
import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import cdist, pdist 
import matplotlib.pyplot as plt

def crossHSIC_test(XX, YY, kernel_X=None, kernel_Y=None,
                    alpha=0.05):
    """
    Cross-HSIC test based on sample-splitting and 
    studentization. 
    """
    n = len(XX) 
    assert len(YY)==n

    if kernel_X is None:
        bwX = median_heuristic(Z=XX)
        kernelX = partial(RBFkernel, bw=bwX)
    else:
        kernelX = kernel_X
    if kernel_Y is None:
        bwY = median_heuristic(Z=YY)
        kernelY = partial(RBFkernel, bw=bwY)
    else:
        kernelY = kernel_Y
    # compute the kenrel matrices 
    K, L = get_K_L_matrices(XX, YY, kernelX, kernelY)
    # comptue the cross-HSIC statistic 
    cHSIC = get_cross_HSIC(XX=XX, YY=YY, kernel_X=kernelX,
                            kernel_Y=kernelY)
    # get the variance estimate 
    var = get_variance(K, L, cHSIC)
    if var <=0:
        print('Cross-HSIC: The variance estimate should not be negative!!')
        var = 1.0 # set it to some  default value
    # get the statistic value 
    stat = cHSIC*sqrt(n/(2*var))
    # get the rejection threshold 
    th = stats.norm.ppf(1-alpha) 
    return 1.0*(stat>th)



def distance_covaraince_DA_test(XX, YY, pX=2,
                                pY=2, alpha=0.05):
    """
    The dimension-agnostic version of the distance-
    covariance test
    """
    n = len(XX) 
    assert len(YY)==n

    kernelX = partial(distance_kernel, p=pX)
    kernelY = partial(distance_kernel, p=pY)

    # compute the kenrel matrices 
    K, L = get_K_L_matrices(XX, YY, kernelX, kernelY)
    # comptue the cross-HSIC statistic 
    cHSIC = get_cross_HSIC(XX=XX, YY=YY, kernel_X=kernelX,
                            kernel_Y=kernelY)
    # get the variance estimate 
    var = get_variance(K, L, cHSIC)
    if var <=0:
        print('DA Dcov: The variance estimate should not be negative!!')
        var = 1.0 # set it to some  default value
    # get the statistic value 
    stat = cHSIC*sqrt(n/(2*var))
    # get the rejection threshold 
    th = stats.norm.ppf(1-alpha) 
    return 1.0*(stat>th)



def distance_covariance_test(XX, YY, pX=2, pY=2, alpha=0.05,
                                return_stat=False):
    """
    The distance covaraince test for independence 
    proposed by Szekely-Rizzo (2008)
    """
    n = len(XX)
    assert len(YY)==n 

    def helper(Z, p):
        n_ = len(Z)
        DZ = cdist(Z, Z, 'minkowski', p=p) ### \ell_p norm 
        DZrow = np.tile(DZ.mean(axis=1), (n_,1)).transpose()
        DZcol = np.tile(DZ.mean(axis=0), (n_, 1))
        DZmean = DZ.mean()
        AZ = DZ - DZrow - DZcol + DZmean
        return AZ, DZmean
    
    A, Amean = helper(XX, pX)
    B, Bmean = helper(YY, pY)
    # compute the distance covariance statistic 
    V2 = (A*B).mean() 
    # compute the S2 statistic 
    S2 = Amean * Bmean 
    # compute the test statistic
    stat = (V2 * n) / S2  
    if return_stat:
        # return only the statistic
        return V2 
    else:
        ####NOTE: In practice, this approach often seems 
        #### to be very conservative in the n<1000 regime #####
        # compute the rejection threshold 
        th = (stats.norm.ppf((1-alpha/2)))**2
        # reject the null if statistic is larger than threhold
        return 1.0*(stat>th)


def HSIC_permutation_test(XX, YY, kernel_X=None, kernel_Y=None, 
                            alpha=0.05, num_perms=50):
    """
    The standard HSIC test calibrated using permutation 
    distribution
    """
    ### no need to handle the kernel_X=None case here, since that 
    ### is taken care of by the get_HSIC function  
    stat_func = partial(get_HSIC, kernel_X=kernel_X, kernel_Y=kernel_Y)
    # obtaint the result of the test (0 or 1 valued output)
    result = permutation_test( XX=XX, YY=YY, stat_func=stat_func,
                                alpha=alpha, num_perms=num_perms)
    return result 

if __name__=='__main__':
    null_expt=False 
    if null_expt: 
        Source = CreateIndependentGaussianSource()
    else:
        Source = CreateDependentGaussianSource(epsilon=1.0)
    N =500
    num_trials=10
    alpha=0.05
    num_perms=20

    NN = np.linspace(10, N, 15, dtype=int)
    PowerCHSIC = np.zeros(NN.shape)
    PowerHSIC = np.zeros(NN.shape)
    PowerDcov = np.zeros(NN.shape)
    PowerDAdcov = np.zeros(NN.shape)

    # define the statistic computing functions 
    dcov_stat_func = partial(distance_covariance_test, return_stat=True)
    
    
    for trial in tqdm(range(num_trials)):
        for i, n in enumerate(NN):
            if n%2: #only even values work 
                n+=1 
            X, Y = Source(n)
            PowerCHSIC[i] += crossHSIC_test(XX=X, YY=Y, alpha=alpha)
            PowerDcov[i] += distance_covariance_test(XX=X, YY=Y, alpha=alpha)
            # PowerDcov[i] += permutation_test(XX=X, YY=Y, stat_func=dcov_stat_func, 
            #                                  alpha=alpha, num_perms=num_perms)
            PowerDAdcov[i] += distance_covaraince_DA_test(XX=X, YY=Y, alpha=alpha)
            PowerHSIC[i] += HSIC_permutation_test(XX=X, YY=Y,alpha=alpha,
                                                    num_perms=num_perms)
        
    PowerCHSIC /= num_trials 
    PowerDcov /= num_trials 
    PowerDAdcov /= num_trials 

    PowerHSIC /= num_trials
    plt.plot(NN, PowerDcov, label='perm-dCov')
    plt.plot(NN, PowerHSIC, label='perm-HSIC')
    plt.plot(NN, PowerCHSIC, label='cross-HSIC')
    plt.plot(NN, PowerDAdcov, label='cross-dCov')

    plt.xlabel('Sample-size (n)', fontsize=13)
    if null_expt: 
        plt.ylabel('Type-I error', fontsize=13)
        plt.title('Type-I error vs Sample Size', fontsize=15)
        plt.plot(NN, np.ones(NN.shape)*alpha, '--')
    else:
        plt.ylabel('Power', fontsize=13)
        plt.title('Power vs Sample Size', fontsize=15)
    plt.legend()
    if null_expt:
        plt.show()
        # plt.savefig('null.png', dpi=450)
    else:
        plt.show()
        # plt.savefig('alternative.png', dpi=450)

    print(f'Cross-HSIC test: {PowerCHSIC.mean():.2f}')
    print(f'Dcov Test: {PowerDcov.mean():.2f}')
    print(f'DA Dcov Test: {PowerDAdcov.mean():.2f}')