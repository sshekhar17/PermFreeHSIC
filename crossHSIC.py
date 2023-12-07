"""
Implementation of the cross-HSIC statistic 
"""

from time import time 
from functools import partial
from random import random
import numpy as np 
import scipy.stats as stats 
from tqdm import tqdm 
import matplotlib.pyplot as plt 

from utils import *
from math import dist, sqrt 


# Compute the cross HSIC statistic 
def get_K_L_matrices(XX, YY, kernel_X=None, kernel_Y=None):
    N = len(XX) 
    n = N//2 
    if (N!=int(2*n)): # i.e. N must be even 
        N -= 1 
        XX = XX[:-1]
        YY = YY[:-1]

    if kernel_X is None:
        kernel_X = RBFkernel
    if kernel_Y is None:
        kernel_X = RBFkernel
    
    K = kernel_X(XX) 
    K[:n, :n] = 0
    K[n:, n:] = 0

    L = kernel_Y(YY) 
    L[:n, :n]=0
    L[n:, n:]=0

    return K, L


def get_T1(K, L, n=None): 
    if n is None:
        N = len(K) 
        assert N%2==0 
        n = N//2 
    # the pointwise product 
    M = K*L 
    # sum the last n columns of the 2D array M 
    Mu = np.sum(M[:, n:], axis=1) 
    # sum the first n terms of the 1D array Mu
    T1 = (1/(n*n))*np.sum(Mu[:n]) 
    return T1
    
def get_T2(K, L, n=None): 
    if n is None:
        N = len(K) 
        assert N%2==0 
        n = N//2 
    # multiply the two kernel matrices 
    KL = np.matmul(K, L) 
    # sum the last n columns of KL 
    term10 = np.sum(KL[:, n:], axis=1)
    # sum the last n elements of the 1D array term10 
    term1 = np.sum(term10[n:]) 
    # get the second term 
    term2 = 0.5*np.trace(KL)
    # get T2 
    assert n>1
    T2 = (1/(n*n*(n-1)))*(term1 - term2)
    return T2

def get_T3(K, L, n=None): 
    if n is None:
        N = len(K) 
        assert N%2==0 
        n = N//2 
    # multiply the two kernel matrices 
    KL = np.matmul(K, L) 
    # sum the first n columns of KL 
    term10 = np.sum(KL[:, :n], axis=1)
    # sum the first n elements of the 1D array term10 
    term1 = np.sum(term10[:n]) 
    # get the second term 
    term2 = 0.5*np.trace(KL)
    # get T3 
    assert n>1
    T3 = (1/(n*n*(n-1)))*(term1 - term2)
    return T3
    


def get_T4(K, L, n=None):
    if n is None:
        N = len(K) 
        assert N%2==0 
        n = N//2 
    # get the KL matrix 
    KL = np.matmul(K, L)

    # get the first term 
    term1 = (np.sum(np.sum(K[:, :n], axis=1)[n:]) *
                np.sum(np.sum(L[:, :n], axis=1)[n:]) )
    # get the second term 
    term2 = np.sum(np.sum(KL[:, n:], axis=1)[n:])
    # get the third term 
    term3 = np.sum(np.sum(KL[:, :n], axis=1)[:n])
    # get the fourth term 
    term4 = 0.5*np.trace(KL)
    # combine the four terms to get T4
    assert n>1
    den = (n*(n-1))**2 # denominator 
    T4 = (1/den)*(term1 - term2 -term3 + term4) 
    return T4

def get_cross_HSIC(XX, YY, kernel_X=None, kernel_Y=None):
    if kernel_X is None:
        kernel_X = RBFkernel
    if kernel_Y is None:
        kernel_Y = RBFkernel
    # some preprocessing 
    N = len(XX)
    if N%2==1:
        XX = XX[:-1]
        YY = YY[:-1] 
        N = N-1 
        assert N>=2 
        # require both XX and YY to have the same number of observations 
        assert(len(YY)==N) 
    # get the kernel matrices 
    K, L = get_K_L_matrices(XX, YY, kernel_X, kernel_Y) 
    # get the four terms 
    T1 = get_T1(K, L)
    T2 = get_T2(K, L)
    T3 = get_T3(K, L)
    T4 = get_T4(K, L) 
    # obtain the cross HSIC statistic 
    cHSIC = T1 - T2 - T3 + T4
    return cHSIC 


def get_w_tilde(K, L):
    N = len(K) 
    n = N//2 
    assert 2*n==N 
    assert len(L)==N
    assert n>2
    # get the KL matrix 
    KL = np.matmul(K, L)
    KL1 = np.sum(KL, axis=1)

    LK = np.matmul(L, K)
    LK1 = np.sum(LK, axis=1)

    K1l = np.sum(K[:, n:], axis=1)
    L1l = np.sum(L[:, n:], axis=1)
    KL1l = np.sum(KL[:, n:], axis=1)

    # w_tilde consists of 6 terms 
    term1 = (n/(2*(n-1))) * np.sum(K*L, axis=1) 
    term2 = (1/(4*(n-1))) * np.trace(KL) * np.ones((int(2*n),)) 
    term3 = (1/(2*(n-1))) * (KL1 + LK1)
    term4 = 1/(2*(n-1)) * K1l * L1l
    term5 = 1/(2*n*(n-1)) * np.sum(KL1l[n:])*np.ones((int(2*n),)) 
    term6 = 1/(4*n*(n-1)) * ( L.sum()*K1l  +  K.sum() * L1l  )

    # get w_tilde 
    w_tilde = term1 + term2 - term3 - term4 - term5 + term6 

    return w_tilde 


def get_variance(K, L, cHSIC):
    N = len(K) 
    n = N//2 
    assert 2*n==N 
    assert len(L)==N
    assert n>2   

    # get w_tilde
    w_tilde = get_w_tilde(K, L) 
    # obtain w by selecing the first n terms of w_tilde
    w = w_tilde[:n]
    # obtain the variance 
    term0 = 4*(n-1)/((n-2)*(n-2)) 
    term1 = (term0/((n-1)*(n-1))) * np.sum(w*w) 
    term2 = term0 * n * cHSIC * cHSIC 

    variance = term1 - term2 

    return variance 

def get_HSIC(XX, YY, kernel_X=None, kernel_Y=None, num_pts_bw=20):
    if kernel_X is None:
        bwX = median_heuristic(Z=XX[:num_pts_bw])
        kernel_X = partial(RBFkernel, bw=bwX)
    if kernel_Y is None:
        bwY = median_heuristic(Z=YY[:num_pts_bw])
        kernel_Y = partial(RBFkernel, bw=bwY)
    # some preprocessing 
    n = len(XX)
    # get the kernel matrices 
    K, L = kernel_X(XX, XX), kernel_Y(YY, YY)
    # define the H matrix 
    H = np.eye(n) - (1/n)*np.ones((n,n)) 
    # compute the empirical HSIC estimate 
    Kc = np.dot(np.dot(H, K), H)
    Lc = np.dot(np.dot(H, L), H)
    hsic_stat = np.sum(Kc.T * Lc) / (n*n)
    return hsic_stat


def get_studentized_cross_hsic(XX, YY, kernel_X=None,
                                kernel_Y=None, num_pts_bw=20):
    if kernel_X is None:
        bwX = median_heuristic(Z=XX[:num_pts_bw])
        kernel_X = partial(RBFkernel, bw=bwX)
    if kernel_Y is None:
        bwY = median_heuristic(Z=YY[:num_pts_bw])
        kernel_Y = partial(RBFkernel, bw=bwY)
    n = len(XX)
    # compute the kenrel matrices 
    K, L = get_K_L_matrices(XX, YY, kernel_X, kernel_Y)
    # comptue the cross-HSIC statistic 
    cHSIC = get_cross_HSIC(XX=XX, YY=YY, kernel_X=kernel_X,
                            kernel_Y=kernel_Y)
    # get the variance estimate 
    var = get_variance(K, L, cHSIC)
    if var <=0:
        print('Cross-HSIC: The variance estimate should not be negative!!')
        var = 1.0 # set it to some  default value
    # get the statistic value 
    stat = cHSIC*sqrt(n/(2*var))
    return stat

def get_studentized_cross_dcov(XX, YY, pX=2, pY=2):
    kernelX = partial(distance_kernel, p=pX)
    kernelY = partial(distance_kernel, p=pY)
    # comptue the studentized hsic with distance kernels 
    stat = get_studentized_cross_hsic(XX=XX, YY=YY, 
                                    kernel_X=kernelX, 
                                    kernel_Y=kernelY)
    return stat


#=====================
# QUARTIC-TIME IMPLEMENTATION 

def getH(XX, YY, kernelX, kernelY):
    # some preprocessing 
    N = len(XX) 
    n = N//2 
    assert n==(N-n) # N is even 
    assert N==len(YY) 
    assert n>0
    # partition the dataset 
    X1, X2 = XX[:n], XX[n:]
    Y1, Y2 = YY[:n], YY[n:]
    # now get the kenrel matrices 
    K = kernelX(X1, X2) 
    L = kernelY(Y1, Y2)

    # defin the inner product 
    # \langle f, \phi(X_i)\psi(Y_j) \rangle 
    def inner_prod(i, j):
         # the term \langle f, \phi(X_i) \psi(Y_i) \rangle
        term0 = (K[i]*L[j]).sum()
        term1 = term0/n
        term2 = ( (K[i].sum()) * (L[j].sum()) - term0 ) / (n*(n-1))
        return  term1 - term2 
    # define the h-function 
    def hfunc(i, j):
        if i==j:
            return 0 
               # the term \langle f, \phi(X_i) \psi(Y_i) \rangle
        term_ii = inner_prod(i, i)
        term_jj = inner_prod(j, j)
        term_ij = inner_prod(i, j)
        term_ji = inner_prod(j, i)
        h = 0.5*(term_ii + term_jj - term_ij -term_ji)
        return h 
    # now compute the H matrix with 
    # H[i, j] = h(Z_i, Z_j) from the note
    H = np.zeros((n, n))
    #TODO: vectorize this 
    for i in range(n):
        for j in range(n):
            H[i][j] = hfunc(i, j) 
    return H

def get_quartic_studentized_cHSIC(XX, YY, kernelX=None, kernelY=None, num_pts_bw):
    if kernelX is None:
        bwX = median_heuristic(XX[:num_pts_bw])
        kernelX = partial(RBFkernel, bw=bwX)
    if kernelY is None:
        bwY = median_heuristic(YY[:num_pts_bw])
        kernelY = partial(RBFkernel, bw=bwY)
    # get the H-matrix 
    H = getH(XX, YY, kernelX, kernelY)
    # sanity check 
    n, n_ = H.shape 
    # H must be a square matrix 
    assert n==n_ 
    # needed for variance computation 
    assert n>=3  
    # compute the cross-HSIC statistic (numerator) from H 
    cHSIC = (H.sum() - np.trace(H))/(n*(n-1))
    # compute the empirical standard deviation 
    H_row_sum = H.sum(axis=1)
    H_row_sum_ = H_row_sum - np.diag(H) 
    var = ( ( H_row_sum_/(n-1) - cHSIC ) ** 2 ).sum() 
    var = 4*(n-1)/((n-2)**2) * var 
    std = sqrt(var)
    stat = sqrt(n)*cHSIC/std 
    return stat, cHSIC, var

    







    
    
#=====================

def PlotNullDistribution(n=1000, num_trials=50, d=5, independent=True, 
                        random_seed=None, progress_bar=True, kernelX=None,
                        kernelY=None, savefig=False, figname=None, 
                        return_stats=False, plot_quartic=False, num_pts_bw=20):
    # set the random seed if needed 
    if random_seed is not None:
        np.random.seed(random_seed)
    # initialize the source for the observations 
    mean0, mean1 = np.zeros((d,)), np.ones((d,))
    if independent: 
        Source = CreateIndependentGaussianSource(mean0=mean0, mean1=mean1) 
    else:
        Source = CreateDependentGaussianSource(mean=mean0) 
    # initialize the variables 
    cHSIC_vals = np.zeros((num_trials,)) # cross-HSIC 
    var_vals = np.zeros((num_trials,)) # the variance terms 
    stat_vals = np.zeros((num_trials,)) # the studentized cross-HSIC stats

    if plot_quartic:
        cHSIC_vals2 = np.zeros((num_trials,)) # cross-HSIC 
        var_vals2 = np.zeros((num_trials,)) # the variance terms 
        stat_vals2 = np.zeros((num_trials,)) # the studentized cross-HSIC stats


    # the main loop
    range_ = tqdm(range(num_trials)) if progress_bar else range(num_trials)
    for i in range_:
        X, Y = Source(n) 
        # initialize the kernels if needed
        if kernelX is None:
            bwX = median_heuristic(X[:num_pts_bw])
            kernelX_ = partial(RBFkernel, bw=bwX)
        else:
            kernelX_=kernelX 
        if kernelY is None:
            bwY = median_heuristic(Y[:num_pts_bw])
            kernelY_ = partial(RBFkernel, bw=bwY)
        else:
            kernelY_=kernelY 
        # compute the kenrel matrices 
        K, L = get_K_L_matrices(X, Y, kernelX_, kernelY_)
        # comptue the cross-HSIC statistic 
        cHSIC = get_cross_HSIC(XX=X, YY=Y, kernel_X=kernelX_,
                                kernel_Y=kernelY_)
        # get the variance estimate 
        var = get_variance(K, L, cHSIC)
        if var <=0:
            print('The variance estimate should not be negative!!')
            var = 1.0 # set it to some  default value
        
        cHSIC_vals[i] = cHSIC 
        var_vals[i] = var 
        stat_vals[i] = cHSIC*sqrt(n/(2*var))

        if plot_quartic:
            s, c, v = get_quartic_studentized_cHSIC(X, Y, kernelX=kernelX_, 
                        kernelY=kernelY_)
            cHSIC_vals2[i] = c 
            var_vals2[i] = v 
            stat_vals2[i] = s

    if plot_quartic:
        errS = np.linalg.norm(stat_vals2-stat_vals, ord=np.inf)
        errC = np.linalg.norm(cHSIC_vals2-cHSIC_vals, ord=np.inf)
        errV = np.linalg.norm(var_vals2-var_vals, ord=np.inf)
        print(f'Deviation b/w stats:{errS}')
        print(f'Deviation b/w cHSIC:{errC}')
        print(f'Deviation b/w var:{errV}')
           

    # plot the figure 
    xx = np.linspace(-10, 10, 1000)
    pp = stats.norm.pdf(xx)
    if not independent: 
        xx = xx + stat_vals.mean()
    
    fig, ax = plt.subplots() 
    ax.hist(stat_vals, density=True, bins=25, label='c$\widehat{HSIC}$ ($n^2$)', alpha=0.4)
    if plot_quartic:
        ax.hist(stat_vals2, density=True, bins=25, label='c$\widehat{HSIC}$ ($n^4$)', alpha=0.4)
    ax.plot(xx, pp, label='N(0,1)', color='k') 
    ax.set_xlabel('Value of studentized c$\widehat{HSIC}$ statistic', 
                    fontsize=13)
    ax.set_ylabel('density', fontsize=13)
    ax.set_title('Distribution of studentized c$\widehat{HSIC}$ statistic', 
                    fontsize=15)
    ax.legend(fontsize=13)

    if savefig:
        if figname is None:
            figname = 'temp.png'
        plt.savefig(figname, dpi=450) 
    if return_stats:
        if plot_quartic:
            return stat_vals, cHSIC_vals, var_vals, stat_vals2, cHSIC_vals2, var_vals2
        else:
            return stat_vals, cHSIC_vals, var_vals

def main_func(seed=1, n=200, d=10, num_trials=500):
    kernelX = None
    kernelY = None
    save_fig=False
    plot_quartic = False

    if plot_quartic:
        stat_vals, cHSIC_vals, var_vals, stat_vals2, cHSIC_vals2, var_vals2 = PlotNullDistribution(n=n, num_trials=num_trials,
                                                                d=d, independent=True,
                                                                kernelX=kernelX,
                                                                kernelY=kernelY,
                                                                return_stats=True, 
                                                                random_seed=seed, 
                                                                plot_quartic=True)
    else:
        stat_vals, cHSIC_vals, var_vals = PlotNullDistribution(n=n, num_trials=num_trials,
                                                                d=d, independent=True,
                                                                kernelX=kernelX,
                                                                kernelY=kernelY,
                                                                return_stats=True, 
                                                                random_seed=seed, 
                                                                plot_quartic=False)

    # Compute the KS distance between the statistic values and standard normal distribution
    norm_cdf = stats.norm.cdf 
    ks_dist, pval = stats.kstest(stat_vals, norm_cdf)
    print(f'\n KS test stat val={ks_dist:.4f}; and p-value {pval:.3f}')
    # # plot the figure 
    # xx = np.linspace(-10, 10, 1000)
    # pp = stats.norm.pdf(xx)

    # fig, ax = plt.subplots() 
    # ax.hist(stat_vals, density=True, bins=20, label='c$\widehat{HSIC}$ ($n^2$)', alpha=0.4)
    # if plot_quartic:
    #     ax.hist(stat_vals2, density=True, bins=10, label='c$\widehat{HSIC}$ ($n^4$)', alpha=0.4)
    # ax.plot(xx, pp, label='N(0,1)', color='k') 
    # ax.set_xlabel('Value of studentized c$\widehat{HSIC}$ statistic', 
    #                 fontsize=13)
    # ax.set_ylabel('density', fontsize=13)
    # ax.set_title('Distribution of studentized c$\widehat{HSIC}$ statistic', 
    #                 fontsize=15)
    # ax.legend(fontsize=13)
    # if save_fig:
    #     plt.savefig('quartic.png', dpi=450)
    # else:
    #     plt.show()

if __name__=='__main__':
    seed = int(time())%10000
    n, d = 200, 200
    main_func(seed=seed, n=n, d=d)