from math import log, sqrt, ceil
import scipy.stats as stats
import numpy as np 
from scipy.spatial.distance import cdist, pdist 
from scipy.stats import t as tdist


def median_heuristic(Z):
    # compute the pairwise distance between the elements of Z 
    dists_ = pdist(Z)
    # obtain the median of these pairwise distances 
    sig = np.median(dists_)
    return sig


def RBFkernel(x, y=None, bw=1.0):
    y = x if y is None else y 
    dists = cdist(x, y, 'euclidean') 
    sq_dists = dists * dists 
    K = np.exp(-sq_dists/(2*bw*bw))
    return K 

def LinearKernel(x, y=None):
    """
        k(x,y) = x^T y  = \sum_{i=1}^d x_i y_i
    """
    y = x if y is None else y 
    K = np.einsum('ji, ki ->jk', x, y) 
    return K 

def PolynomialKernel(x, y=None, p=2, c=1.0):
    """
        k(x,y) = (c + x^Ty)**p 
    """
    L = LinearKernel(x, y) 
    K = (c + L)**p 
    return K 

def RationalQuadraticKernel(x, y=None, c=1.0):
    y = x if y is None else y 
    dists = cdist(x, y, 'euclidean') 
    sq_dists = dists * dists 
    K = 1 - sq_dists/(sq_dists + c)
    return K 

def ExponentialKernel(x, y=None, bw=1.0):
    y = x if y is None else y 
    dists = cdist(x, y, 'euclidean') 
    K = np.exp(-dists/bw)
    return K 



def distance_kernel(XX, YY=None, p=2, zx=None, zy=None, q=1):
    YY = XX if YY is None else YY
    if len(XX.shape)==1:
        XX = XX.reshape((-1,1))
        YY = YY.reshape((-1,1))
    def helper(X_, z):
        if z is None:
            Z = np.zeros(X_.shape)
        elif isinstance(z, np.ndarray):
            assert len(z)==X_.shape[1] 
            if len(z.shape)>1:
                z = z.reshape((-1,))
            Z = np.tile(z, (X_.shape[0], 1))
        # compute the distance kernel 
        # 2k(x,y) = d(x, z0) +  d(y, z0) - d(x, y)
        DXX = cdist(X_, X_,'minkowski', p=p)**q
        DXZ = cdist(X_, Z, 'minkowski',p=p)**q
        DZX = cdist(Z, X_,'minkowski', p=p)**q
        K_ = 0.5*(DXZ + DZX - DXX)
        return K_ 
    # compute the X-kernel matrix 
    KX = helper(X_=XX, z=zx)
    # compute the Y-kernel matrix 
    KY = helper(X_=YY, z=zy)
    # product kernel matrix 
    K = KX * KY 
    return K 


def CreateIndependentTSource(df1=1, loc1=0, scale1=1.0, ndims1=1,
                                df2=1, loc2=0, scale2=1.0, ndims2=1):    
    def Source(n):
        X_ = tdist.rvs(df=df1, loc=loc1, scale=scale1, size=int(n*ndims1))
        X = X_.reshape((n, ndims1))
        Y_ = tdist.rvs(df=df2, loc=loc2, scale=scale2, size=int(n*ndims2))
        Y = Y_.reshape((n, ndims2))
        return X, Y
    return Source




def CreateIndependentGaussianSource(mean0=None, mean1=None, cov0=None, cov1=None, 
                                    ndims=10):
    # set the default values of the mean vectors 
    if mean0 is None:
        mean0 = np.zeros((ndims,))
    if mean1 is None:
        mean1 = np.ones((ndims,))
    d0, d1 = len(mean0), len(mean1)
    # set the default values of the covariance matrices
    cov0 = np.eye(d0) if cov0 is None else cov0
    cov1 = np.eye(d1) if cov1 is None else cov1
    # define the source function 
    def Source(n):
        X = np.random.multivariate_normal(mean=mean0, cov=cov0, size=n)
        Y = np.random.multivariate_normal(mean=mean1, cov=cov1, size=n)
        return X, Y
    # return the source function 
    return Source 

def CreateDependentGaussianSource(mean=None, cov=None, func=None, 
                                    ndims=10, epsilon=0.1, p=2): 
    """
    returns a Source function that generates X, Y 
    where X is multivariate normal and 
    Y = epsilon * (X**p) + (1-epsilon) * (Y**p)
    """
    # set the dfault value of the mean vector
    if mean is None: 
        mean = np.zeros((ndims,))
    # set the default vlaue of the covariance matrix 
    if cov is None:
        cov = np.eye(len(mean))
    # set the default vlaue of the function relating the two samples 
    if func is None: 
        func = lambda x, z: epsilon*(x**p)  + (1-epsilon)*(z**p)
    # define the source function 
    def Source(n): 
        # generate the first sample 
        X = np.random.multivariate_normal(mean=mean, cov=cov, size=n) 
        Z = np.random.multivariate_normal(mean=mean, cov=cov, size=n)
        Y = func(X, Z) 
        return X, Y
    # return the source function 
    return Source 

def permutation_test(XX, YY, stat_func=None, alpha=0.05, 
                        num_perms=200, return_thresh=False): 
    if stat_func is None:
        raise Exception("No stat_func passed!!!")
    n = len(XX)
    assert n==len(YY)
    vals = np.zeros((num_perms,)) 
    
    for i in range(num_perms):
        permX = np.random.permutation(n) 
        permY = np.random.permutation(n) 
        XX_, YY_ = XX[permX], YY[permY] 
        vals[i] = stat_func(XX_, YY_) 

    vals = np.sort(vals) 
    th_idx = min(max(0, ceil((1-alpha)*num_perms)-1), num_perms-1)
    thresh = vals[th_idx] 
    if return_thresh:
        return thresh 
    else:
        stat = stat_func(XX, YY) 
        return 1.0*(stat>thresh)


def compute_bootstrap_std(XX, YY, stat_func, stat_kwargs, B=100):
    """
        Compute the bootstrap estimate of the standard deviation of 
        a statistic computed by stat_func on data (XX, YY)
    """
    if stat_kwargs is None:
        stat_kwargs = {}
    Stat = np.zeros((B,))
    n = len(XX)
    assert n==len(YY)
    for trial in range(B): 
        # get the sample with replacement 
        idx = np.random.choice(a=n, size=n)
        XX_, YY_ = XX[idx], YY[idx]
        Stat[trial] = stat_func(XX_, YY_, **stat_kwargs)
    # compute the bootstrap variance 
    m = Stat.mean() 
    std = sqrt( np.sum((Stat - m)**2) / B )
    return std 


def initialize_kernel(kernel_type='RBF+median'):
    if kernel_type=='RBF':
        kernel_X, kernel_Y = RBFkernel, RBFkernel
    elif kernel_type=='RBF+median':
        # when kernels are set to None, all 
        # functions for computing statistics 
        # use the RBF kernel with bandwidth chosen 
        # via the median heuristic.
        # so we don't need to do anything here
        kernel_X, kernel_Y= None, None
    elif kernel_type=='Linear':
        kernel_X = LinearKernel
        kernel_Y = LinearKernel
    elif kernel_type=='Polynomial':
        kernel_X = PolynomialKernel
        kernel_Y = PolynomialKernel
    elif kernel_type=='RationalQuadratic':
        kernel_X = RationalQuadraticKernel
        kernel_Y = RationalQuadraticKernel
    elif kernel_type=='Exponential':
        kernel_X = ExponentialKernel
        kernel_Y = ExponentialKernel
    else:
        # print('Warning: kernel type {} not recognized')
        # print('Falling back to "RBF+median" kernel')
        kernel_X = None 
        kernel_Y = None 
    return kernel_X, kernel_Y


def generate_seeds(n=10, start=1234):
    if start==0:
        start = 1
    func = lambda s: (s*2022)%10000 
    seeds = [start]
    for i in range(n-1):
        seed = func(seeds[i]) 
        seeds.append(seed)
    for s in seeds: 
        print(s)



def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def predict_power(Power, alpha=0.05, factor=2):
    """
        Predict the power of MMD test using the power of cMMD test
        \Phi( z_{\alpha} + \sqrt{2}( \Phi^{-1}(Power) - z_{\alpha}))
    """
    normal = stats.norm
    z_a = normal.ppf(alpha) 
    sqrt_factor = sqrt(factor)
    temp = z_a + sqrt_factor*np.array([
        (normal.ppf(p) - z_a) for p in Power
    ])
    power = normal.cdf(temp) 
    return power



if __name__=='__main__':
    generate_seeds(n=10)    