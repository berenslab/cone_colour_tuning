import numpy as np

# for bootstrapping
from sklearn.utils import resample

import statsmodels
import scipy as scp



"""
helper functions
"""

# define normalizing function
def normalize(data, mode=None):
    """
    normalizes data with respect to mode
    :data: shape: (n,tpts)
    :param mode: in ['z_score', 'std_norm', 'max_division']
    
    """
    assert mode in ['z_score', 'std_norm', 'max_division']
    if mode=='z_score':
        return (data - data.mean(axis=1)[:, None])/(data.std(axis=1)[:,None])
    elif mode=='std_norm':
        return data /(data.std(axis=1)[:,None])
    elif mode=='max_division':
        return data /(np.abs(data).max(axis=1)[:,None])
    
    
def get_bootstrapped_confs(values, p_conf=95,data_per_iteration=0.75,n_iterations=1000, verbose=False):
    """
    returns p_conf-confidence intervals for mean and data mean
    ---
    values: 1d array
    p_conf: confidence level in [0,100]
    data_per_iteration: fraction of data to use in each iteration in [0,1]
    returns: percentiles and data mean
    """
    
    n_size = int(len(values) * data_per_iteration)
    if verbose:
        print('size of resampled data for bootstrap:', n_size,'out of', len(values))
    
    bt_means = np.zeros(n_iterations)

    # run bootstrap
    for i in range(n_iterations):
        # get one batch
        bt_data = resample(values, n_samples=n_size, random_state=i)
        bt_means[i] = np.mean(bt_data) 
    percentiles = [np.percentile(bt_means,(100-p_conf)/2) ,
                   np.mean(values),
                   np.percentile(bt_means, 100-(100-p_conf)/2) 
                    ]
    return percentiles
    

    
    
"""
this implementation for binomial conf intervals is inspired by:
https://stackoverflow.com/questions/13059011/is-there-any-python-function-library-for-calculate-binomial-confidence-intervals
"""

def cdf_binP(N,p,x1,x2):
    """
    cdf over [x1, x2] (including bounds)
    """
    # stats.binom.cdf includes the upper and lower bound
    return scp.stats.binom.cdf(x2,N,p) - scp.stats.binom.cdf(x1-1,N,p)

def calcBin_conf_interval(k, N, cl = 95):
    '''
    Calculate the exact confidence interval for a binomial proportion
    k : succesful trials
    N : total number of throws
    cl: confidence level
    ---

    ''' 
    precision = 10**-5
    
    k = float(k)
    N = float(N)
    #Set the confidence bounds (here symmetric)
    vTU = (100 - float(cl))/2 #upper 
    vTL = vTU # lower

    p_hat = k/N #estimator for p 
    
    # calculate lower bound
    if(k==0):
            lower_bound = 0.0
    else:
        v = p_hat/2 #initial guess 
        v_lower = 0
        v_upper = p_hat
        q = vTL/100

        # iterate over precision for lower_bound
        while((v_upper-v_lower) > precision):
            # original:
            if(cdf_binP(N, v, k+1, N) > q):
                    v_upper = v
                    v = (v_lower+v)/2
            else:
                    v_lower = v
                    v = (v+v_upper)/2
      
        lower_bound = v
    
    # calculate upper bound
    if(k==N):
            upper_bound = 1.0
    else:
        v = (1+p_hat)/2 # initial guess
        v_lower =p_hat
        v_upper = 1
        q = vTU/100
        while((v_upper-v_lower) > precision):
            if(cdf_binP(N, v, 0, k-1) < q):
                    v_upper = v
                    v = (v_lower+v)/2
            else:
                    v_lower = v
                    v = (v+v_upper)/2
        upper_bound = v
    return (lower_bound, upper_bound)