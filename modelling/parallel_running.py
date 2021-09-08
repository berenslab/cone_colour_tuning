from hc_model import run_model_hc
from hc_model import extract_params
from hc_model import get_param_labels

import numpy as np

def run_model_x(params, o, mode, N=200, return_H_in=False):
    """
    helper to extract parameters
    and set params 0 if negative
    """
    
    w,dc,a = extract_params(params,mode)
    
    # set to 0 if negative
    w[w<0] = 0 
    dc[dc<0] = 0 
    a[a<0] = 0 
    
    if return_H_in:
        k_fit, H_in = run_model_hc(o,w,dc,a,N=N,
                             return_H_in=return_H_in)
        return k_fit, H_in #, baseline_store
        
    else:
        k_fit = run_model_hc(o,w,dc,a,N=N,
                             return_H_in=return_H_in)
        return k_fit #, baseline_store



def mse(fit, data):
    """
    compute sqrt(mse)!!!
    """
    #mse = 1/(np.shape(data)[0]*np.shape(data)[1]) * np.sum((data - fit)**2)
    mse = 1/(len(data)) * np.sum((data - fit)**2)
    #if mse >20:
    #    mse = 100
    if np.isnan(mse):
        #raise Warning('mse is Nan. Check your model.')
        return 10e3
    else:
        return min(np.sqrt(mse), 10e3)

def runparallelsamples(sample,wave_blocked, mode, wave_normal,return_H_in=False):
    """
    returns model output and loss for one given parameterset
    """
    if return_H_in:
        model_out, H_in =  run_model_x(sample, wave_blocked, mode=mode, return_H_in=return_H_in)
        loss = mse(model_out.flatten(), wave_normal.flatten())
        return model_out, loss, H_in
    
    else:
        model_out =  run_model_x(sample, wave_blocked, mode=mode, return_H_in=return_H_in)
        loss = mse(model_out.flatten(), wave_normal.flatten())
        return model_out, loss