import numpy as np

"""
functions 
"""


def get_H_in(w, waves):
    """
    w: weights
    waves: input waves from cones
    ---
    input to HC_i(lambda)
    """
    H_in = np.dot(w,waves)
    return H_in


def get_hc_out_normfrac(w):
    """
    calculate the normalizing fraction term for HC output
    """
    '''
    sum_1 = np.sum(w, axis=1)
    hc_out_frac = (w.T/sum_1).T
    '''
    hc_out_frac=np.zeros(np.shape(w))
    sum_1 = np.sum(w, axis=1)
    for i in range(3):
        if sum_1[i]==0:
            hc_out_frac[i] = np.zeros(4)
        else:
            hc_out_frac[i] = w[i]/sum_1[i]
    return hc_out_frac


def mult_mat_vec(M,v):
    """
    multiplying Matrix column_i with v_i
    """
    M_out = np.zeros(np.shape(M))
    for i in range(len(v)):
        M_out[i] = M[i] * v[i]
        
    return M_out
        
    
def get_baseline(w,dc,a, N=100, return_H_in=False):
    """
    iteration to get baseline (of k's) [line 258 ff]
    ---
    w: connectivity matrix
    dc: dark currents
    a: synaptical feedback strength HC->cone 
    ---
    N: number of iterations
    """
    # initialization of baseline:
    baseline = np.copy(dc) # c synaptic (output strength per cone)
    darkC = np.copy(dc) # c synaptic (output strength per cone)

    # calculate normalizing fraction for output of HC 
    hc_out_normfrac = get_hc_out_normfrac(w)

    baseline_store = np.zeros((N,4))
    H_input_store = np.zeros((N,3))

    for i in range(N):

        # HC input is input_weights*baseline
        H_input = np.dot(w,baseline)
        H_input_store[i] = H_input
       
        if i ==0:
            H_input_store[-1:-3] = H_input


        current = darkC - (0.2 + min( (i/N)*0.8*1,0.8 ))  * (np.dot(hc_out_normfrac.T, H_input*a) 
                                 + np.dot(hc_out_normfrac.T, H_input_store[i-1]*a)
                                 + np.dot(hc_out_normfrac.T, H_input_store[i-2]*a)) /3
        

        # update baseline and correct for negative values
        baseline = np.copy(current) 
        baseline[baseline<0]=0
        baseline_store[i] = baseline 
    
    if return_H_in:
        return baseline, baseline_store, H_input_store[-1]
    else:
        return baseline, baseline_store#, H_input_store



def get_output_spectra(o, w, dc, a, N=100, return_H_in=False):
    """
    iteration to get output spectra k [line 304 ff]
    """
    
    curve_len = np.shape(o)[1]

    # initialize cone output with recorded HC blocking data
    current = np.copy(o)

    # calculate normalizing fraction for output of HC 
    hc_out_normfrac = get_hc_out_normfrac(w)

    # raw cone output (opsin*dark_current)
    cone_out_raw = mult_mat_vec(o,dc)

    current_store = np.zeros((N,4,len(o[0])))
    H_in_store = np.zeros((N,3,curve_len))
    
    for i in range(N):
        H_in = np.dot(w,current)
        H_in_store[i] = H_in
        if i ==0:
            H_in_store[-1:-5] = H_in
        
        current = cone_out_raw - (0.2 + min( (i/N)*0.8*1,0.8 )) *  ( np.dot(hc_out_normfrac.T, mult_mat_vec(H_in,a))
                                           + np.dot(hc_out_normfrac.T, mult_mat_vec(H_in_store[i-1],a))
                                           + np.dot(hc_out_normfrac.T, mult_mat_vec(H_in_store[i-2],a)))/3
        current_store[i] = current
        
    if return_H_in:
        return  current, current_store, H_in_store
    
    else:
        return current, current_store
    


"""
put everything together 
"""

def run_model_hc(o,w,dc,a, N=100, return_H_in=False):
    """
    run baseline and output_spectra, and rescales to baseline
    ---
    o: opsin curves
    w: connectivity pattern shape (3,4)
    dc: dark currents, shape=(4)
    a: HC output strenth, shape=(3)    
    """
    
    curve_len = np.shape(o)[1]
    
    if return_H_in:
        baseline, baseline_store, H_in_base = get_baseline(w,dc,a, N=N, return_H_in=return_H_in)
        k_fit, k_fit_store, H_in_store = get_output_spectra(o, w,dc,a, N=N, return_H_in=return_H_in)
    
    else:    
        baseline, baseline_store = get_baseline(w,dc,a, N=N, return_H_in=return_H_in)
        k_fit, k_fit_store = get_output_spectra(o, w,dc,a, N=N, return_H_in=return_H_in)
    
    # shift the output to simulated baseline
    k_fit_shift = k_fit - np.repeat(baseline,curve_len).reshape(4,curve_len)
    
    # scale st they have minimum at -1
    #scale = np.min(k_fit_shift, axis=1)
    
    # scale st absolut value of 1
    scale = np.max(np.abs(k_fit_shift), axis=1)
    k_fit_scaled = np.array([k_fit_shift[i]/np.abs(scale[i]) for i in range(4)])
        
    # shift H_in
    if return_H_in:
        H_in = H_in_store[-1] - np.repeat(H_in_base,curve_len).reshape(3,curve_len)
        #print(H_in_base)
        return k_fit_scaled, H_in
    
    else:
        return k_fit_scaled #, baseline_store
#######################

def extract_params(params, mode, dc=np.ones(4)):
    """
    puts 1d parameter list into correct shape, depending on mode
    ---
    mode: indicating which HC is in the model:
    [['HC0', 'HC1', 'HC2'],
     ['HC0', 'HC1'],
     ['HC0', 'HC2'],
     ['HC1', 'HC2'],
     ['HC0'],
     ['HC1'],
     ['HC2'],
     ['special']]
    
    """
    
    # all HC
    if 'HC0' in mode and 'HC1' in mode and 'HC2' in mode:
        w00, w01, w02, w03, w11, w12, w13, w22, w23, a0,a1,a2 = params
        w = np.zeros((3,4)) # HC,cone
        w[0,0] = w00
        w[0,1] = w01
        w[0,2] = w02
        w[0,3] = w03
        w[1,1] = w11
        w[1,2] = w12
        w[1,3] = w13
        w[2,2] = w22
        w[2,3] = w23
    
    # HC0 and HC1
    elif 'HC0' in mode and 'HC1' in mode and not 'HC2' in mode:
        w00, w01, w02, w03, w11, w12, w13, a0, a1 = params
        w = np.zeros((3,4)) # HC,cone
        w[0,0] = w00
        w[0,1] = w01
        w[0,2] = w02
        w[0,3] = w03
        w[1,1] = w11
        w[1,2] = w12
        w[1,3] = w13
        #w[2,2] = w22
        #w[2,3] = w23
        a2 = 0
        
    # HC0 and HC2
    elif 'HC0' in mode and not 'HC1' in mode and 'HC2' in mode:
        w00, w01, w02, w03, w22, w23, a0, a2 = params
        w = np.zeros((3,4)) # HC,cone
        w[0,0] = w00
        w[0,1] = w01
        w[0,2] = w02
        w[0,3] = w03
        #w[1,1] = w11
        #w[1,2] = w12
        #w[1,3] = w13
        w[2,2] = w22
        w[2,3] = w23
        a1 = 0

    # HC1 and HC2
    elif not 'HC0' in mode and 'HC1' in mode and 'HC2' in mode:
        w11, w12, w13, w22, w23, a1,a2 = params
        w = np.zeros((3,4)) # HC,cone
        #w[0,0]= w00
        #w[0,1] = w01
        #w[0,2] = w02
        #w[0,3] = w03
        w[1,1] = w11
        w[1,2] = w12
        w[1,3] = w13
        w[2,2] = w22
        w[2,3] = w23
        a0 = 0 
        
    # HC0
    elif 'HC0' in mode and not 'HC1' in mode and not 'HC2' in mode:
        w00, w01, w02, w03, a0= params
        w = np.zeros((3,4)) # HC,cone
        w[0,0] = w00
        w[0,1] = w01
        w[0,2] = w02
        w[0,3] = w03
        #w[1,1] = w11
        #w[1,2] = w12
        #w[1,3] = w13
        #w[2,2] = w22
        #w[2,3] = w23
        a1 = 0
        a2 = 0
        
    # HC1
    elif not 'HC0' in mode and 'HC1' in mode and not 'HC2' in mode:
        w11, w12, w13, a1 = params
        w = np.zeros((3,4)) # HC,cone
        #w[0,0]= w00
        #w[0,1] = w01
        #w[0,2] = w02
        #w[0,3] = w03
        w[1,1] = w11
        w[1,2] = w12
        w[1,3] = w13
        #w[2,2] = w22
        #w[2,3] = w23
        a0 = 0  
        a2 = 0
    
    # HC2
    elif not 'HC0' in mode and not'HC1' in mode and 'HC2' in mode:
        w22, w23, a2 = params
        w = np.zeros((3,4)) # HC,cone
        #w[0,0]= w00
        #w[0,1] = w01
        #w[0,2] = w02
        #w[0,3] = w03
        #w[1,1] = w11
        #w[1,2] = w12
        #w[1,3] = w13
        w[2,2] = w22
        w[2,3] = w23
        a0 = 0 
        a1 = 0
        
    # HC0 and HC1 fully connected, others 'normal'
    elif mode == 'special':
        w00, w01, w02, w03, w10, w11, w12, w13, w22, w23, a0,a1,a2 = params
        w = np.zeros((3,4)) # HC,cone
        w[0,0] = w00
        w[0,1] = w01
        w[0,2] = w02
        w[0,3] = w03
        w[1,0] = w10
        w[1,1] = w11
        w[1,2] = w12
        w[1,3] = w13
        w[2,2] = w22
        w[2,3] = w23
    
    else:
        raise Warning('mode not implemented')
    
    a = np.zeros(3)
    a[0] = a0
    a[1] = a1
    a[2] = a2
    
    
    # fixed here
    #dc = np.ones(4)
    #dc[0] = dc0
    #dc[1] = dc1
    #dc[2] = 0.5
    #dc[3] = 2
    
    return w,dc,a



def get_param_labels(mode, dc=np.ones(4)):
    """
    returns list of str. for parameter names, depending on mode
    ---
    mode: indicating which HC is in the model:
    [['HC0', 'HC1', 'HC2'],
     ['HC0', 'HC1'],
     ['HC0', 'HC2'],
     ['HC1', 'HC2'],
     ['HC0'],
     ['HC1'],
     ['HC2']
     ['special']]
    
    """
    
    # all HC
    if 'HC0' in mode and 'HC1' in mode and 'HC2' in mode:
        labels = ['w00', 'w01', 'w02','w03', 'w11', 'w12', 'w13', 'w22', 'w23', 'a0','a1','a2']
       
    # HC0 and HC1
    elif 'HC0' in mode and 'HC1' in mode and not 'HC2' in mode:
        labels = ['w00', 'w01', 'w02','w03', 'w11', 'w12', 'w13', 'a0','a1']

    # HC0 and HC2
    elif 'HC0' in mode and not 'HC1' in mode and 'HC2' in mode:
        labels = ['w00', 'w01', 'w02','w03', 'w22', 'w23', 'a0','a2']

    # HC1 and HC2
    elif not 'HC0' in mode and 'HC1' in mode and 'HC2' in mode:
        labels = ['w11', 'w12', 'w13', 'w22', 'w23', 'a1','a2']
      
    # HC0
    elif 'HC0' in mode and not 'HC1' in mode and not 'HC2' in mode:
        labels = ['w00', 'w01', 'w02','w03', 'a0',]
       
    # HC1
    elif not 'HC0' in mode and 'HC1' in mode and not 'HC2' in mode:
         labels = [ 'w11', 'w12', 'w13', 'a1']
      
    # HC2
    elif not 'HC0' in mode and not'HC1' in mode and 'HC2' in mode:
        labels = ['w22', 'w23', 'a2']
        
     # HC0 and HC1 fully connected, others 'normal'
    elif mode == 'special':
        labels = ['w00', 'w01', 'w02', 'w03', 'w10', 'w11', 'w12', 'w13', 'w22', 'w23', 'a0','a1','a2']
       
    else:
        raise Warning('mode not implemented')
    
    return labels