import spinmob as s
import mcphysics as mc
import numpy as np
import os

def fit(path, peak_domain=None, A0=None, b0=None, sigma0=None):
    """ Function to fit Gaussian to compton scattering data. 
    
    Parameters:
    ----------
    path: string, path of the directory containing the data files or path of a single data file
    peak_domain (optional): list of 2 int, domain of peak under consideration;
    A0 (optional): int, initial guess for parameter A;
    b0 (optional): int, initial guess for parameter b;
    sigma0 (optional): int, initial guess for parameter sigma;
    
    Returns:
    --------
    fit_parameters: array, parameters obtained through fitting [A, A.std, b, b.std, sigma, sigma.std];  
    """
    # Find the number of runs or if it's a single file
    try:
        runs = len(os.listdir(path))
    except:
        runs = 1
    f = s.data.fitter() # initiate fitter object
    
    # initiate parameter arrays
    A = np.zeros(runs) 
    A_std = np.zeros(runs)
    b = np.zeros(runs)
    b_std = np.zeros(runs)
    sigma = np.zeros(runs)
    sigma_std = np.zeros(runs)

    if peak_domain is not None:
        if not isinstance(peak_domain, list) or len(peak_domain) != 2:
            return TypeError('peak_domain must be a list of size 2.')

    if runs==1: # single run
        data = mc.data.load_chns(path) # load the single .chn data file
        f.set_functions(f = 'A * exp(-0.5*((x - b)/sigma)**2)/(sigma*sqrt(2*pi))', p = 'A='+str(A0)+',b='+str(b0)+',sigma='+str(sigma0))
        if peak_domain is None:
            f.set_data(xdata = data['Channel'], ydata = data['Counts'], xlabel='Channel', ylabel='Counts')
        else:
            f.set_data(xdata = data['Channel'][peak_domain[0]:peak_domain[1]], ydata = data['Counts'][peak_domain[0]:peak_domain[1]], xlabel='Channel', ylabel='Counts')
        f.fit() # fit to data
        A = f.get_fit_results()['A']
        A_std = f.get_fit_results()['A.std']
        b = f.get_fit_results()['b']
        b_std = f.get_fit_results()['b.std']
        sigma = f.get_fit_results()['sigma']
        sigma_std = f.get_fit_results()['sigma.std']

    else: # multiple runs
        data_set = mc.data.load_chns_directory(path)
        for i in range(0,runs):
            data = data_set[i]
            f.set_functions(f = 'A * exp(-0.5*((x - b)/sigma)**2)/(sigma*sqrt(2*pi))', p = 'A='+str(A0)+',b='+str(b0)+',sigma='+str(sigma0))
            if peak_domain is None:
                f.set_data(xdata = data['Channel'], ydata = data['Counts'], xlabel='Channel', ylabel='Counts')
            else:
                print(data)
                f.set_data(xdata = data['Channel'][peak_domain[0]:peak_domain[1]], ydata = data['Counts'][peak_domain[0]:peak_domain[1]], xlabel='Channel', ylabel='Counts')
            f.fit() # fit to data
            A[i] = f.get_fit_results()['A']
            A_std[i] = f.get_fit_results()['A.std']
            b[i] = f.get_fit_results()['b']
            b_std[i] = f.get_fit_results()['b.std']
            sigma[i] = f.get_fit_results()['sigma']
            sigma_std[i] = f.get_fit_results()['sigma.std']
    fit_parameters = [A, A_std, b, b_std, sigma, sigma_std]
    return fit_parameters

#__________________________________________________________________________________
