#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 14:19:49 2020

@author: ingrid
"""

def curve_fitting(input, prop_clones):

    import pandas as pd
    import numpy as np
    
    data = pd.read_csv(input, index_col=0)
    if data.shape[0] < 1000:
        raise ValueError('The CurveFitting script has to be run on the Probs_for_Sample.csv file, which you get after specifying a sample size of interest using the sample = x option in the GenRis script.')
    numgenos = pd.to_numeric(list(data.columns))
    fitting = pd.DataFrame(numgenos)
    mean = np.mean(data, axis = 0)
    fitting['mean'] = mean.values
    perc = np.percentile(data.T, [16, 84], axis = 1) # 68% conf. interval = standard deviation
    perc = perc.T
    perc = pd.DataFrame(perc)
    fitting = pd.merge(fitting, perc, right_index=True, left_index=True)
    fitting = fitting.rename(columns={"0_x": "numgenos", "0_y": "lower_conf", 1:"upper_conf"})
    
    
    lower_conf = fitting[['numgenos', 'lower_conf']]
    lower_conf == 0
    lower_conf = lower_conf[~(lower_conf == 0).any(axis=1)]
    
    
    
    # curve fitting
    import scipy as sp
    from scipy.optimize import curve_fit
    def power_law(x, a, b):
        return a*np.power(x, b)
    # Fit the power-law data
    pars1, cov1 = curve_fit(f=power_law, xdata=fitting['mean'], ydata=fitting['numgenos'], p0=[0, 0], bounds=(-np.inf, np.inf))
    pars2, cov2 = curve_fit(f=power_law, xdata=lower_conf['lower_conf'], ydata=lower_conf['numgenos'], p0=[0, 0], bounds=(-np.inf, np.inf)) 
    pars3, cov3 = curve_fit(f=power_law, xdata=fitting['upper_conf'], ydata=fitting['numgenos'], p0=[0, 0], bounds=(-np.inf, np.inf))

    
    # Create figure and add axes object
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update(mpl.rcParamsDefault)
    fitting = pd.DataFrame(np.sort(fitting.values, axis=0), index=fitting.index, columns=fitting.columns)
    
    fig, ax = plt.subplots()
    ax.plot(fitting['mean'], power_law(fitting['mean'], *pars1), linewidth=2, color='black')
    ax.plot(fitting['lower_conf'], power_law(fitting['lower_conf'], *pars2), linestyle='--', linewidth=1, color='black')
    ax.plot(fitting['upper_conf'], power_law(fitting['upper_conf'], *pars3), linestyle='--', linewidth=1, color='black')
    ax.set(xlabel = "probability of picking clones", ylabel = "number of genotypes")
    fig.savefig("FittedCurve_ProbNumgenos.pdf")
    
    
    # calculate R square
    f = power_law
    xdata=fitting['mean']
    ydata=fitting['numgenos']
    popt, pcov = curve_fit(f, xdata, ydata)
    residuals = ydata-f(xdata, *popt)
    ss_res = np.sum(residuals**2) 
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R^2 of the curve fitted through the means is', r_squared,'.')

    
    
    # use function to calculate genotype number for specific clone picking probability
    est_mean = power_law(prop_clones, *pars1)
    est_low_conf = power_law(prop_clones, *pars2) 
    est_high_conf = power_law(prop_clones, *pars3) 
    print('With a proportion of', prop_clones, 'clones, you can expect around', est_mean, '(',est_low_conf, '< σ <', est_high_conf,') distinct genotypes in your sample.')

if __name__ == "__main__":
    curve_fitting()
