#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 14:19:49 2020

@author: ingrid
"""

def curve_fitting(input, prop_iMLG):

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
    fitting = fitting.rename(columns={"0_x": "numgenos", "0_y": "low_perc", 1:"up_perc"})
        
    lowPerc = fitting[['numgenos', 'low_perc']]
    lowPerc = lowPerc[lowPerc['low_perc'] != 0] # removes numgenos with low_perc=0
    
    upPerc = fitting[['numgenos', 'up_perc']]
    upPerc = upPerc[upPerc['up_perc'] != 0] # removes numgenos with up_perc=0
    
    
    # curve fitting
    import scipy as sp
    from scipy.optimize import curve_fit
    def power_law(x, a, b):
        return a*np.power(x, b)
    # Fit the power-law data
    pars1, cov1 = curve_fit(f=power_law, xdata=fitting['mean'], ydata=fitting['numgenos'], p0=[0, 0], bounds=(-np.inf, np.inf))
    if len(lowPerc.index) > 2: # curve fitting is not possible with too few data points
        pars2, cov2 = curve_fit(f=power_law, xdata=lowPerc['low_perc'], ydata=lowPerc['numgenos'], p0=[0, 0], bounds=(-np.inf, np.inf))
    if len(upPerc.index) > 2: # curve fitting is not possible with too few data points
        pars3, cov3 = curve_fit(f=power_law, xdata=upPerc['up_perc'], ydata=upPerc['numgenos'], p0=[0, 0], bounds=(-np.inf, np.inf))
    
    # Create figure and add axes object
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update(mpl.rcParamsDefault)    
    lowPerc = lowPerc.sort_values(by=['low_perc'], ascending=False) # sort values to get smooth curve
    upPerc = upPerc.sort_values(by=['up_perc'], ascending=False)
    
    fig, ax = plt.subplots()
    ax.plot(fitting['mean'], power_law(fitting['mean'], *pars1), linewidth=2, color='black')
    if len(lowPerc.index) > 2:
        ax.plot(lowPerc['low_perc'], power_law(lowPerc['low_perc'], *pars2), linestyle='--', linewidth=1, color='black')
    else: # the lower precentile for most/all data points is 0
        plt.axvline(x=0, linestyle='--', linewidth=1, color='black')
    if len(upPerc.index) > 2:
        ax.plot(upPerc['up_perc'], power_law(upPerc['up_perc'], *pars3), linestyle='--', linewidth=1, color='black')     
    ax.set(xlabel = "probability of picking identical MLGs", ylabel = "number of genotypes")
    fig.savefig("FittedCurve_ProbNumgenos.pdf")
    
    
    # calculate R square
    f = power_law
    xdata=fitting['mean']
    ydata=fitting['numgenos']
    popt, pcov = curve_fit(f, xdata, ydata)
    residuals = ydata-f(xdata, *popt)
    ss_res = np.sum(residuals**2) 
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = round(1 - (ss_res / ss_tot), 3)
    print('R^2 of the curve fitted through the means is', r_squared,'.')
    
        
    # use function to calculate genotype number for specific clone picking probability
    est_mean = round(power_law(prop_iMLG, *pars1), 2)
    if len(lowPerc.index) > 2:
        est_low_perc = round(power_law(prop_iMLG, *pars2), 2)
    else:
        est_low_perc = 0
    if len(upPerc.index) > 2:
        est_high_perc = round(power_law(prop_iMLG, *pars3), 2)
    else:
        est_high_perc = "unknown"
    print('With a proportion of', prop_iMLG, 'identical MLGs, you can expect around', est_mean, '(',est_low_perc, '< σ <', est_high_perc,') distinct genotypes in your sample.')
    if len(lowPerc.index) < 3 or len(upPerc.index) < 3:
        print('The probability of picking identical MLGs across the chosen number of genotypes is generally too close to 0 to accurately estimate lower and/or upper percentiles.')

if __name__ == "__main__":
    curve_fitting()
