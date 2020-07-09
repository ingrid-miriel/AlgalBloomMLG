#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:27:03 2020

@author: ingrid
"""
def function(numgenos, days, growthrate, sd, samplesize):
    import numpy as np
    if np.size(days) != 1 or np.size(growthrate) != 1:
        raise ValueError('Days and/or growthrate have to be defined each by exactly one value!')
    # define variables
    elif np.size(numgenos) > 1 and np.size(samplesize) > 1 and np.size(sd) == 1:
        numgenos = np.array(numgenos)
        samplesize = np.array(samplesize)
        sd = np.array([sd])
        x = samplesize
        y = numgenos
        x_title = "sample size"
        y_title = "number of genotypes"
        z = "NumgenosSamplesize"
    elif np.size(sd) > 1 and np.size(samplesize) > 1 and np.size(numgenos) == 1:
        sd = np.array(sd)
        numgenos = np.array([numgenos])
        samplesize = np.array(samplesize)
        x = samplesize
        y = sd
        x_title = "sample size"
        y_title = "SD of growth rate"
        z = "SamplesizeSDgrowth"
    elif np.size(sd) > 1 and np.size(numgenos) > 1 and np.size(samplesize) == 1:
        numgenos = np.array(numgenos)
        samplesize = np.array([samplesize])
        sd = np.array(sd)
        x = sd
        y = numgenos
        x_title = "SD of growth rate"
        y_title = "number of genotypes"
        z = "NumgenosSDgrowth"
    else:
        raise ValueError('Two variables (numgenos, sd or samplesize) need to be defined by a range of values! No more variables and no less!')

    import numpy as np
    import matplotlib.pyplot as plt
    import random
    import pandas as pd
    import seaborn as sns
    
    growthrate = growthrate + 1
    startamount = 1 # initial occurrence of all genotypes
    meanprobs3 = [] # set up final probability matrix

    # set up population matrix
    for a in range(len(numgenos)):
        pops=np.zeros((numgenos[a], days))
        pops[:,0] = startamount
    
        # set up growth rate matrix
        grow = np.zeros((numgenos[a], days))
        grow[:,1:days] = growthrate
        
        # set up probability matrices for intermediate loops
        meanprobs1 = []
        meanprobs2 = []
        
        for b in range(len(sd)):
            vargrow = np.random.normal(loc=growthrate, scale=sd[b], size=numgenos[a]) # create range of different growth rates from normal distribution
            vargrow = np.reshape(vargrow, (numgenos[a],1))                          
            grow[:,1:days] = np.tile(vargrow, (1, days-1)) # create array of equal size as pops filled with growthrates
        
            # loop to multiply 2 arrays (genotypes & growth rate) with each other
            for c in range(1,days):
                pops[:,c] = pops[:,c-1]*grow[:,c]
    
            # calculate proportion of each genotype in final population    
            finalprop = pops[:,days-1]/sum(pops[:,days-1])
        
            # Repeated picking (500 times) of clones with different sample sizes
            stats = []
            picks = [None] * len(samplesize)
            uniqpicks = [None] * len(samplesize)
            prop = [None] * len(samplesize)
            for d in range(0, 500):
                for e in range(len(samplesize)):
                    picks[e] = random.choices(population=range(0,numgenos[a]), # pick cells
                               weights=finalprop,
                               k=samplesize[e])
                    uniqpicks[e] = np.unique(picks[e], return_counts=True) # identify unique cells and clones
                    temp = samplesize[e] - ((uniqpicks[e][1] == 1).sum()) # count the occurrence of non-unique genotypes
                    prop[e] = temp/samplesize[e] # calculate proportion of clones in each sample (sample size)
                stats.append(prop.copy()) # append the proportions from each picking step to growing list
            temp2 = np.mean(stats, axis = 0) # calculate mean proportion of clones across 500 repeats
            meanprobs1.append(temp2) # append the mean proportions for each sd to growing list
        if np.size(samplesize) == 1:
            meanprobs2 = np.array(meanprobs1).T.tolist() # transpose
        if np.size(samplesize) > 1:
            meanprobs2 = meanprobs1
        if np.size(numgenos) > 1:
            meanprobs3.append(meanprobs2) # append mean proportions for each  population size (numgenos) to growing list
        if np.size(numgenos) == 1:
            meanprobs3 = meanprobs2
    
    # transform into numpy array            
    df = pd.DataFrame(np.array(meanprobs3).reshape(len(y), len(x)), index = y, columns = x)
    df.to_csv(z+'_GenRis_out.csv')
    df_percent = df * 100
    df_percent = df_percent.round()
    df_percent = df_percent.astype(int)
    
    #plot as heatmap
    sns.set(font_scale=0.6)
    fig, ax = plt.subplots()
    ax = sns.heatmap(df_percent, annot=True, fmt="d", linewidth=.2)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set(xlabel = x_title, ylabel = y_title)
    fig.savefig(z+"_Heatmap.pdf")