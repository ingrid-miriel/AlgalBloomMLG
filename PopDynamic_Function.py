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
    elif np.size(numgenos) > 1 and np.size(samplesize) > 1 and np.size(sd) == 1: # if more than one value in each argument
        return function1(numgenos, days, growthrate, sd, samplesize)
    elif np.size(sd) > 1 and np.size(samplesize) > 1 and np.size(numgenos) == 1:
        return function2(numgenos, days, growthrate, sd, samplesize)
    elif np.size(sd) > 1 and np.size(numgenos) > 1 and np.size(samplesize) == 1:
        return function3(numgenos, days, growthrate, sd, samplesize)    
    else:
        raise ValueError('Two variables (numgenos, sd or samplesize) need to be defined by a range of values! No more variables and no less!')

# NumgenosSamplesizePermutation
def function1(numgenos, days, growthrate, sd, samplesize):

    import numpy as np
    import matplotlib.pyplot as plt
    import random
    import pandas as pd
    import seaborn as sns
    
    numgenos = np.array(numgenos)
    samplesize = np.array(samplesize)
    startamount = 1 # initial occurrence of all genotype
    
    # set up population matrix
    samplesize_numgenos_probs = [None] * len(numgenos)
    for a in range(len(numgenos)):
        pops=np.zeros((numgenos[a], days))
        pops[:,0] = startamount
    
        # set up growth rate matrix
        grow = np.zeros((numgenos[a], days))
        grow[:,1:days] = growthrate
        vargrow = np.random.normal(loc=growthrate, scale=sd, size=numgenos[a]) # create range of different growth rates from normal distribution
        vargrow = np.reshape(vargrow, (numgenos[a],1))                          
        grow[:,1:days] = np.tile(vargrow, (1, days-1)) # create array of equal size as pops filled with growthrates
    
    
        # loop to multiply 2 arrays (genotypes & growth rate) with each other
        for b in range(1,days):
            pops[:,b] = pops[:,b-1]*grow[:,b]
    
        # calculate proportion of each genotype in final population    
        finalprob = pops[:,days-1]/sum(pops[:,days-1])
        
        # Repeated picking of clones with different sample sizes
        samplesize_stats = []
        picks = [None] * len(samplesize)
        uniqpicks = [None] * len(samplesize)
        prob = [None]*len(samplesize)
        for c in range(0, 500):
            for d in range(len(samplesize)):
                picks[d] = random.choices(population=range(0,numgenos[a]), # pick cells
                               weights=finalprob,
                               k=samplesize[d])
                uniqpicks[d] = np.unique(picks[d], return_counts=True) # identify unique cells and clones
                temp = samplesize[d] - ((uniqpicks[d][1] == 1).sum()) # counts the occurrence of non-unique genotypes
                prob[d] = temp/samplesize[d] # calculate proportion of clones in sample
            samplesize_stats.append(prob.copy())
            
        samplesize_numgenos_probs[a] = np.mean(samplesize_stats, axis = 0) # calculate mean proportion of clones across permutations
       
    df = pd.DataFrame(samplesize_numgenos_probs, index=numgenos, columns = samplesize) # reshape list into dataframe
    df.to_csv('NumgenosSamplesize_out.csv')
    df_percent = df * 100
    df_percent = df_percent.round()
    df_percent = df_percent.astype(int)
    
    #plot as heatmap
    sns.set(font_scale=0.6)
    fig, ax = plt.subplots()
    ax = sns.heatmap(df_percent, annot=True, fmt="d", linewidth=.2)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set(xlabel = "sample size", ylabel = "number of genotypes")
    fig.savefig("Test_NumgenosSamplesize_Heatmap.pdf")


# SamplesizeSDgrowPermutation
def function2(numgenos, days, growthrate, sd, samplesize):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    import pandas as pd
    import seaborn as sns

    sd = np.array(sd)
    samplesize = np.array(samplesize)
    startamount = 1 # initial occurrence of all genotype
    
    # set up population matrix
    pops = np.zeros((numgenos, days))
    pops[:,0] = startamount
    
    # set up growth rate matrix
    grow = np.zeros((numgenos, days))
    
    # permutate through different sd
    samplesize_sd_probs = [None] * len(sd)
    for a in range(len(sd)):
        vargrow = np.random.normal(loc=growthrate, scale=sd[a], size=numgenos) # create range of different growth rates from normal distribution
        vargrow = np.reshape(vargrow, (numgenos,1))                          
        grow[:,1:days] = np.tile(vargrow, (1, days-1)) # create array of equal size as pops filled with growthrates
    
        # loop to multiply 2 arrays (genotypes & growth rate) with each other
        for b in range(1,days):
            pops[:,b] = pops[:,b-1]*grow[:,b]
    
        # calculate proportion of each genotype in final population    
        finalprob = pops[:,days-1]/sum(pops[:,days-1])
        
        # Repeated picking of clones with different sample sizes
        samplesize_stats = []
        picks = [None] * len(samplesize)
        uniqpicks = [None] * len(samplesize)
        prob = [None]*len(samplesize)
        for c in range(0, 500):
            for d in range(len(samplesize)):
                picks[d] = random.choices(population=range(0,numgenos), # pick cells
                               weights=finalprob,
                               k=samplesize[d])
                uniqpicks[d] = np.unique(picks[d], return_counts=True) # identify unique cells and clones
                prob[d] = (samplesize[d] - ((uniqpicks[d][1] == 1).sum())) /samplesize[d] # calculate proportion of clones in sample
            samplesize_stats.append(prob.copy())
            
        samplesize_sd_probs[a] = np.mean(samplesize_stats, axis = 0) # calculate mean proportion of clones across permutations
       
    df = pd.DataFrame(samplesize_sd_probs, index=sd, columns = samplesize) # reshape list into dataframe
    df.to_csv('SamplesizeSDgrow_out.csv')
    df_percent = df * 100
    df_percent = df_percent.round()
    df_percent = df_percent.astype(int)
    
    # plot as heatmap
    sns.set(font_scale=0.6)
    fig, ax = plt.subplots()
    ax = sns.heatmap(df_percent, annot=True, fmt="d", linewidth=.2)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set(xlabel = "sample size", ylabel = "SD of growth rates")
    fig.savefig("Test_SampleSizeSDgrowth_Heatmap.pdf")
    
    
# NumgenosSDgrowthPermutation
def function3(numgenos, days, growthrate, sd, samplesize):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    import pandas as pd
    import seaborn as sns
    
    numgenos = np.array(numgenos)
    sd = np.array(sd)
    startamount = 1 # initial occurrence of all genotype

    # set up population matrix
    numgenos_sd_probs=np.zeros((len(numgenos), len(sd)))
    
    for a in range(len(numgenos)):
        pops=np.zeros((numgenos[a], days))
        pops[:,0] = startamount
    
        # set up growth rate matrix
        grow = np.zeros((numgenos[a], days))
        grow[:,1:days] = growthrate
        for b in range(len(sd)):
            vargrow = np.random.normal(loc=growthrate, scale=sd[b], size=numgenos[a]) # create range of different growth rates from normal distribution
            vargrow = np.reshape(vargrow, (numgenos[a],1))                          
            grow[:,1:days] = np.tile(vargrow, (1, days-1)) # create array of equal size as pops filled with growthrates
            
             # loop to multiply 2 arrays (genotypes & growth rate) with each other
            for c in range(1,days):
                pops[:,c] = pops[:,c-1]*grow[:,c]
    
            # calculate proportion of each genotype in final population    
            finalprob = pops[:,days-1]/sum(pops[:,days-1])
            
            # Repeated picking of clones
            sdgrow_stats = []
            for c in range(0, 500):
                picks = random.choices(population=range(0,numgenos[a]), weights=finalprob, k=samplesize)
                uniqpicks = np.unique(picks, return_counts=True)  
                prob = (samplesize - ((uniqpicks[1] == 1).sum())) / samplesize # calculate proportion of clones in sample
                sdgrow_stats.append(prob.copy())
    
            numgenos_sd_probs[a, b] = np.mean(sdgrow_stats, axis = 0) # calculate mean proportion of clones across permutations
            
    df = pd.DataFrame(numgenos_sd_probs, index=numgenos, columns = sd) # reshape list into dataframe
    df.to_csv('NumgenosSDgrow_out.csv')
    df_percent = df * 100
    df_percent = df_percent.round()
    df_percent = df_percent.astype(int)
    
    # plot as heatmap
    sns.set(font_scale=0.6)
    fig, ax = plt.subplots()
    ax = sns.heatmap(df_percent, annot=True, fmt="d", linewidth=.2)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set(xlabel = "SD of growth rates", ylabel = "number of genotypes") 
    fig.savefig("Test_NumgenosSDgrowth_Heatmap.pdf")
