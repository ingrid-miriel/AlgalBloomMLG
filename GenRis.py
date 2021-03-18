#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:27:03 2020

@author: ingrid
"""
def genris(numgenos, days, growthrate, sd, samplesize, mincell, maxcell, **kwargs): 
    import numpy as np
    if np.size(days) != 1 or np.size(growthrate) != 1:
        raise ValueError('Days and/or growthrate have each to be defined by exactly one value!')
    # define variables
    elif np.size(numgenos) > 1 and np.size(samplesize) > 1 and np.size(sd) == 1:
        numgenos = np.array(numgenos)
        samplesize = np.array(samplesize)
        sd = np.array([sd])
        array = np.zeros((len(numgenos), len(samplesize), 10))
        x = samplesize
        y = numgenos
        x_title = "sample size"
        y_title = "number of genotypes"
        z = "NumgenosSamplesize"
    elif np.size(sd) > 1 and np.size(samplesize) > 1 and np.size(numgenos) == 1:
        sd = np.array(sd)
        numgenos = np.array([numgenos])
        samplesize = np.array(samplesize)
        array = np.zeros((len(sd), len(samplesize), 10))
        x = samplesize
        y = sd
        x_title = "sample size"
        y_title = "σ of growth rate"
        z = "SamplesizeSDgrowth"
    elif np.size(sd) > 1 and np.size(numgenos) > 1 and np.size(samplesize) == 1:
        numgenos = np.array(numgenos)
        samplesize = np.array([samplesize])
        array = np.zeros((len(numgenos), len(sd), 10))
        sd = np.array(sd)
        x = sd
        y = numgenos
        x_title = "σ of growth rate"
        y_title = "number of genotypes"
        z = "NumgenosSDgrowth"
    else:
        raise ValueError('Exactly two variables (numgenos, sd or samplesize) need to be defined by a range of values!')
    if np.size(sd) > 1 and 'sample' in kwargs:
        raise ValueError('The estimation of genotype richness for a specific sample size requires ranges of number of genotypes (numgenos) and samplesizes!')

    import numpy as np
    import matplotlib.pyplot as plt
    import random
    import pandas as pd
    import seaborn as sns
    
    if 'sample' in kwargs:
        [[num]] = np.where(samplesize == kwargs['sample'])
        sampleprobs2 = pd.DataFrame(columns = numgenos)

    for i in range(0, 10):
        sampleprobs = np.zeros((100, len(numgenos)))
        meanprobs3 = [] # set up final probability matrix
        # set up population matrix
        for a in range(len(numgenos)):        
            pops=np.zeros((numgenos[a], 4))
            startamount = []
            
            for j in range(0, numgenos[a]):
                startamount.append(random.randint(mincell, maxcell))      
            pops[:, 0] = startamount
            
            # set up probability matrices for intermediate loops
            meanprobs1 = []
            meanprobs2 = []
            
            # assign growthrates to genotypes, calculate proportion of genotypes after exponential growth
            for b in range(len(sd)):
                pops[:, 1] = np.array(np.random.normal(loc=growthrate, scale=sd[b], size=numgenos[a]))
                pops[:, 2] = pops[:, 0] * np.exp(pops[:, 1] * days)
                pops[:, 3] = pops[:, 2]/sum(pops[:, 2])
            
                # Repeated picking (100 times) of clones with different sample sizes
                stats = []
                picks = [None] * len(samplesize)
                uniqpicks = [None] * len(samplesize)
                prop = [None] * len(samplesize)
                for d in range(0, 100):
                    for e in range(len(samplesize)):
                        picks[e] = random.choices(population=range(0, numgenos[a]), # pick cells
                                   weights=pops[:, 3],
                                   k=samplesize[e])

                        uniqpicks[e] = np.unique(picks[e], return_counts=False) # identify distinct genotypes
                        prop[e] = (samplesize[e] - len(uniqpicks[e]))/samplesize[e] # calculate proportion of clones in each sample (sample size)       
                    stats.append(prop.copy()) # append the proportions from each picking step to growing list
                if 'sample' in kwargs:
                    sampleprobs[:, a] = [row[num] for row in stats]
                meanprobs1.append(np.mean(stats, axis = 0)) # calculate mean proportion of clones across 100 repeats and append it for each sd to growing list
            if np.size(samplesize) == 1:
                meanprobs2 = np.array(meanprobs1).T.tolist() # transpose
            if np.size(samplesize) > 1:
                meanprobs2 = meanprobs1
            if np.size(numgenos) > 1:
                meanprobs3.append(meanprobs2) # append mean proportions for each  population size (numgenos) to growing list
            if np.size(numgenos) == 1:
                meanprobs3 = meanprobs2
        if 'sample' in kwargs:
            sampleprobs1 = pd.DataFrame(sampleprobs, columns = y)
            sampleprobs2 = sampleprobs2.append(sampleprobs1, ignore_index=True)
        # transform into numpy array
        out = pd.DataFrame(np.array(meanprobs3).reshape(len(y), len(x)))
        array[:, :, i] = out # fill 3D array
        print((i+1)*10,'%')
    if 'sample' in kwargs:
        sampleprobs2.to_csv('Probs_for_Sample.csv')
    df = np.mean(array, axis=2) # calculate means of 100 repetions
    df = pd.DataFrame(df, index = y, columns = x)
    df.to_csv(z+'_GenRis.csv')
    
    df_percent = df * 100
    df_percent = df_percent.round()
    df_percent = df_percent.astype(int)
    
    #plot as heatmap
    sns.set(font_scale=0.6)
    fig, ax = plt.subplots()
    ax = sns.heatmap(df_percent, annot=True, fmt="d", linewidth=.2, cbar_kws={'label': 'probability of picking clones (%)'})
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set(xlabel = x_title, ylabel = y_title)
    fig.savefig(z+"_Heatmap.pdf")
    
    
if __name__ == "__main__":
    genris()
