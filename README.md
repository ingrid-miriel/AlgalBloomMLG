# PickingClones
Simulating the probabiliy of picking clones and unique genotypes in phytoplankton populations

You can run this simulation in python with the following commands:

> import os

> os.getcwd() # check your current working directory

> os.chdir('path/to/folder/with/GenRiS.py') # set working directory where the GenRiS.py file is located and where the output files will be generated


> from GenRiS import function

> function(numgenos = v, days = w, growthrate = x, sd = y, samplesize = z)

The GenRiS.py file needs to be stored in the same folder, where you want to create your output files.
The simulation perumatates through a range of different population sizes (number of genotypes = numgenos), samplesizes (number of cells that will be picked) and standard deviations of the normally distributed growth rate function (sd). Larger sd will result in larger differences between growth rates of individual genotypes. Ranges have to be expressed as lists (e.g. numgenos = [200, 400, 600, 800, 1000, 1500, 2000, 3000]), while discrete values can be filled in directly.
Only two variables with ranges of values are allowed per run! 
The growth rate has to be expressed in divisions per day and the variable "day" refers to time period of interest, in which the algal population grows exponentially.

Examples:

> from GenRiS import function

> function(numgenos = [200, 400, 600, 800, 1000, 1500, 2000, 3000], days = 60, growthrate = 0.2, sd = 0.025, samplesize = [50, 100, 150, 200, 300, 400, 500, 1000])


> from GenRiS import function

> function(numgenos = 2000, days = 60, growthrate = 0.2, sd = [0.021, 0.023, 0.025, 0.027, 0.029], samplesize = [50, 100, 150, 200, 300, 400, 500, 1000])
