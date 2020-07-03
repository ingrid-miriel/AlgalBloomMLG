# PickingClones
simulating the probabiliy to pick clones and unique genotypes in phytoplankton populations

You can run this simulation in python with the following commands:

> from PopDynamic_Function import function
> function(numgenos = v, days = w, growthrate = x, sd = y, samplesize = z)

The PopDynamics_Function.py file needs to be stored in the same folder, where you want to create your output files.
The simulation perumatates through a range of different population sizes (numgenos), samplesizes and standard deviations of the normally distributed growth rate function (sd). Ranges have to be expressed as lists (e.g. numgenos = [200, 400, 600, 800, 1000, 1500, 2000, 3000]), while discrete values can be filled in directly.
Only two variables with ranges of values are allowed per run! 

Example:

> from PopDynamic_Function import function
> function(numgenos = [200, 400, 600, 800, 1000, 1500, 2000, 3000], days = 60, growthrate = 1.2, sd = 0.025, samplesize = [50, 100, 150, 200, 300, 400, 500, 1000])

> from PopDynamic_Function import function
> function(numgenos = 2000, days = 60, growthrate = 1.2, sd = [0.021, 0.023, 0.025, 0.027, 0.029], samplesize = [50, 100, 150, 200, 300, 400, 500, 1000])
