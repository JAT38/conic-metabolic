"""
Conic flux balance analysis, code for second part of Section 5.1 in:

Conic analysis of nonlinear metabolic networks
J.A. Taylor, A. Rapaport, D. Dochain


Copyright (C) 2023, Joshua A. Taylor
taylor.a.josh@gmail.com
"""

#import packages
import numpy as np
import matplotlib.pyplot as plt

#load data and models
from ExEColiData import Cases
from FBA import *
from CFBA import *

caseno = 1
Q = Cases(caseno)

nd = 10 #number of instances to solve
sc = 0.5 #increment for lower limit on biomass

#store optimal concentrations for each problem instance
solCFBAzcit = np.zeros(nd)
solCFBAzglu = np.zeros(nd)
solCFBAzo2 = np.zeros(nd)

#indices for metabolites of interest
mindcit = Q['sig'][Q['MM'][2]]
mindglu = Q['sig'][Q['MM'][11]]
mindo2 = Q['sig'][Q['MM'][16]]

for i in range(nd):
    Q['vconlim'] = sc*(i+1)/nd #Xi in example, lower limit on biomass flux

    solCFBA = CFBA(Q) #solve CFBA

    solCFBAzcit[i] = solCFBA['z'][mindcit]
    solCFBAzglu[i] = solCFBA['z'][mindglu]
    solCFBAzo2[i] = solCFBA['z'][mindo2]

#plotting
ddd = sc/nd
rnd = np.multiply(sc/nd,range(1,nd+1)) #data for x axis

fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(rnd, solCFBAzcit, 'k-',label='Citrate')
ax.plot(rnd, solCFBAzglu, 'k--',label='Glucose')
ax.plot(rnd, solCFBAzo2, 'k:',label='Oxygen')

ax.set(xlim=(sc/nd, sc), xticks=np.arange(sc/nd, sc, ddd),
       ylim=(-0.05, 1.2), yticks=np.arange(0, 2, 1))

plt.ylabel('Concentration (mM)',fontsize='large')
plt.xlabel(r'$\Xi$',fontsize='large')
plt.legend(fontsize='large')
plt.show()
