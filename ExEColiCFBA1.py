"""
Conic flux balance analysis, code for first part of Section 5.1 in:

Conic analysis of nonlinear metabolic networks
J.A. Taylor, A. Rapaport, D. Dochain


Copyright (C) 2023, Joshua A. Taylor
taylor.a.josh@gmail.com
"""

#import packages
import numpy as np
import matplotlib.pyplot as plt
import time

#load data and models
from ExEColiData import Cases
from FBA import *
from CFBA import *

caseno = 1
Q = Cases(caseno)

sol0 = CFBA(Q) #run CFBA to get solvers going so time comparison is fair

startb = time.time() #time solvers
solFBA = FBA(Q) #run FBA
print("FBA time: %s" % (time.time() - startb))

starta = time.time()
solCFBA = CFBA(Q) #run CFBA
print("CFBA time: %s" % (time.time() - starta))

print()

#sensitivities to V^max and K^m
sensmu = solCFBA['dFdmu']
sensk = solCFBA['dFdk']

#plotting
nq = len(sensmu)

bw = 0.25
vind1 = range(1,nq+1)
vind2 = [x + 1.5*bw for x in vind1]

plt.subplots(figsize=(6, 4))
plt.bar(vind1,sensmu,width=0.25, color='black',edgecolor='black',label=r'$\frac{d\mathcal{F}}{d\mu}$')
plt.bar(vind2,sensk,width=0.25, color='white',edgecolor='black',label=r'$\frac{d\mathcal{F}}{d\kappa}$')
plt.ylabel('Sensitivity',fontsize='large')
plt.xlabel('Reaction number',fontsize='large')
plt.xticks(range(1,nq+1))
plt.legend(fontsize='large')
plt.show()