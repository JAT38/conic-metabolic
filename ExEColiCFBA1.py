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

Q['vconlim'] = 0.2

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

#sensitivities to V^max and K^m
sensmuc = sensmu[[2,11,16]]
senskc = solCFBA['dFdk'][[2,11,16]]

#plotting
nqc = len(sensmuc)

bwc = 0.3
vind1c = range(1,nqc+1)
vind2c = [x + 1.5*bwc for x in vind1c]
vind3c = [x + 1.5*bwc/2 for x in vind1c]

plt.subplots(figsize=(6, 4))
plt.bar(vind1c,sensmuc,width=0.3, color='black',edgecolor='black',label=r'$\frac{d\mathcal{F}}{dV^{\max}}$')
plt.bar(vind2c,senskc,width=0.3, color='white',edgecolor='black',label=r'$\frac{d\mathcal{F}}{dK^{m}}$')
plt.ylabel('Sensitivity',fontsize='large')
plt.xlabel('Reaction',fontsize='large')
plt.xticks(vind3c, ['3 (ACONTa)','12 (GLCpts)','17 (O2t)'])
plt.legend(fontsize='large')
plt.show()