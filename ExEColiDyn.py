"""
Dynamic conic flux balance analysis, code for Section 5.2 in:

Conic analysis of nonlinear metabolic networks
J.A. Taylor, A. Rapaport, D. Dochain


Copyright (C) 2023, Joshua A. Taylor
taylor.a.josh@gmail.com
"""

# load packages
import numpy as np
import matplotlib.pyplot as plt
import time
import cobra
from cobra.io import load_model

# data and model files
from ExEColiData import Cases
from DCFBA import *

caseno = 1
Q = Cases(caseno)

nT = 10 #number of problem instances to solve
sT = 50 #time horizon of each instance is i*sT, for i=1,...,nT
times = np.zeros(nT) #record time taken
periods = np.zeros(nT) #number of periods


for i in range(nT):
    print()
    print('Iteration:')
    print(i)
    print()
    Q['tau'] = sT*(i+1) #time horizon
    Q['bmdeath'] = 0*np.ones(Q['tau'])
    start = time.time() #time solver
    solDCFBA = DCFBA(Q)
    times[i] = time.time() - start
    periods[i] = Q['tau']

#solution for last iteration, with longest time horizon
z = solDCFBA['z']
x = solDCFBA['x']
v = solDCFBA['v']

#indices of concentrations to plot
mindcit = Q['sig'][Q['MM'][2]]
mindglnc = Q['sig'][Q['MM'][8]]
mindglu = Q['sig'][Q['MM'][11]]

#concentrations to plot
zcit = z[mindcit,0:Q['tau']]
zglu = z[mindglu,0:Q['tau']]
zglnc = z[mindglnc,0:Q['tau']]

#indices of reactions to plot
vindcit = Q['MM'][2]
vindglun = Q['MM'][8]

#reactions to plot
vcit = v[vindcit,:]
vglun = v[vindglun,:]

#data for x axis
tt = np.arange(0, Q['dt'] * Q['tau'], Q['dt'])

#plotting
fig1, ax1 = plt.subplots(figsize=(6, 4))

ax1.plot(tt, x[0:Q['tau']], 'k:',label='Biomass')
ax1.plot(tt, zcit, 'k-',label='Citrate')
ax1.plot(tt, zglu, 'k--',label='Glucose')
ax1.plot(tt, zglnc, 'k-.',label='Glutamine')

ax1.set(xlim=(0, Q['dt'] * (Q['tau']-1)))

plt.xlabel('Hours',fontsize='large')
plt.ylabel('Concentration (mM)',fontsize='large')
plt.legend(fontsize='large')


fig2, ax2 = plt.subplots(figsize=(6, 4))

ax2.plot(periods,times,'k-')
plt.xlabel('Periods',fontsize='large')
plt.ylabel('Time (s)',fontsize='large')
ax2.set(xlim=(periods[0], periods[nT-1]))
plt.show()


ax3 = plt.subplot(211)
ax3.set(xlim=(0, Q['dt'] * (Q['tau']-1)))
ax3.plot(tt, x[0:Q['tau']], 'k:',label='Biomass')
ax3.plot(tt, zcit, 'k-',label='Citrate')
ax3.plot(tt, zglu, 'k--',label='Glucose')
ax3.plot(tt, zglnc, 'k-.',label='Glutamine')
plt.ylabel('Concentration (mM)',fontsize='large')
plt.legend(fontsize='large')
plt.tick_params('x', labelbottom=False)


ax4 = plt.subplot(212)
ax4.set(xlim=(0, Q['dt'] * (Q['tau']-1)))
ax4.plot(tt, vcit, 'k-',label='ACONTa')
ax4.plot(tt, vglun, 'k--',label='GLUN')
plt.ylabel('Flux (mmol/g dw/h)',fontsize='large')
plt.xlabel('Hours',fontsize='large')
plt.legend(fontsize='large')

plt.show()