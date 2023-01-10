"""
Conic minimal cut set analysis, code for Section 5.3 in:

Conic analysis of nonlinear metabolic networks
J.A. Taylor, A. Rapaport, D. Dochain


Copyright (C) 2023, Joshua A. Taylor
taylor.a.josh@gmail.com
"""

# load packages 
import numpy as np
from cobra.io import load_model

# load data and model files
from ExEColiData import Cases
from CMCS import *
from CMCSlin import *

caseno = 1
Q = Cases(caseno) #import data for example

solCMCS = CMCS(Q) #solve for CMCS
solCMCSlin = CMCSlin(Q) #solve linear approximation for CMCS

print()

# Display the CMCS found in each case
model = load_model("textbook") #e coli model built into COBRAPy
print('CMCS found by MISOCP')
print('Reactions')
for i in range(Q['r']):
    if abs(solCMCS['binv'][i]) > 0.9:
        print(['BiGG ID: ',model.reactions[i].id])
        print(['Value of dual variable: ',solCMCS['delta'][i]])

print('Metabolites')
for i in range(Q['m']):
    if abs(solCMCS['binz'][i]) > 0.9:
        print(['BiGG ID: ',model.metabolites[i].id])
        print(['Value of dual variable: ',solCMCS['gamma'][i]])

print()


print('CMCS found by MILP')
print('Reactions')
for i in range(Q['r']):
    if abs(solCMCSlin['binv'][i]) > 0.9:
        print(['BiGG ID: ',model.reactions[i].id])
        print(['Value of dual variable: ',solCMCSlin['delta'][i]])

print('Metabolites')
for i in range(Q['m']):
    if abs(solCMCSlin['binz'][i]) > 0.9:
        print(['BiGG ID: ',model.metabolites[i].id])
        print(['Value of dual variable: ',solCMCSlin['gamma'][i]])
