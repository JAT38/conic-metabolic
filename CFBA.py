def CFBA(Q):
    """
    Conic flux balance analysis, SOCP, as described in:
    
    Conic analysis of nonlinear metabolic networks
    J.A. Taylor, A. Rapaport, D. Dochain
    

    Copyright (C) 2023, Joshua A. Taylor
    taylor.a.josh@gmail.com
    """
    
    import numpy as np
    import cvxpy as cp
        
    v = cp.Variable(Q['r']) #reaction rates

    z = cp.Variable(Q['m']) #metabolite concentrations
    
    obj = cp.Maximize(Q['d'] @ z) #objective
        
    lcon = [
        Q['vL'] <= v, v <= Q['vU'], 
        Q['zL'] <= z, z <= Q['zU'], 
        Q['S'] @ v == 0
        ] #linear constraints
    
    lconv = [ Q['c'] @ v >= Q['vconlim'] ] #lower limit on c@v. for the COBRAPy textbook example, this is the biomass flux

    soccon = [
        cp.SOC(Q['mu'][i] * z[Q['sig'][Q['MM'][i]]] - Q['k'][i] * v[Q['MM'][i]] + Q['mu'][i] * Q['k'][i],
               cp.vstack((Q['mu'][i] * z[Q['sig'][Q['MM'][i]]] , Q['k'][i] * v[Q['MM'][i]] , Q['mu'][i] * Q['k'][i]))) for i in range(Q['MMlen'])
        ] #SOC constraints for the Michaelis-Menten kinetics
  
    prob = cp.Problem(obj, lcon + lconv + soccon) #form problem

    prob.solve(solver=cp.GUROBI,verbose=True) #solve the problem
    # prob.solve(solver=cp.GUROBI) #solve the problem
    
    
    #store SOC dual variables
    Phi = np.zeros((3,Q['MMlen']))
    phi = np.zeros(Q['MMlen'])
    for i in range(Q['MMlen']):
        temp = soccon[i].dual_value
        phi[i] = temp[0]
        Phi[:,[i]] = temp[1]

    #sensitivities to mu and k
    dFdmu = np.zeros(Q['MMlen'])
    dFdk = np.zeros(Q['MMlen'])
    for i in range(Q['MMlen']):
        dFdmu[i] = Phi[0,i] * z.value[Q['sig'][Q['MM'][i]]] + Phi[2,i] * Q['k'][i]\
            + phi[i] * (Q['k'][i] + z.value[Q['sig'][Q['MM'][i]]])
        dFdk[i] = Phi[1,i] * v.value[Q['MM'][i]] + Phi[2,i] * Q['mu'][i]\
            + phi[i] * (Q['mu'][i] - v.value[Q['MM'][i]])    
        

    sol = {'obj':obj.value,'v':v.value,'z':z.value,'dL':lcon[0].dual_value,\
           'dU':lcon[1].dual_value,'gL':lcon[2].dual_value,\
           'gU':lcon[3].dual_value,'lambda':lcon[4].dual_value,\
           'phi':phi,'Phi':Phi,'dFdmu':dFdmu,'dFdk':dFdk,'dv':lconv[0].dual_value}
            
    return sol