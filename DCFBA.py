def DCFBA(Q):
    """
    Dynamic conic flux balance analysis, SOCP, as described in:
    
    Conic analysis of nonlinear metabolic networks
    J.A. Taylor, A. Rapaport, D. Dochain
    

    Copyright (C) 2023, Joshua A. Taylor
    taylor.a.josh@gmail.com
    """
    
    # load packages
    import numpy as np
    import cvxpy as cp
        
    v = cp.Variable((Q['r'],Q['tau'])) #reaction rates
    
    z = cp.Variable((Q['m'],Q['tau']+1)) #metabolite concentrations

    x = cp.Variable(Q['tau']+1) #metabolite concentrations

    obj = cp.Maximize(x @ np.ones(Q['tau']+1)) #objective
    
    dyncon = [
            z[:,t+1] - z[:,t] == Q['dt'] * Q['S'] @ v[:,t] for t in range(Q['tau'])
        ] #dynamic constraints

    dynconx = [
            x[t+1] - x[t] == Q['dt'] * (v[Q['bmind'],t] - Q['bmdeath'][t] * x[t]) for t in range(Q['tau'])
        ] #biomass evolution

    ICcon = [
        x[0] == Q['x0'], z[:,0] == Q['z0']
        ] #initial conditions

    lcon1 = [
            x[t]*Q['vL'] <= v[:,t] for t in range(Q['tau'])
        ] #lower reaction limits

    lcon2 = [
            v[:,t] <= x[t]*Q['vU'] for t in range(Q['tau'])
        ] #upper reaction limits

    lcon3 = [
            Q['zL'] <= z[:,t] for t in range(Q['tau']+1)
        ] #lower metabolite limits

    lcon4 = [
            z[:,t] <= Q['zU'] for t in range(Q['tau']+1)
        ] #upper metabolite limits

    soccon = [
        cp.SOC(Q['mu'][i] * z[Q['sig'][Q['MM'][i]],t] - Q['k'][i] * v[Q['MM'][i],t] + Q['mu'][i] * Q['k'][i]*x[t],
               cp.vstack((Q['mu'][i] * z[Q['sig'][Q['MM'][i]],t] , Q['k'][i] * v[Q['MM'][i],t] , Q['mu'][i] * Q['k'][i]*x[t]))) for i in range(Q['MMlen']) for t in range(Q['tau'])
        ] #SOC constraints for reactions with single Michaelis-Menten kinetics, approximated by the Contois function
            

    prob = cp.Problem(obj, lcon1 + lcon2 + lcon3 + lcon4 + dyncon + dynconx + ICcon + soccon) #form problem

    prob.solve(solver=cp.GUROBI,verbose=True) #solve the problem
    # prob.solve(solver=cp.GUROBI) #solve the problem
    
    sol = {'obj':obj.value,'v':v.value,'z':z.value,'x':x.value}
            
    return sol