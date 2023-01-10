def FBA(Q):
    """
    Flux balance analysis LP
    
    Conic analysis of nonlinear metabolic networks
    J.A. Taylor, A. Rapaport, D. Dochain
    
    
    Copyright (C) 2023, Joshua A. Taylor
    taylor.a.josh@gmail.com
    """
    
    import numpy as np
    import cvxpy as cp
        
    v = cp.Variable(Q['r']) #reaction rates
    
    obj = cp.Maximize(Q['c'] @ v) #objective
    
    con = [
        Q['vL'] <= v, v <= Q['vU'], 
        Q['S'] @ v == 0
        ] #constraints
    
    prob = cp.Problem(obj, con) #form problem

    prob.solve(solver=cp.GUROBI) #solve the problem

    sol = {'obj':obj.value,'v':v.value,'dualvL':con[0].dual_value,'dualvU':con[1].dual_value}
        
    return sol