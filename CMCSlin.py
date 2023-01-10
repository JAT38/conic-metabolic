def CMCSlin(Q):
    """
    Solve MILP for conic minimal cut set, as described in:
    
    Conic analysis of nonlinear metabolic networks
    J.A. Taylor, A. Rapaport, D. Dochain
    

    Copyright (C) 2023, Joshua A. Taylor
    taylor.a.josh@gmail.com
    """ 
    #load package
    import cvxpy as cp
        
    pv = cp.Variable(1) #reaction target dual

    pz = cp.Variable(1) #metabolite target dual
    
    delta = cp.Variable(Q['r']) #dual for v=0
    
    binv = cp.Variable(Q['r'], boolean=True) #binary variables for counting nonzero elements of v
    
    gamma = cp.Variable(Q['m']) #dual for z=0

    binz = cp.Variable(Q['m'], boolean=True) #binary variables for counting nonzero elements of z
    
    lam = cp.Variable(Q['m']) #stochiometric dual
        
    alpha = cp.Variable(Q['MMlen']) #linear approx to MM dual
    
    obj = cp.Minimize(Q['wv']*cp.sum(binv) + Q['wz']*cp.sum(binz)) #objective with just delta and gamma

    #disjunctive constraints
    lindisjrev = [
        cp.abs(delta[i]) <= Q['bigM']*binv[i] for i in Q['rev']
    ]
    lindisjirrev = [
        delta[i] <= Q['bigM']*binv[i] for i in Q['irrev']
    ]

    lindisjirrev0 = [
        delta[i] >= 0 for i in Q['irrev']
    ]

    lindisjgam = [
        gamma <= cp.multiply(Q['bigM'],binz)
    ]

    lindisjalpha = [
        alpha[i] <= cp.multiply(Q['bigM'],binz[Q['sig'][Q['MM'][i]]]) for i in range(Q['MMlen'])
    ]

    #most of the linear constraints in the dual system
    lincon = [
        Q['vt'] * pv + Q['zt'] * pz == 1,
        -Q['M'].T @ cp.multiply(alpha , Q['eta']) + gamma == cp.multiply(pz,Q['Tz']),
        pv >= 0, pz >=0,
        alpha >= 0, gamma >= 0,
        #constraint to keep knockout reaction/metabolite out of cut set
        cp.multiply(Q['Tv'],delta) == 0,
        cp.multiply(Q['Tz'],gamma) == 0
        ] #linear constraints

    linconMM = [
        alpha[i] + Q['S'][:,Q['MM'][i]].T @ lam + delta[Q['MM'][i]] == Q['Tv'][Q['MM'][i]] * pv  for i in range(Q['MMlen'])
        ]
    
    linconnoMM = [
        Q['S'][:,i].T @ lam + delta[i] == Q['Tv'][i] * pv for i in Q['noMM']
        ]


    prob = cp.Problem(obj, lindisjgam + lindisjalpha + lindisjrev + lindisjirrev + lindisjirrev0 + lincon + linconMM + linconnoMM) #form problem

    # prob.solve(solver=cp.GUROBI) #solve the problem
    prob.solve(solver=cp.GUROBI,verbose=True) #solve the problem
        
    sol = {'obj':obj.value,'pv':pv.value,'pz':pz.value,'lam':lam.value,\
            'delta':delta.value,'gamma':gamma.value,'alpha':alpha.value,\
            'binv':binv.value,'binz':binz.value}
            
    return sol