def CMCS(Q):
    """
    Solve MISOCP for conic minimal cut set, as described in:
    
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
    
    Phi = cp.Variable((3,Q['MMlen'])) #SOC dual
    
    phi = cp.Variable(Q['MMlen']) #SOC dual
    
    eps = cp.Variable(Q['MMlen']) #for SOC norm

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
        eps[i] <= cp.multiply(Q['bigM'],binz[Q['sig'][Q['MM'][i]]]) for i in range(Q['MMlen'])
    ]

    #most of the linear constraints in the dual system
    lincon = [
        (Phi[2,:]-phi.T) @ cp.multiply(Q['k'],Q['mu']) + Q['vt'] * pv + Q['zt'] * pz == 1,
        Q['M'].T @ cp.multiply(Phi[0,:].T - phi , Q['mu']) + gamma == cp.multiply(pz,Q['Tz']),
        pv >= 0, pz >=0,
        phi >= 0, gamma >= 0,
        #constraint to keep knockout reaction/metabolite out of cut set
        cp.multiply(Q['Tv'],delta) == 0,
        cp.multiply(Q['Tz'],gamma) == 0
        ] #linear constraints

    linconMM = [
        (Phi[1,i] + phi[i]) * Q['k'][i] + Q['S'][:,Q['MM'][i]].T @ lam + delta[Q['MM'][i]] == Q['Tv'][Q['MM'][i]] * pv  for i in range(Q['MMlen'])
        ]
    
    linconnoMM = [
        Q['S'][:,i].T @ lam + delta[i] == Q['Tv'][i] * pv for i in Q['noMM']
        ]

    linconeps = [
        cp.norm(phi[i],1) + cp.norm(Phi[:,i],2) <= eps[i]  for i in range(Q['MMlen'])
        ]

    soccon = [
        cp.SOC(phi[i],Phi[:,i]) for i in range(Q['MMlen'])
        ] #SOC constraints

    prob = cp.Problem(obj, lindisjgam + lindisjalpha + lindisjrev + lindisjirrev + lindisjirrev0 + lincon + linconMM + linconnoMM + linconeps + soccon) #form problem

    prob.solve(solver=cp.GUROBI,verbose=True) #solve the problem
    # prob.solve(solver=cp.GUROBI) #solve the problem

    sol = {'obj':obj.value,'pv':pv.value,'pz':pz.value,'lam':lam.value,\
           'delta':delta.value,'gamma':gamma.value,'phi':phi.value,'Phi':Phi.value,\
            'binv':binv.value,'binz':binz.value}
            
    return sol