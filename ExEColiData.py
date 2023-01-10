def Cases(caseno):
    """
    Data for examples in:

    Conic analysis of nonlinear metabolic networks
    J.A. Taylor, A. Rapaport, D. Dochain
    
    The data is from from COBRAPy and Table 1 in:

    Application of dynamic flux balance analysis to an industrial Escherichia coli fermentation
    A. L. Meadows, R. Karnik, H. Lam, S. Forestell, B. Snedecor
    Metabolic Engineering 12 (2) (2010) 150-160.


    Copyright (C) 2023, Joshua A. Taylor
    taylor.a.josh@gmail.com
    """
    
    #load packages
    import numpy as np
    import cobra
    from cobra.io import load_model
            
    if caseno == 1:
        
        model = load_model("textbook")   #e coli model built into COBRAPy     
                   
        S = cobra.util.create_stoichiometric_matrix(model) #stoichiometric matrix
                
        m = S.shape[0] #number of metabolites
        r = S.shape[1] #number of reactions
            
        c = np.zeros(r) #objective reaction vector
        d = np.zeros(m) #objective metabolite vector
        dd = -1 #objective coefficient to be assigned to certain metabolites later

        #list of reactions with MM kinetics and their associated metabolites
        react1 = ['ACKr','ACt2r','ACONTa','FUM','FRD7','FUMt2_2','FORt2','FORti','GLUN','GLNabc','GLUSy','GLCpts','SUCCt2_2','SUCCt3','SUCDi','SUCOAS','O2t'] #reactions with single MM
        metab1 = ['ac_c','ac_e','cit_c','fum_c','fum_c','fum_e','for_e','for_e','gln__L_c','gln__L_e','gln__L_c','glc__D_e','succ_e','succ_c','succ_c','succ_c','o2_e'] #metabolites in same order           
        
        MMlen = len(react1) #number of reactions with MM kinetics
        
        k = np.ones(MMlen) #vector of Michaelis constants
        mu = np.ones(MMlen) #vector of maximum reaction rates

        # parameters from Meadows et al (2010)

        k[0] = 6       
        k[1] = 6
        k[8] = 0.1
        k[9] = 0.1
        k[10] = 0.1
        k[11] = 0.53
        k[16] = 0.0024

        mu[0] = 2
        mu[1] = 2
        mu[6] = 2.2
        mu[7] = 2.2
        mu[11] = 10
        mu[16] = 20
   
        # upper and lower reaction limit vectors
        vL = np.zeros(r)
        vU = np.zeros(r)
        rev = []
        irrev = []

        for i in range(r):
            (vL[i],vU[i]) = model.reactions[i].bounds #take values from COBRAPy
            if vL[i] < -0.01: #reversible reactions
                rev = np.append(rev,i)
            else: #irreversible reactions
                irrev = np.append(irrev,i)

        rev = rev.astype(int)
        irrev = irrev.astype(int)
           
        # upper and lower metabolite limits
        zL = np.zeros(m)
        zU = zL + 1000*np.ones(m)

        vconlim = 0.1 # lower limit for c@v

        print()
        print('Reactions in objective:')
        for i in range(r):
            c[i] = model.reactions[i].objective_coefficient
            if c[i] > 0:
                bmind = i #index of biomass reaction for COBRAPy textbook example
                print(model.reactions[i].id)

        print()

        MM = [] #list of indices of reactions with MM kinetics
        
        print('Reactions with MM kinetics:')
        for j in range(MMlen):
            for i in range(r): #find indices of reactions with MM kinetics
                if react1[j] == model.reactions[i].id:
                    MM = np.append(MM,i)
                    print(model.reactions[i].name)
                            
        MM = MM.astype(int) #indices of reactions with MM kinetics
        
        print()
        
        sig = np.zeros(r) #element k is index of metabolite in MM kinetics for reaction k

        print('Associated metabolites:')
        for j in range(MMlen):            
            for i in range(m): #find metabolites indices
                if metab1[j] == model.metabolites[i].id:
                    sig[MM[j]] = i
                    d[i] = dd #assign value in objective vector
                    print(model.metabolites[i].name)
            
        sig = sig.astype(int)  
     
        print()   
     
        #Matrix M, as defined in paper                   
        M = np.zeros((MMlen,m))
                        
        for i in range(MMlen):
            M[i,sig[MM[i]]] = 1            

        noMM = np.setdiff1d(range(r),MM) #reactions without MM
                           

        ################################################################
        ################################################################
        #Parameters for dynamic CFBA
        dt = 0.001 #time step
        tau = 500 #number of time steps
        x0 = 1 #initial biomass concentration
        z0 = np.zeros(m) #initial concentrations of other metabolites
        bmdeath = 1 + np.sin((10/tau) * np.arange(tau)) #rate of biomass death

        ################################################################
        ################################################################
        #Parameters for CMCS analysis
        targets = ['FRD7'] #reaction to knockout
        dTv = len(targets)
        vt = 1
        Tv = np.zeros(r)
        for j in range(dTv):
            for i in range(r): #find indices of target reactions
                if targets[j] == model.reactions[i].id:
                    Tv[i] = 1

        # target region for metabolites (none in this example)
        zt = 0
        Tz = np.zeros(m)
        dTz = 1

        bigM = 100 #big M parameter for disjunctive constraints
        #weights for objective terms
        wz = 1
        wv = 1
        eta = np.divide(mu,k) #parameter for linear approximation to MM



    ################################################################
    ################################################################
    #put all parameters in Q
    Q = {'m':m,'r':r,'S':S,'c':c,'d':d,'vL':vL,'vU':vU,'zL':zL,'zU':zU,
             'mu':mu,'k':k,'sig':sig,'M':M,'vconlim':vconlim,
             'vt':vt,'zt':zt,'Tv':Tv,'Tz':Tz,'dTv':dTv,'dTz':dTz,'bigM':bigM,
             'MM':MM,'MMlen':MMlen,'noMM':noMM,'eta':eta,'wz':wz,'wv':wv,
             'rev':rev,'irrev':irrev,
             'dt':dt,'tau':tau,'bmind':bmind,'x0':x0,'z0':z0,'bmdeath':bmdeath
             }
            
    return Q