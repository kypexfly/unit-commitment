""" ECONOMIC LOAD DISPACH BY:
* ITERATIVE LAMBDA METHOD w/o LOSSES
* BINARY SEARCH METHOD w/o LOSSES
* ITERATIVE WITH LOSSES (IMPORTANT: POWER LIMIT ISN'T CONSIDERED)

USED LIBRARIES: numpy, sympy (diff, symbols, poly)

CODE BY: Ricardo Chu Zheng, 2022
GITHUB: https://github.com/kypexfly
REFERENCE: Power Generation, Operation, and Control, Allen J. Wood, Bruce F. Wollenberg, Gerald B. Sheblé (2013) """

import numpy as np
from sympy import diff, symbols, poly

def binary_search(data:np.ndarray, PD:float):
    
    NGEN = data.shape[0]
    PMIN = data[:, 3]
    PMAX = data[:, 4]
    PGEN = 0
    if(PD < np.sum(PMIN)):
        return("PD is less than PMIN. ELD cannot be satisfied, choose another value of PD.")
    elif(PD > np.sum(PMAX)):
        return("PD is larger than PMAX. ELD cannot be satisfied, choose another value of PD.")
    a = data[:, 0]
    b = data[:, 1]
    dCmin = 2*a * PMIN + b
    dCmax = 2*a * PMAX + b
    minLambda = np.min(dCmin)
    maxLambda = np.max(dCmax)
    deltaLambda = (maxLambda - minLambda)/2
    Lambda = minLambda + deltaLambda
    iter = 0

    while abs(PGEN - PD) > 1e-6:
        P = (Lambda - b)/(2*a)
        for i in range(NGEN):
            if Lambda < dCmin[i]:
                P[i] = PMIN[i]
            elif Lambda > dCmax[i]:
                P[i] = PMAX[i]
        PGEN = np.sum(P)
        deltaLambda = deltaLambda/2
        if PGEN > PD: Lambda = Lambda - deltaLambda
        elif PGEN < PD: Lambda = Lambda + deltaLambda
        iter = iter + 1
    
    return P

def iterative_lambda(data:np.ndarray, PD:float):
    NGEN = data.shape[0]
    PMIN = data[:, 3]
    PMAX = data[:, 4]
    
    if(PD < np.sum(PMIN)):
        return("PD is less than PMIN. ELD cannot be satisfied, choose another value of PD.")
    elif(PD > np.sum(PMAX)):
        return("PD is larger than PMAX. ELD cannot be satisfied, choose another value of PD.")
    
    a = data[:, 0]
    b = data[:, 1]
    Lambda = np.max(b)
    P = np.zeros(NGEN)
    err = PD
    d_err = -sum(1/(2*a))
    iter = 0
    while abs(err) > 1e-6:
        P = (Lambda - b)/(2*a)
        for i in range(NGEN):
            if P[i] < PMIN[i]:
                P[i] = PMIN[i]
            elif P[i] > PMAX[i]:
                P[i] = PMAX[i]
        err = PD - sum(P)
        Lambda = Lambda - err/d_err
        iter = iter + 1

    return P

def iterative_w_losses(data:np.ndarray, PD:float, Ploss:str):
    """IMPORTANT: POWER LIMIT ISN'T CONSIDERED"""
    NGEN = data.shape[0]
    PMIN = data[:, 3]
    PMAX = data[:, 4]
    
    if(PD < np.sum(PMIN)):
        return("PD is less than PMIN. ELD cannot be satisfied, choose another value of PD.")
    elif(PD > np.sum(PMAX)):
        return("PD is larger than PMAX. ELD cannot be satisfied, choose another value of PD.")
    
    Pvars = symbols([f"P{i}" for i in range(NGEN)])
    a = data[:, 0]
    b = data[:, 1]
    P = np.random.randint(low=PMIN, high=PMAX, size=NGEN)
    err = PD
    iter = 0
    A = np.zeros((NGEN + 1, NGEN + 1))
    diag_index = np.diag_indices(NGEN)
    A[diag_index] = 2*a
    A[-1,:-1] = 1
    B = np.append(np.array(b*-1), PD + float(poly(Ploss).eval(list(P))))
    while(err > 1e-6):
        Pold = P
        Pdict = dict(zip(Pvars, Pold))
        for i in range(NGEN):
            positions = diff(Ploss, Pvars[i]).as_poly().args[1:]
            P2 = [Pdict[i] for i in positions]
            A[:-1, -1][i] = poly(diff(Ploss, Pvars[i])).eval(P2) - 1
            
        positions = poly(Ploss).args[1:]
        P2 = [Pdict[i] for i in positions]
        B[-1] = PD + poly(Ploss).eval(list(P))
        result = np.linalg.solve(A, B)
        P = result[:NGEN]
        # Lambda = result[-1]
        err = np.sum((P-Pold)**2)
        iter = iter + 1
    
    return P