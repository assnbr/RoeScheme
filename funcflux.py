# -*- coding: utf-8 -*-

""" FUNCTION FLUX - IDEAL GAS"""

import numpy as np
from idealgas import p_idgas


def func_flux_ig (gamma,q):
    
    # Primitive variables
    r=q[0]
    u=q[1]/r
    E=q[2]/r
    p=p_idgas(E,gamma,r,u)
    
    # Flux vector
    F0 = np.array(r*u)
    F1 = np.array(r*u**2+p)
    F2 = np.array(u*(r*E+p))
    flux=np.array([F0,F1,F2])
    
    return (flux)

