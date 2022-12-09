# -*- coding: utf-8 -*-

import numpy as np
from idealgas import F_idgas, p_idgas
from funcflux import func_flux_ig



######################################################################
#                   1-D Euler System of Equations
#                         Roe Scheme Solver
#                      Ideal gas based - Gamma
######################################################################


def flux_roe_ig (gamma,q,dx,a,nx):

    # Compute primitive variables and enthalpy
    r=q[0]
    u=q[1]/r
    E=q[2]/r
    p=p_idgas(E,gamma,r,u)
    htot = F_idgas(gamma,p,r,u)[2]
    
    # Initialize Roe flux
    Phi=np.zeros((3,nx-1))
    
    for j in range (0,nx-1):
    
        # Compute Roe averages
        R=np.sqrt(r[j+1]/r[j])                       # R_{j+1/2}
        #rmoy=R*r[j]                                 # {hat rho}_{j+1/2} - not used
        umoy=(R*u[j+1]+u[j])/(R+1)                   # {hat U}_{j+1/2}
        hmoy=(R*htot[j+1]+htot[j])/(R+1)             # {hat H}_{j+1/2}
        amoy=np.sqrt((gamma-1.0)*(hmoy-0.5*umoy*umoy))     # {hat a}_{j+1/2}
        
        # Auxiliary variables used to compute P_{j+1/2}^{-1}
        alph1=(gamma-1)*umoy*umoy/(2*amoy*amoy)
        alph2=(gamma-1)/(amoy*amoy)

        # Compute vector (W_{j+1}-W_j)
        wdif = q[:,j+1]-q[:,j]
        
        # Compute matrix P^{-1}_{j+1/2}
        Pinv = np.array([[0.5*(alph1+umoy/amoy), -0.5*(alph2*umoy+1/amoy),  alph2/2],
                        [1-alph1,                alph2*umoy,                -alph2 ],
                        [0.5*(alph1-umoy/amoy),  -0.5*(alph2*umoy-1/amoy),  alph2/2]]);
                
        # Compute matrix P_{j+1/2}
        P    = np.array([[ 1,              1,              1              ],
                        [umoy-amoy,        umoy,           umoy+amoy      ],
                        [hmoy-amoy*umoy,   0.5*umoy*umoy,  hmoy+amoy*umoy ]]);
        
        # Compute matrix Lambda_{j+1/2}
        lamb = np.array([[ abs(umoy-amoy),  0,              0                 ],
                        [0,                 abs(umoy),      0                 ],
                        [0,                 0,              abs(umoy+amoy)    ]]);
                      
        # Compute Roe matrix |A_{j+1/2}|
        A=np.dot(P,lamb)
        A=np.dot(A,Pinv)
        
        # Compute |A_{j+1/2}| (W_{j+1}-W_j)
        Phi[:,j]=np.dot(A,wdif)
        
    # Compute Phi=(F(W_{j+1}+F(W_j))/2-|A_{j+1/2}| (W_{j+1}-W_j)/2
    F = func_flux_ig(gamma,q)
    Phi=0.5*(F[:,0:nx-1]+F[:,1:nx])-0.5*Phi
    dF = (Phi[:,1:-1]-Phi[:,0:-2])
    
    return (dF)