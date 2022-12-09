#! /usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
from matplotlib import rc

from idealgas import F_idgas, p_idgas
from initialconditions import initial_cond
from fluxroe import flux_roe_ig
from plotting import prfls_matrix, prfl_graph_gen, trnd_graph_gen

rc('font', family='serif')
rc('lines', linewidth=1.5)
rc('font', size=14)



######################################################################
#                     1-D Euler System of Equations
#                            Roe Scheme Solver
#         Ideal gas based on ratio of specific heats Gamma
######################################################################


#===============================================================
# General Parameters
#===============================================================

gamma  = 1.4                # Ratio of specific heats
CFL    = 0.5                # Courant Number
ncells = 400                # Number of cells
x_ini =0; x_fin = 1         # Limits of computational domain
dx = (x_fin-x_ini)/ncells   # Step size
nx = ncells+1               # Number of points (boundary sections)
x = np.linspace(x_ini+dx/2,x_fin,nx)     # Mesh
plt_trend = 240             # customized trend plot position


#===============================================================
# Initial conditions & End time
#===============================================================

t  = 0

# The Sod's shock tube problem: IC=1
IC=1
r0 = np.zeros(nx)
u0 = np.zeros(nx)
p0 = np.zeros(nx)
halfcells = int(ncells/2)
p0,u0,r0,tEnd = initial_cond(IC,p0,u0,r0,halfcells)
TimeRange=np.array([[0]])

a0=F_idgas(gamma,p0,r0,u0)[0]
E0=F_idgas(gamma,p0,r0,u0)[1]
q=np.array([r0,r0*u0,r0*E0])      # Vector of conserved variables

# Plot of initial conditions
prfl_graph_gen(x,x_ini,x_fin,r0,u0,p0,E0,t)
profiles_r=np.array([r0.copy()]).T
profiles_u=np.array([u0.copy()]).T
profiles_p=np.array([p0.copy()]).T
profiles_E=np.array([E0.copy()]).T


#===============================================================
# Solver Loop
#===============================================================

it = 0
a  = a0
dt=CFL*dx/max(abs(u0)+a0)     # Using the system's largest eigenvalue

while t < tEnd:

    q0 = q.copy()
    dF = flux_roe_ig(gamma,q0,dx,a,nx)
    
    q[:,1:-2] = q0[:,1:-2]-dt/dx*dF
    q[:,0]=q0[:,0]; q[:,-1]=q0[:,-1]     # Dirichlet BCs
    
    # Compute primary variables
    rho=q[0]
    u=q[1]/rho
    E=q[2]/rho
    p=p_idgas(E,gamma,rho,u)
    a=F_idgas(gamma,p,rho,u)[0]
    if min(p)<0: print ('Negative pressure found !')
    
    # Update/correct time step
    dt=CFL*dx/max(abs(u)+a)
    
    # Update time and iteration counter
    t=t+dt; it=it+1
    TimeRange=np.concatenate((TimeRange,np.array([[t]])))
    
    # Data for trend plot
    profiles_r=prfls_matrix(rho,profiles_r)
    profiles_u=prfls_matrix(u,profiles_u)
    profiles_p=prfls_matrix(p,profiles_p)
    profiles_E=prfls_matrix(E,profiles_E)
    
    # Profile Plot
    if it%40 == 0:
        #print (it) 
        prfl_graph_gen(x,x_ini,x_fin,rho,u,p,E,t)

# Trend Plot
trnd_graph_gen(TimeRange,profiles_r,ncells,plt_trend,'rho (Density)','time (s)','density (kg/m³)')
trnd_graph_gen(TimeRange,profiles_u,ncells,plt_trend,'u (Gas Velocity)','time (s)','velocity (m/s)')
trnd_graph_gen(TimeRange,profiles_p,ncells,plt_trend,'p (Absolute Pressure)','time (s)','pressure (Pa-a)')
trnd_graph_gen(TimeRange,profiles_E,ncells,plt_trend,'E (Total Specific Energy)','time (s)','total spec energy (J/kg)')

