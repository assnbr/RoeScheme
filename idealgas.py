# -*- coding: utf-8 -*-

import numpy as np

"""
IDEAL GAS - Caloric Equation of State

Notes:
- Gas ideal is P/Pcr<0.05 or T/Tcr>2.0
- An ideal gas with constant specific heats is referred to as a perfect gas.
- For air, the specific heat does not vary more than approximately 12% over a few hundred Kelvin in the range T > 270K.
- We must always take note that the value of the adiabatic index (Gama) must always be greater than one. For example, air it has k = 1.3 for hot-standard and k = 1.4 for cold-standard, respectively.
- The values for Cv and Cp as a function of temperature have been tabulated for various gases in the back of some textbooks, e.g., Table A-20 in Moran et al., 8th ed.
- See Book2 for: Thermal Equation of State.
- See Book3 for: Molar Specific Heat for some gases.
"""

# Beta - number of degrees of freedom (for simplicity, assume that Beta is constant)
# Beta = 3 for monatomic gas (Helium, Neon, Argon or Krypton)
# Beta = 5 for diatomic gas
# Beta = 6 for multiatomic gas
# Gama - ratio of specific heats - adiabatic index [-]
# Gama=1+2/Beta
#   air: Gama=1.40 for cold-standard
#   isothermal gas: Gama=1.0



""" Thermodynamic Properties of Ideal Gas """

# aig - adiabatic speed of sound, isentropic processes [m/s]
# eig - specific internal energy [J/kg]
# Gamaig - adiabatic index [-]
# hig - specific enthalpy [J/kg]
# prig - absolute pressure for a calorically ideal gas [Pa]
# rhoig - density [kg/m³]

def PVT_idealgas(Gamaig,prig,rhoig):
    #prig=(Gamaig-1)*rhoig*eig
    eig=prig/((Gamaig-1)*rhoig)
    aig=np.sqrt(Gamaig*prig/rhoig)
    hig=eig+prig/rhoig
    return (aig,eig,hig)



""" Properties of the flux with ideal gas """

# ag - adiabatic speed of sound, isentropic processes [m/s]
# Etg - total specific energy [J/kg]
# Gama - adiabatic index [-]
# Htg - total specific enthalpy [J/kg]
# prg - absolute pressure for a calorically ideal gas [Pa]
# rhog - density [kg/m³]
# ug - gas velocity [m/s]

def F_idgas(Gama,prg,rhog,ug):
    ag,eg,hg=PVT_idealgas(Gama,prg,rhog)
    Etg=eg+0.5*ug**2            
    Htg=Etg+prg/rhog
    return (ag,Etg,Htg)

def p_idgas(Etg,Gama,rhog,ug):
    prg=(Gama-1)*rhog*(Etg-0.5*ug**2)
    return (prg)


