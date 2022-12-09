# -*- coding: utf-8 -*-

"""
STIFFENED GAS
Stiffened Gas Equation of State

- This equation of state is based on the ideal gas law, to which a factor is added to reduce the compressibility.
- This EoS can represent liquid-like fluids in addition to gases.
"""


""" Thermodynamic Properties of Stiffened Gas """

# Gama - adiabatic index [-]
# Cp - heat capacity at constant pressure [J/kg.K]
# p_inf - pressure representing the molecular attraction between molecules [Pa]
# e_sg - specific internal energy [J/kg]
# p_sg - absolute pressure [Pa]
# Grun_sg - First Grüneisen coefficient for stiffened gas EoS [-]


def PVT_stiffenedgas(Cp_sg, Gama_sg, p_inf, rho_sg, T_sg):
    # index 0: gas
    # index 1: liquid
    e_sg=Cp_sg*T_sg/Gama_sg+p_inf/rho_sg
    p_sg=(Gama_sg-1)*rho_sg*Cp_sg*T_sg/Gama_sg-p_inf
    Grun_sg=Gama_sg-1
    return (e_sg, Grun_sg, p_sg)

