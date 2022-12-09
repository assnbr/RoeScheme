# -*- coding: utf-8 -*-


""" INITIAL CONDITIONS """

#===============================================================
# Choose the case
# 6 cases of initial conditions are available
#===============================================================

def initial_cond (ICnd,Pp0,Uu0,Rr0,halfclls):
    
    if ICnd == 1:
        print ("Configuration 1, Sod's Problem")
        Pp0[:halfclls] = 1.0  ; Pp0[halfclls:] = 0.1;
        Uu0[:halfclls] = 0.0  ; Uu0[halfclls:] = 0.0;
        Rr0[:halfclls] = 1.0  ; Rr0[halfclls:] = 0.125;
        tEnd_sim = 0.20;
    elif ICnd== 2:
        print ("Configuration 2, Left Expansion and right strong shock")
        Pp0[:halfclls] = 1000.; Pp0[halfclls:] = 0.1;
        Uu0[:halfclls] = 0.0  ; Uu0[halfclls:] = 0.0;
        Rr0[:halfclls] = 3.0  ; Rr0[halfclls:] = 0.2;
        tEnd_sim = 0.01;
    elif ICnd == 3:
        print ("Configuration 3, Right Expansion and left strong shock")
        Pp0[:halfclls] = 7.   ; Pp0[halfclls:] = 10.;
        Uu0[:halfclls] = 0.0  ; Uu0[halfclls:] = 0.0;
        Rr0[:halfclls] = 1.0  ; Rr0[halfclls:] = 1.0;
        tEnd_sim = 0.10;
    elif ICnd == 4:
        print ("Configuration 4, Shocktube problem of G.A. Sod, JCP 27:1, 1978")
        Pp0[:halfclls] = 1.0  ; Pp0[halfclls:] = 0.1;
        Uu0[:halfclls] = 0.75 ; Uu0[halfclls:] = 0.0;
        Rr0[:halfclls] = 1.0  ; Rr0[halfclls:] = 0.125;
        tEnd_sim = 0.17;
    elif ICnd == 5:
        print ("Configuration 5, Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997")
        Pp0[:halfclls] = 3.528; Pp0[halfclls:] = 0.571;
        Uu0[:halfclls] = 0.698; Uu0[halfclls:] = 0.0;
        Rr0[:halfclls] = 0.445; Rr0[halfclls:] = 0.5;
        tEnd_sim = 0.15;
    elif ICnd == 6:
        print ("Configuration 6, Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997")
        Pp0[:halfclls] = 10.33; Pp0[halfclls:] = 1.0;
        Uu0[:halfclls] = 0.92 ; Uu0[halfclls:] = 3.55;
        Rr0[:halfclls] = 3.857; Rr0[halfclls:] = 1.0;
        tEnd_sim = 0.09;
    
    return (Pp0,Uu0,Rr0,tEnd_sim)

