# -*- coding: utf-8 -*-

""" PLOTTING """

import numpy as np
import matplotlib.pyplot as plt


# Profile Graphic Generator

def prfl_graph_gen (x,x_ini,x_fin,rho,u,p,E,tm):

    fig,axes = plt.subplots(nrows=4, ncols=1)
    plt.subplot(4, 1, 1)
    plt.title('t='+str(np.around(tm,4))+'s')
    plt.plot(x, rho, 'k-')
    plt.ylabel('$rho$',fontsize=16)
    plt.tick_params(axis='x',bottom=False,labelbottom=False)
    plt.grid(True)

    plt.subplot(4, 1, 2)
    plt.plot(x, u, 'r-')
    plt.ylabel('$u$',fontsize=16)
    plt.tick_params(axis='x',bottom=False,labelbottom=False)
    plt.grid(True)

    plt.subplot(4, 1, 3)
    plt.plot(x, p, 'b-')
    plt.ylabel('$p$',fontsize=16)
    plt.tick_params(axis='x',bottom=False,labelbottom=False)
    plt.grid(True)

    plt.subplot(4, 1, 4)
    plt.plot(x, E, 'g-')
    plt.ylabel('$E$',fontsize=16)
    plt.grid(True)
    plt.xlim(x_ini,x_fin)
    plt.xlabel('x',fontsize=16)
    plt.subplots_adjust(left=0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(top=0.95)
    
    plt.show()
    #fig.savefig("fig_Sod_Roe_it"+str(it)+".pdf", dpi=300)


# Profiles Matrix

def prfls_matrix (profiles,mtr_profiles):
    
    profiles=np.array([profiles.copy()]).T
    mtr_profiles=np.concatenate([mtr_profiles,profiles],axis=1)
    
    return (mtr_profiles)


# Trend Graphic Generator

def trnd_graph_gen (TimeRng,Mprofile,nclls,plt_trnd,g_title,x_label,y_label):
    
    plt.plot(TimeRng, Mprofile[0,:],label='cell #0')
    plt.plot(TimeRng, Mprofile[plt_trnd,:],label='cell #'+str(plt_trnd))
    plt.plot(TimeRng, Mprofile[nclls,:],label='cell #'+str(nclls))
    plt.title(g_title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.show()





