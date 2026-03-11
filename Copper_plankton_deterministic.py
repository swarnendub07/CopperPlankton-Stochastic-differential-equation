# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 15:54:14 2026

@author: sbanerj1
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# time interval
T=1000
t_span=(0,T)



def copper_pankton_fish(t,y):
    P=y[0]
    Z=y[1]

# parameters 

    K=2.00              #   (0.1-5) Scenedesmus carrying capacity
    r=0.5               #   Scenedesmus intrinsic rate of natural increase
    a=0.4               #   max intake rate of Daphnia                             %1.8
    kP=0.6              #   half saturation constant of Daphnia                    %0.8
    chi=0.6             #   Daphnia conversion efficiency
    d=0.05              #   Daphnia natural mortality rate
    F=0.1               #   Fish predation rate of X
    kZ=0.5
    i=0.03

# Copper internal concentration 
# Cu=6-10 %% disease free equilibrium oscillation

    Cu=14.3             #  external copper concentration
    kmP=20              #  Scenedesmus maximum intake rate                    (to be estimated)
    kmZ=15              #  Daphnia maximal intake rate
    kcP=6               #  Scenedesmus half saturation constant               (to be estimated)
    kcZ=7               #  Daphnia half saturation constant                   (to be estimated)
    keP=1               #  Scenedesmus constant loss rate                     (to be estimated)
    keZ=1               #  Daphnia constant loss rate                          


# Copper effects

    v_r=4               #  Scenedesmus growth's deficiency EC50
    u_r=50              #  Scenedesmus growth's toxicity EC50
    d_r=5               #  Copper effect on Scenedesmus growth
    b_r=2               #  Copper effect on Scenedesmus growth

    v_a=5               #  Daphnia predation's deficiency EC50
    u_a=16.8            #  Daphnia predations toxicity EC50
    d_a=5               #  Copper effect on Daphnia predation
    b_a=1               #  Copper effect on Daphnia predation

    p_m=0.021           #  Copper response coefficient for Daphnia mortality


    v_m=5
    u_m=16.8
    d_m=5
    b_m=1 
    
    Cp_Cu=Cu*kmP/(Cu+kcP)*1/keP                                             
    Cz_Cu=((Cu*kmZ)/(Cu+kcZ)+chi*a*P/(kP+P)*Cp_Cu)*1/keZ 

    Cp_vr=v_r*kmP/(v_r+kcP)*1/keP                        
    Cp_ur=u_r*kmP/(u_r+kcP)*1/keP 

    Cp_va=v_a*kmP/(v_a+kcP)*1/keP                                              
    Cp_ua=u_a*kmP/(u_a+kcP)*1/keP 

    Cp_vm=v_m*kmP/(v_m+kcP)*1/keP                                              
    Cp_um=u_m*kmP/(u_m+kcP)*1/keP                                              

    Cz_va=((v_a*kmZ)/(v_a+kcZ)+chi*a*P/(kP+P)*Cp_va)*1/keZ                        
    Cz_ua=((u_a*kmZ)/(u_a+kcZ)+chi*a*P/(kP+P)*Cp_ua)*1/keZ     

    Cz_vm=(v_m*kmZ/(v_m+kcZ)+chi*a*P/(kP+P)*Cp_vm)*1/keZ                       
    Cz_um=(u_m*kmZ/(u_m+kcZ)+chi*a*P/(kP+P)*Cp_um)*1/keZ 
    
    Cu_r=-1+np.tanh(d_r*(Cp_Cu-Cp_vr))-np.tanh(b_r*(Cp_Cu-Cp_ur))                    
    Cu_a=0.5*np.tanh(d_a*(Cz_Cu-Cz_va))-0.5*np.tanh(b_a*(Cz_Cu-Cz_ua))                  
    Cu_d=1+p_m*Cz_Cu
    Cu_m=0.5*np.tanh(d_m*(Cz_Cu-Cz_vm))-0.5*np.tanh(b_m*(Cz_Cu-Cz_um))

    dPdt=Cu_r*r*P*(1-P/K)-Cu_a*a*P*Z/(kP+P)+i*(K-P)
    dZdt=Cu_a*chi*a*P*Z/(kP+P)-Cu_d*d*Z-Cu_m*F*Z**2/(Z**2+kZ**2)

    return [dPdt, dZdt]

# initial conditions
y01=[0.5, 3]
y02=[1, 0.5]


sol1=solve_ivp(copper_pankton_fish, t_span, y01)
sol2=solve_ivp(copper_pankton_fish, t_span, y02)

# extract results

t1=sol1.t
P1=sol1.y[0]
Z1=sol1.y[1]

t2=sol2.t
P2=sol2.y[0]
Z2=sol2.y[1]

plt.figure(1)
plt.plot(t1, P1, color='darkgreen', label="Phytoplankton (IC1)")
plt.plot(t1, Z1, color='darkred', label="Zooplankton (IC1)")
plt.plot(t2, P2, color='limegreen', label="Phytoplankton (IC2)")
plt.plot(t2, Z2, color='red', label="Zooplankton (IC2)")

plt.xlabel("Time")
plt.ylabel("Population density")
plt.title("Copper-plankton deterministic model")
plt.legend()


