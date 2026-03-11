# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 23:24:02 2026

@author: sbanerj1
"""
# %reset -f  : run in the console to clear all variables

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(50)


#PARAMETERS
# Population dynamics
# scaled parameter

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


T=20000
N1=20000 
dt=T/N1
dW1= np.sqrt(dt)*np.random.randn(N1)
dW2= np.sqrt(dt)*np.random.randn(N1)
#W= cumsum(dW)   
x0=0.5
y0=3

sigma1=0.04 
sigma2=0.04


R=1; Dt= R*dt; L=N1/R;                          # L EM steps of size Dt= R*dt

L=int(L)
xem= np.zeros(L)                               # preallocate for efficiency
P= x0

yem= np.zeros(L)                               # preallocate for efficiency
Z=y0

for j in range(0,L):
    
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
    
    Winc1=np.sum(dW1[R*j:R*(j+1)])
    P=P+Dt*Cu_r*r*P*(1-P/K)-Dt*Cu_a*a*P*Z/(kP+P)+Dt*i*(K-P)+sigma1*Winc1*P
    
    #if(x<1e-20)
    #    x=0;

    
    xem[j]=P
   
    Winc2=np.sum(dW2[R*j:R*(j+1)])
    Z=Z+Dt*Cu_a*chi*a*P*Z/(kP+P)-Dt*Cu_d*d*Z-Dt*Cu_m*F*Z**2/(Z**2+kZ**2)+sigma2*Winc2*Z
   
    #  if(y<0)
    #     y=0;

    
    yem[j]=Z
   
t=np.arange(0,T+Dt,Dt)

plt.figure(1)

plt.plot(t,np.concatenate(([x0], xem)))
plt.plot(t,np.concatenate(([y0], yem)))
plt.xlim(0, 10000)
plt.title("Copper-plankton stochastic model")
plt.xlabel("Time")
plt.ylabel("Population density")

plt.figure(2)
plt.plot(np.concatenate(([x0],xem)),np.concatenate(([y0],yem)))
plt.title("Copper-plankton stochastic model")
plt.xlabel("Phytoplankton density")
plt.ylabel("Zooplankton density")