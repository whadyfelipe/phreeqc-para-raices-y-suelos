#suelos eq richards no lineal con 1 sola especie quimica con raices y succion de agua
#con succion de nutrientes
#ec de concentracion lineal con el termino F no lineal
#modelo de raices muy simple aproximado 

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 11:10:54 2020

@author: whadyimac
"""
import os
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.optimize import newton
from numpy.linalg import solve
from scipy import interpolate
import matplotlib.pyplot as plt


D0=1e3 #cm2/day
K0=51 #cm
r0=2
r1=6.4
l00=15.7
#l01=2.7
K1=8 #cm paper
#Ks=5 #cm/day  FUE CONVERTIDA EN FUNCION DE z 
Km=5.8e-3 #micro mol/cm3
fi=0.4 #porosity
b=239 #buffer for phosfate no units
d=2 #nutirent impedance phosfate no units
Df=1e-5 #nutrient diffusivity phosphate cm2/day
m=0.5 #non-saturated water flow paramter
W=0.05
Wdim=D0/K0*W
rho=0.1
rhodim=rho*(D0*Km)/K0
Pc=0.232e5 #Pa
a0=5e-2 #cm root radio
a=a0
a1=a0 #asumida
L1=8 #cm
beta=np.arccos(0.5)
kr=7.85e-10*1e4*86400*1e-6/(2*np.pi*a) #cm^4*dia^-1*Pa^-1
kz0=1.875e-10*1e8*86400*1e-6 #cm^4*dia^-1*Pa^-1
kz1=4.6e-14*1e8*86400*1e-6 #cm^4*dia^-1*Pa^-1
Pz0=-1e6 #Pa pressure at the surface
Fm=3.26e-6*86400 #micro mol*cm^-2*day^-1
L=50 #cm
L0=L
lb=l00 #cm asumido
ln=0.7
la=15 #cm
zmax=100 #cm total simulation depth
ns=7 #number of chemical species to be transported
n=200 #nodes
tfinal=30 #days
dx=zmax/(n-1) #cm
dt=0.5 #days
ntsteps=int(tfinal/dt) #time steps

zcoo=np.linspace(0,zmax,n)
times=np.linspace(0,tfinal,ntsteps)      
nodes=list(range(n))
set1=[0]
set2=list(range(1,n-1))
set3=[n-1]

def Ks(z):
#TABLA PARA INTERPOLAR LOS VALORES DE Ks VARIABLES CON LA PROFUNDIDAD    
    x=np.array([0,20,40,60,80,100])
    Ks_values=np.array([5,4,3,2,3,4])
    Ks_interp=interpolate.interp1d(x,Ks_values)
    y=Ks_interp(z)
    return y
    
S=np.zeros(n)
def P(S):
    if S>0 and S<1:
        y=-Pc*(S**(-1/m)-1)**(1-m)
    if S<=0:
        Saux=1e-4
        yp=(1-m)*(-1/m)*Saux**(-1/m-1)*(Saux**(-1/m)-1)
        y=(S-Saux)*yp+P(Saux)
    if S>=1:
        y=0
    return y

def psi1(z):
    if z>=lb and z<=L0-la:
        y=1/ln
    else:
        y=0
    return y

def D(S):
    if S>0 and S<1:
        y=S**(1/2-1/m)*((1-S**(1/m))**(-m)+(1-S**(1/m))**m-2)
    if S>=1:
        d=0.999999
        y=dD(d)*(S-d)+D(d)
    if S<=0:
        y=0
    return y

def dD(S):
    if S>0 and S<1:
        y=S**(1/2-1/m)*(S**(-1+1/m)*(1-S**(1/m))**(-1-m)-S**(-1+1/m)*(1-S**(1/m))**(-1+m))+\
        (1/2-1/m)*S**(-1/2-1/m)*(-2+(1-S**(1/m))**(-m)+(1-S**(1/m))**m)
    if S<=0:
        y=0
    if S>=1:
        xe=0.999999
        d=1e-6
        yp=(dD(xe)-dD(xe-d))/d
        y=dD(xe)+yp*(S-xe)
    return y
    

def k(S):
    if S>0 and S<=1:
        y=S**(1/2)*(1-(1-S**(1/m))**m)**2
    if S<=0:
        y=0
    if S>1: y=1    
    return y

def dk(S):
    if S>0 and S<1:
        y=2*S**(-1/2+1/m)*(1-S**(1/m))**(-1+m)*(1-(1-S**(1/m))**m)+\
        (1-(1-S**(1/m))**m)**2/(2*S**(1/2))
    if S<=0:
        y=0
    if S>=1:
        xe=0.999999
        d=1e-6
        yp=(dk(xe)-dk(xe-d))/d
        y=dk(xe)+yp*(S-xe)
        
    return y    

S=np.zeros(n)
S0=np.zeros(n)   #initial condtion for S
X0=np.zeros(n) #initialguess for optimization for S
bnds=np.zeros((n,2))
Smat=np.zeros((n,ntsteps+1))
cmat=np.zeros((n,ntsteps+1))
u=np.zeros(n)
Pr=np.zeros(n)

def Fw(S):
#solution of root pressure equation
    n0=100
    APr=np.zeros((n0,n0))
    BPr=np.zeros(n0)
    dx0=L0/(n0-1)
    z0=np.linspace(0,L0,n0)
    set10=[0]
    set20=list(range(1,n0-1))
    set30=[n0-1]    
    F=np.zeros(n0)
    interpS=interpolate.interp1d(zcoo,S)
    S0=interpS(z0)
    for i in range(n0):
        if i in set20:
            APr[i,i+1]=-kz0*(1/dx0**2)
            APr[i,i]=-kz0*(-2/dx0**2)+\
            (2*np.pi*a*kr+(2*np.pi*a*kr*kz1)**(1/2)*psi1(z0[i]))
            APr[i,i-1]=-kz0*(1/dx0**2)
            BPr[i]=(2*np.pi*a*kr+(2*np.pi*a*kr*kz1)**(1/2)*psi1(z0[i]))*P(S0[i])
        if i in set10:
            APr[i,i]=1
            BPr[i]=Pz0
        if i in set30:
            APr[i,i]=1/dx0
            APr[i,i-1]=-1/dx0
            BPr[i]=0
    Pr=solve(APr,BPr)
    for i in range(n0):
        if S0[i]>0.01:
            F[i]=(2*np.pi*a*kr+(2*np.pi*a*kr*kz1)**(1/2)*psi1(z0[i]))/\
            (np.pi*(a0+K1*np.cos(beta))**2)*(P(S0[i])-Pr[i]) 
        else:
            F[i]=0
    return z0,F

def Fc0(S,c,z,t):
    if z<=L0:
        lambda0=Fm*a0/(Df*fi**(1+d)*S**(1+d)*Km)
        alfa0=4*np.exp(-np.euler_gamma)*Df*fi**(1+d)*S**(1+d)/(a0**2*(b+fi*S))
#        Lc0=lambda0/2*np.log(alfa0*t+alfa0*K0/r0*np.log(1-z/K0)+1)
        Lc0=lambda0/2*np.log(alfa0*t+1)
        F0=2*np.pi*a0*(2*Fm*c/Km)/(1+c/Km+Lc0+np.sqrt(4*c/Km+(1-c/Km+Lc0)**2))
    else:
        F0=0
    return F0

def fc1(S,c,z,t):
    if z>lb and z<=L0:
        tau=-(K0/r0)*np.log(1-lb/K0)
        lambda1=Fm*a1/(Df*fi**(1+d)*S**(1+d)*Km)
        alfa1=4*np.exp(-np.euler_gamma)*Df*fi**(1+d)*S**(1+d)/(a1**2*(b+fi*S))        
#        Lc1=lambda1/2*np.log(alfa1*t+alfa1*K0/r0*np.log(1-z/K0)+\
#        alfa1*K1/r1*np.log(1-(z-zprime)/(K1*np.cos(beta)))+1)
        Lc1=lambda1/2*np.log(alfa1*(t-tau)+1)
        f1=2*np.pi*a1/np.cos(beta)*(2*Fm*c/Km)/(1+c/Km+Lc1+np.sqrt(4*c/Km+(1-c/Km+Lc1)**2))
    else:
        f1=0
    return f1

#definition of nonlinear RIchards equation system
def Frichards(S):
    f=np.zeros(n)
    Fw_table=Fw(S)
    Fw_node=np.zeros(n)
    Fw_spline=interpolate.interp1d(Fw_table[0],Fw_table[1])
    for i in range(n):
        if zcoo[i]<=L0:
            Fw_node[i]=Fw_spline(zcoo[i])
        else:
            Fw_node[i]=0
    for i in range(n):
        if i in set2:
            dS=(S[i+1]-S[i-1])/(2*dx)
            mu=D0*dD(S[i])*dS-Ks(zcoo[i])*dk(S[i])
            if -mu>=0: #simple upwind
                dS=(S[i]-S[i-1])/dx
#                mu=D0*dD(S[i])*dS-Ks*dk(S[i])
            else:
                dS=(S[i+1]-S[i])/dx
#                mu=D0*dD(S[i])*dS-Ks*dk(S[i])
            f[i]=fi/dt*(S[i]-S0[i])-D0*D(S[i])*(S[i+1]-2*S[i]+S[i-1])/dx**2-\
            mu*dS+Fw_node[i]+k(S[i])*(Ks(zcoo[i+1])-Ks(zcoo[i-1]))/(2*dx)
        if i in set1:
            f[i]=-D0*D(S[i])/dx*(S[i+1]-S[i])+Ks(zcoo[i])*k(S[i])-Wdim
        if i in set3:
            f[i]=-D0*D(S[i])/dx*(S[i]-S[i-1])+Ks(zcoo[i])*k(S[i])
#    fmag=np.linalg.norm(f)**2        
    return f

#definition of nonlinear concentration equation system
def Fconcent(c,*args):
    f=np.zeros(n)
    Fw_node=np.zeros(n)
    Fw_spline=interpolate.interp1d(Fw_table[0],Fw_table[1])
    for i in range(n):
        if zcoo[i]<=L0:
            Fw_node[i]=Fw_spline(zcoo[i])
        else:
            Fw_node[i]=0
    for i in range(n):
        if i in set2:
            Fctot=(Fc0(S[i],c[i],zcoo[i],t)+\
            psi1(zcoo[i])*fc1(S[i],c[i],zcoo[i],t))/(np.pi*(a0+K1*np.cos(beta))**2)    
#            Fctot=Fc0(S[i],c[i],zcoo[i],t)/(np.pi*a0**2)
            dS=(S[i+1]-S[i-1])/(2*dx)
            mu=D0*D(S[i])*dS+Df*fi**(d+1)*(d+1)*S[i]**d*dS-Ks(zcoo[i])*k(S[i])
#            u=-D0*D(S[i])*dS+Ks*k(S[i])
            if -mu>=0: #upwind
                dc=(c[i]-c[i-1])/dx
            else:
                dc=(c[i+1]-c[i])/dx
#eqn conveccion difusion siplificada
#            dc=(c[i+1]-c[i-1])/(2*dx)           
#            f[i]=(c[i]-c0[i])/dt-\
#            Df*(c[i+1]-2*c[i]+c[i-1])/dx**2-mu*dc-c[i]*Fw_node[i]+Fctot
#eqn convection del paper con alto buffer
#            dc=(c[i+1]-c[i-1])/(2*dx)          
#            f[i]=(b+fi*S[i])*(c[i]-c0[i])/dt-\
#            Df*fi**(d+1)*S[i]**(d+1)*(c[i+1]-2*c[i]+c[i-1])/dx**2-mu*dc\
#            -c[i]*Fw_node[i]+Fctot
#eqn convection del paper con alto buffer simplificadacomo dicen en el paper
#            dc=(c[i+1]-c[i])/(dx)            
            f[i]=(fi*S[i])*(c[i]-c0[i])/dt-\
            Df*fi**(d+1)*S[i]**(d+1)*(c[i+1]-2*c[i]+c[i-1])/dx**2-mu*dc\
                -c[i]*Fw_node[i]+Fctot
        if i in set1:
#collocation condition
            f[i]=-Df*fi**(1+d)*S[i]**(1+d)*(c[i+1]-c[i])/dx+\
            Wdim*c[i]-rhodim
#dirichlet condition
#            f[i]=c[i]-0.66*Km
#finite volume 
#            Fctot=(Fc0(S[i],c[i],zcoo[i],t)+\
#            psi1(zcoo[i])*fc1(S[i],c[i],zcoo[i],t)*dx/2)/(np.pi*(a0+K1*np.cos(beta))**2)    
#            dS=(S[i+1]-S[i])/dx
#            mu=D0*D(S[i])*dS+Df*fi**(d+1)*(d+1)*S[i]**d*dS-Ks*k(S[i])
#            f[i]=(b+fi*S[i])*(c[i]-c0[i])/dt*dx/2+\
#            (c[i]*Wdim-rhodim)-Df*fi**(d+1)*S[i]**(d+1)*(c[i+1]-c[i])/dx\
#            -c[i]*Fw_node[i]*dx/2+Fctot*dx/2   
        if i in set3:
            f[i]=(c[i]-c[i-1])/dx
    return f

#solution of richards equations system
S0[:]=0.75
X0[:]=0.5 #this gugess is important it may not converge
Smat[:,0]=S0
for i in range(n):
    u[i]=-D0*D(S0[i])*(S0[i]-S0[i-1])/dx+Ks(zcoo[i])*k(S0[i])


c0=np.zeros(n)
c=np.zeros(n)

c0[:]=0.66*Km
cmat[:,0]=c0

for j in range(ntsteps):
    t=j*dt
    print('Time= ', t) 
    L0=l00+K0*(1-np.exp(-r0*t/K0))
    tau=-(K0/r0)*np.log(1-lb/K0)
    if t>=tau:
        L1=K1*(1-np.exp(-r1*(t-tau)/K1))
#    for i in range(n):
#        bnds[i]=(1e-6,0.999999)
#solve as optimization 
#    res = minimize(Frichards, X0, method='TNC', bounds=bnds)            
#    S=res.x
#    print('Optimization=',res.success)            
#solve with generic fsolve with options
    res=fsolve(Frichards,X0,xtol=1e-4,full_output=1,col_deriv=1,\
    band=(1,1),epsfcn=1e-6)
    print('Solution code for Richards equation',res[2],res[3])
    F_residual=np.linalg.norm(res[1]['fvec'])
    print('fvec_norm=',F_residual)
    S=res[0]
#solve with root for comparison
#    res=root(Frichards, X0, method='lm',options={"xtol":1e-2,"ftol":1e-2})
#    print('Solution code for Richards equation',res.success,res.status,res.message)
#    S=res.x
#edit negative values after solving
    for i in range(n):
        if S[i]<=0.01 : S[i]=0.01  
    S0=S
    X0=S0
    Smat[:,j+1]=S
#recover the velocity values just in case they're needed
    for i in range(n):
        u[i]=-D0*D(S[i])*(S[i]-S[i-1])/dx+Ks(zcoo[i])*k(S[i])
#Recover the Fw values at the root nodes        
    Fw_table=Fw(S)
#solution of concentration equations for k species
    res=fsolve(Fconcent,c0,args=(t),full_output=1,col_deriv=1,band=(1,1))
    print('Solution code for concentration equation',res[2],res[3])
    print('=======================================================')
    c=res[0]
    for i in range(n):
        if c[i]<=0 : c[i]=0  #edit negative values
    c0=c
    cmat[:,j+1]=c
        
#save results
if not os.path.exists('resultados_richards'):
    os.mkdir('resultados_richards')

#np.savetxt('resultados_richards/c.txt',cmat)        
np.savetxt('resultados_richards/S.txt',Smat)        
np.savetxt('resultados_richards/t.txt',times)        
    
plt.figure()  
for i in range(0,ntsteps,1):
    plt.plot(zcoo,cmat[:,i]/Km)
plt.ylabel('c')
plt.xlabel('y')        
    

plt.figure()  
for i in range(0,ntsteps,1):
    plt.plot(zcoo,Smat[:,i])
plt.ylabel('S')
plt.xlabel('y')   

