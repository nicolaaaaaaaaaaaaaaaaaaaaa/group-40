import numpy as np
import math as m
import pandas as pd
import os
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def update_a(a_old,ac,CT,F,phi,solidity,Cn):
    if a_old<=ac :
            a_new= CT/(4*F*(1-a_old))
    else :
        if ac==1/3 :
            a_et = a_et_calc(CT,F,a_old)
            a_new = 0.1*a_et + (1-0.1)*a_old
        if ac==0.2 :
            K = K_calc(F,phi,solidity,Cn)
            a_new = 1 + K/2 *(1- 2*ac) - 1/2 * np.sqrt((K*(1-2*ac)+2)**2+4*(K*ac**2-1))
    return(a_new)

def update_a_prime(a_prime_old,F,phi,solidity,Ct):       
    return(Ct * solidity * (1+a_prime_old) / (4*F*np.sin(phi)*np.cos(phi)))

def K_calc(F,phi,solidity,Cn):
     K = (4*F*np.sin(phi)**2)/(solidity*Cn)
     return(K)

def a_et_calc(CT,F,a_old):
     a_et = CT / (4*F*(1-1/4 * (5 - 3 *a_old)*a_old))
     return(a_et)




def TSR_calc(V0,R,omega):
    TSR = omega * R / V0
    return(TSR)

def V0_calc(TSR,R,omega):
    V0 = omega * R / TSR
    return(V0)

def omega_calc(TSR,V0,R):
    omega = TSR * V0 / R
    return(omega)

def R_calc(TSR,V0,omega):
    R = TSR * V0 / omega
    return(R)



def solidity_calc(c,B,r):
    solidity=c*B/(2*np.pi*r)
    return(solidity)



def Vrel_calc(V0,a,a_prime,omega,r):
    Vrel=np.sqrt((V0-a*V0)**2+((a_prime+1)*omega*r)**2)
    return(Vrel)

def pn_calc(rho,Vrel,c,Cn):
    pn=0.5*rho*Vrel**2*c*Cn
    return(pn)

def pt_calc(rho,Vrel,c,Ct):
    pt=0.5*rho*Vrel**2*c*Ct
    return(pt)

def Cn_calc(Cl,Cd,phi):
    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
    return(Cn)

def Ct_calc(Cl,Cd,phi):
    Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
    return(Ct)

def CT_calc(a=None,solidity=None,phi=None,Cn=None,rho=None,V0=None,R=None,T=None,Tmax=None):
    """
    a,solidity,phi,Cn OR T,Tmax OR T,rho,V0,R
    """
    if not a==None:
        CT=(1-a)**2*solidity*Cn/np.sin(phi)**2
    elif not Tmax==None:
        CT=T/Tmax
    else:
        CT=T/Tmax_calc(rho,V0,R)
    return(CT)

def CP_calc(P,rho,V0,R,Pmax=None):
    """
    P,rho,V0,R OR P,Pmax
    """
    if Pmax==None:
        CP=P/Pmax_calc(rho,V0,R)
    else :
        CP=P/Pmax
    return(CP)

def Tmax_calc(rho,V0,R):
    Tmax=0.5*rho*V0**2*np.pi*R**2
    return(Tmax)

def Pmax_calc(rho,V0,R):
    Pmax=0.5*rho*V0**3*np.pi*R**2
    return(Pmax)



def F_calc(B,R,r,phi):  
    F = 2/np.pi * np.arccos(np.exp(-B/2 * (R-r)/(r * np.sin(abs(phi)))))
    return(F)




def deg(angle):
    return(angle*180/np.pi)

def rad(angle):
    return(angle*np.pi/180)

def rpm(omega):
    N=30*omega/np.pi
    return(N)

def rad_per_sec(N):
    omega=np.pi*N/30
    return(omega)

def phi_calc(a=None,a_prime=None,R=None,TSR=None,r=None,alpha=None,theta=None):
    """
    a,R,TSR,r OR alpha,theta
    """
    if not a==None :
        phi=np.arctan(((1-a)*R)/((1+a_prime)*TSR*r))
    else :
        phi = alpha + theta
    return(phi)

def beta_calc(theta,theta_p):
    beta = theta - theta_p
    return(beta)

def theta_calc(phi=None,alpha=None,beta=None,theta_p=None):
    """
    phi,alpha OR beta,theta_p
    """
    if not alpha==None and not phi==None:
        theta = phi - alpha
    elif not beta==None and not theta_p==None:
        theta = beta + theta_p
    return(theta)

def alpha_calc(phi,theta):
    alpha= phi - theta
    return(alpha)

def theta_p_calc(theta,beta):
    theta_p = theta - beta
    return(theta_p)

def dCP_calc(TSR,B,a,Ct,phi,c,R):
    dCP = TSR*B*(1-a)**2*Ct / (2*np.pi*np.sin(phi)**2) * c/R
    return(dCP)

def get_bladedata(value,type,blade=pd.DataFrame()):
    """
    type = 'r', 'c' or 'beta'
    """
    
    if blade.empty :
        blade=pd.read_csv('bladedat.txt',delimiter='\t',header=None,names=['r','beta','c','t_perc'])
        
    r=      np.interp(value,blade[type],blade['r'])
    c=      np.interp(value,blade[type],blade['c'])
    beta=   np.interp(value,blade[type],blade['beta'])
    t_perc= np.interp(value,blade[type],blade['t_perc'])
        
    return(r,c,beta,t_perc)

