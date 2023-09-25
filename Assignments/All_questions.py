""" Ashes questions : cgrinde@dtu.dk"""
""" 7 to 9 pages"""

""" Importation of libraries """

import numpy as np
import math as m
import pandas as pd
import os
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from compute import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ac=1/3
# ac=0.2
eps=1e-10
max_iteration=10000
B=3
rho=1.225
R=89.17

blade=pd.read_csv('bladedat.txt',delimiter='\t',header=None,names=['r','beta','c','t_perc'])
    
list_radius=list(blade['r'])


#os.chdir('Assignment_data')

files=['FFA-W3-241.txt','FFA-W3-301.txt','FFA-W3-360.txt','FFA-W3-480.txt','FFA-W3-600.txt','cylinder.txt']

#Initializing tables    
cl_tab=np.zeros([105,6])
cd_tab=np.zeros([105,6])
cm_tab=np.zeros([105,6])
aoa_tab=np.zeros([105,])
#Readin of tables. Only do this once at startup of simulation
for i in range(np.size(files)):
    aoa_tab[:],cl_tab[:,i],cd_tab[:,i],cm_tab[:,i] = np.loadtxt(files[i], skiprows=0).T

# Thickness of the airfoils considered
# NOTE THAT IN PYTHON THE INTERPOLATION REQUIRES THAT THE VALUES INCREASE IN THE VECTOR!

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;

def force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])
    

    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cm_tab[:,i])
    
    #Interpolate to current thickness:
    cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    cm=np.interp (thick,thick_prof,cm_aoa[0,:])


    return cl, cd, cm 


"""
PLOT FUNCTIONS
"""

def plot_contour(fig,x,y,z,name=None,colorscale=None,coloring ='heatmap', showlabels = True, labelfont_size=14, labelfont_color='black'):
    fig.add_trace(
        go.Contour(
            z=z,
            x=x, # horizontal axis
            y=y,
            name=name,
            colorscale=colorscale,
            contours=dict(
            coloring =coloring,
            showlabels = showlabels,
            labelfont = dict(size = labelfont_size,color = labelfont_color)
            )))
    
def plot_point(fig, x, y, name=None,text=None,marker_color='black',marker_size=10,text_position='top center',secondary_y=False):
    fig.add_trace(
        go.Scatter(
            mode='markers+text',
            x=[x], # horizontal axis
            y=[y], # vertical axis
            text=text,
            textposition=text_position,
            marker=dict(color=marker_color, size=marker_size),
            name=name,
            showlegend=False
        ),secondary_y=secondary_y)

def plot_line(fig,x,y,name=None,line_color=None,secondary_y=False,dash=None,mode=None):
    fig.add_trace(
        go.Scatter(
            x=x, # horizontal axis
            y=y, # vertical axis
            line_color = line_color,
            line=dict(dash=dash),
            name=name,
            mode=mode
            ), secondary_y=secondary_y)

def plot_layout(fig,x_title,y_title,graph_title,secondary_y_title=None,secondary_y=False,gridcolor='gray',width=800,height=600):
    fig.update_layout(title_text=graph_title,plot_bgcolor='rgba(0, 0, 0, 0)',width=1000,height=600)
    fig.update_xaxes(title_text=x_title,gridcolor=gridcolor,zerolinecolor=gridcolor)
    fig.update_yaxes(title_text=y_title,gridcolor=gridcolor,zerolinecolor=gridcolor)
    if secondary_y:
        fig.update_yaxes(title_text=secondary_y_title,gridcolor=gridcolor,zerolinecolor=gridcolor,secondary_y=True)

def display_results(name_file,compare_file=None,
                    x_names=None,y_names=None,curve_names=None,
                    compare_name='Base',
                    x_titles=None,y_titles=None,graph_titles=None,
                    colors=None):

    Res = pd.read_excel(name_file)
    if compare_file : Res_Comp = pd.read_excel(compare_file)
    Ashes = pd.read_excel('Ashes_global_values.xlsx')
    for i in range(len(y_names)):
        print(colors[i])
        fig = make_subplots()
        for j in range(len(y_names[i])):
            plot_line(fig,x=Res[x_names[i]],y=Res[y_names[i][j]],name=curve_names[i][j],line_color=colors[i][j])
            if compare_file : plot_line(fig,x=Res_Comp[x_names[i]],y=Res_Comp[y_names[i][j]],name=curve_names[i][j]+" "+compare_name,dash='dot',line_color=colors[i][j])
            if y_names[i][j] in Ashes.columns :
                plot_line(fig,x=Ashes[x_names[i]],y=Ashes[y_names[i][j]],name=curve_names[i][j]+" Ashes",mode='markers',line_color='black')
        plot_layout(fig,x_title=x_titles[i],y_title=y_titles[i],graph_title=graph_titles[i])
        fig.show()


"""
BEM and load calculations
"""

def bem_algorithm(r, R, solidity, TSR, B, beta=None, alpha=None, theta_p=None, Cl=None, Cd=None, t_perc=None, ac=1/3, eps=1e-10, max_iteration=10000):
    
    alpha_b=alpha
    a_new, a_prime_new = 0, 0
    k, count = 0, 0

    while k==0:

        count+=1      
        a_old, a_prime_old = a_new, a_prime_new
        
        phi = phi_calc(a=a_old,a_prime=a_prime_old,R=R,TSR=TSR,r=r)
        
        if not t_perc==None:
            if alpha_b==None:
                theta = theta_calc(phi=phi,beta=beta,theta_p=theta_p)
                alpha = alpha_calc(phi=phi,theta=theta)
            Cl,Cd,Cm = force_coeffs_10MW(alpha*180/np.pi,t_perc,aoa_tab,cl_tab,cd_tab,cm_tab)
        
        Cn = Cn_calc(Cl,Cd,phi)
        Ct = Ct_calc(Cl,Cd,phi)
        
        F = F_calc(B,R,r,phi)

        CT= CT_calc(a=a_old,solidity=solidity,phi=phi,Cn=Cn)

        a_new = update_a(a_old,ac,CT,F,phi,solidity,Cn)
        a_prime_new = update_a_prime(a_prime_old,F,phi,solidity,Ct)
        
        if abs(a_new - a_old) < eps and abs(a_prime_new - a_prime_old) < eps:
            k=1
        if count>max_iteration :
            k=1
            print('not converged')
    
    BEM=pd.DataFrame([[a_new,a_prime_new,Cn,Ct,Cl,Cd,CT,F,phi,alpha,count]],columns=['a','a_prime','Cn','Ct','Cl','Cd','CT','F','phi','alpha','iterations'])
    return(BEM)

def get_local_loads(r,B,V0,omega,rho,R,TSR,theta_p,solidity,beta,t_perc,c,ac,eps,max_iteration):
    
    BEM=bem_algorithm(r=r,R=R,
                    solidity=solidity, 
                    TSR=TSR, B=B,
                    beta=beta, theta_p=theta_p, 
                    t_perc=t_perc,
                    ac=ac, eps=eps, max_iteration=max_iteration)
    
    a,a_prime,Cn,Ct=BEM.at[0,'a'],BEM.at[0,'a_prime'],BEM.at[0,'Cn'],BEM.at[0,'Ct']
    
    Vrel=Vrel_calc(V0,a,a_prime,omega,r)
    
    pn=pn_calc(rho,Vrel,c,Cn)
    pt=pt_calc(rho,Vrel,c,Ct)
    
    P=omega*B*pt*r
    T=B*pn

    return(pn,pt,T,P)

def get_local_loads_over_radius(list_radius,blade,B,V0,omega,rho,R,TSR,theta_p,ac,eps,max_iteration):
    list_res=[]
    
    for r in list_radius[:-1] :
        
        r,c,beta,t_perc=get_bladedata(value=r,type='r',blade=blade)
        
        beta=rad(beta)
        
        solidity=solidity_calc(c,B,r)
        
        pn,pt,T,P=get_local_loads(r,B,V0,omega,rho,R,TSR,theta_p,solidity,beta,t_perc,c,ac,eps,max_iteration)
        
        list_res.append([r,pn,pt,T,P])

    list_res.append([list_radius[-1],0,0,0,0])
    return(np.array(list_res))

def get_global_loads(list_radius,blade,B,V0,omega,rho,R,TSR,theta_p,ac,eps,max_iteration):
    
    list_res=get_local_loads_over_radius(list_radius,blade,B,V0,omega,rho,R,TSR,theta_p,ac,eps,max_iteration)
    
    P=integrate.simps(y=list_res[:,4],x=list_res[:,0])
    T=integrate.simps(y=list_res[:,3],x=list_res[:,0])
    CP=CP_calc(P,rho,V0,R)
    CT=CT_calc(T=T,rho=rho,V0=V0,R=R)
    
    return(T,P,CT,CP)

def get_global_loads_over_TSR_and_pitch(list_TSR,list_theta_p,list_V0,list_radius,R,B,rho,ac,eps,max_iteration): 
    
    X,Y=np.meshgrid(list_TSR,list_theta_p)
    CP_dist,CT_dist,T_dist,P_dist=np.zeros((len(X),len(X[0]))),np.zeros((len(X),len(X[0]))),np.zeros((len(X),len(X[0]))),np.zeros((len(X),len(X[0])))

    for k1 in range (len(X)):
        
        for k2 in range (len(X[0])):

            TSR=X[k1,k2]
            theta_p=Y[k1,k2]
            omega=omega_calc(TSR,list_V0[k2],R)
            
            T,P,CT,CP=get_global_loads(list_radius,blade,B,list_V0[k2],omega,rho,R,TSR,theta_p,ac,eps,max_iteration)
            
            CP_dist[k1,k2], CT_dist[k1,k2], P_dist[k1,k2], T_dist[k1,k2] = CP, CT, P, T
           
    return(CP_dist,CT_dist,P_dist,T_dist)
                               
def local_cp_over_chord(list_chord,TSR_max,theta_p_max,r,R,B,ac,eps,max_iteration): 

    theta_p_max=rad(theta_p_max)

    r,c,beta,t_perc=get_bladedata(value=r,type='r',blade=blade)
    
    beta=rad(beta)

    solidity=solidity_calc(c,B,r)
    
    BEM=bem_algorithm(r,R,solidity, TSR_max, B, beta=beta, theta_p=theta_p_max, t_perc=t_perc, ac=ac, eps=eps, max_iteration=max_iteration)
    alpha,Cl,Cd=BEM.at[0,'alpha'],BEM.at[0,'Cl'],BEM.at[0,'Cd']
    
    list_dCP=[]

    for c in list_chord :
        
        solidity=solidity_calc(c,B,r)
        
        BEM=bem_algorithm(r,R,solidity, TSR_max, B, theta_p_max, alpha=alpha, Cl=Cl, Cd=Cd, ac=ac, eps=eps, max_iteration=max_iteration)
        
        a,Ct,phi=BEM.at[0,'a'],BEM.at[0,'Ct'],BEM.at[0,'phi']

        dCP=dCP_calc(TSR_max,B,a,Ct,phi,c,R)

        theta=theta_calc(phi=phi,alpha=alpha)
        beta=beta_calc(theta,theta_p_max)

        list_dCP.append([c,beta,dCP])
    
    return(np.array(list_dCP))

def compute_TSR_list(list_V0,V0_rated,TSR_max,omega_max,R):
    list_TSR=[]
    for V0 in list_V0:
        if V0<V0_rated:
            list_TSR.append(TSR_max)
        else : 
            list_TSR.append(TSR_calc(V0,R,omega_max))
    return(np.array(list_TSR))

def loads_as_function_of_V0(list_V0,list_theta_p,CP_dist,CT_dist,P_dist,T_dist,P_max):
    list_result=[]
    for j in range(len(list_V0)):
        V0=list_V0[j]
        CP_V0=0
        for i in range(len(list_theta_p)):
            if P_dist[i,j]<P_max and CP_dist[i,j]>CP_V0:
                CP_V0, CT_V0, P_V0, T_V0 = CP_dist[i,j],CT_dist[i,j], P_dist[i,j], T_dist[i,j]
                theta_p_V0 = list_theta_p[i]
        list_result.append([V0,theta_p_V0,T_V0,P_V0,CT_V0,CP_V0])
    Res=pd.DataFrame(list_result,columns=['Wind speed','Pitch angle','T','P','CT','CP'])
    Res['Pitch angle']*=180/np.pi
    Res.to_excel('Loads_as_function_of_V0.xlsx',index=False)
    return(Res)


"""
Questions
"""

def question_one(list_TSR,list_theta_p,list_radius,V0,R,B,rho,ac,eps,max_iteration):
    
    list_V0=[V0]*len(list_TSR)
    CP_dist,CT_dist,P_dist,T_dist=get_global_loads_over_TSR_and_pitch(list_TSR,rad(list_theta_p),list_V0,list_radius,R,B,rho,ac,eps,max_iteration)
    
    indice_max = np.argmax(CP_dist)
    i, j = np.unravel_index(indice_max, CP_dist.shape)
    Cp_max = CP_dist[i,j]
    TSR_max = list_TSR[j]
    theta_p_max = list_theta_p[i]

    print(f"CP_max ={Cp_max:.3f} \nTSR_max = {TSR_max:.3f} \ntheta_p_max = {theta_p_max:.3f}")

    fig = make_subplots()
    plot_contour(fig,list_TSR,list_theta_p,CP_dist,name=" Power coefficient (CP)")
    plot_point(fig, TSR_max, theta_p_max, text=f"CP_max ={Cp_max:.3f}")
    plot_layout(fig,x_title="Tip speed ratio",y_title="Pitch angle (°)",graph_title="Power coefficient",width=800,height=600)
    fig.show()

def question_two(list_chord,TSR_max, theta_p_max,r,R,B,ac,eps,max_iteration):
    list_cp=local_cp_over_chord(list_chord=list_chord,
                    TSR_max=TSR_max,
                    theta_p_max=theta_p_max,
                    r=r,R=R,
                    B=B,
                    ac=ac,eps=eps,max_iteration=max_iteration)

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    plot_line(fig,list_cp[:,0],list_cp[:,1]*180/np.pi,name="Beta",secondary_y=True)
    plot_line(fig,list_cp[:,0],list_cp[:,2],name="Local CP")
    i=np.argmax(list_cp[:,2])
    plot_point(fig,x=list_cp[i,0],y=list_cp[i,2],text=f"Local CP_max ={list_cp[i,2]:.3f}")
    plot_point(fig,x=list_cp[i,0],y=list_cp[i,1]*180/np.pi,text=f"Beta ={list_cp[i,1]*180/np.pi:.3f}",secondary_y=True)
    plot_layout(fig,x_title="Chord (m)",y_title="Local power coefficient",graph_title="Local CP and Beta as function of chord",secondary_y_title="Beta (°)",secondary_y=True)
    
    fig.show()

def question_3(TSR_max,P_max, CP_max,rho,R):
    V0_rated = np.cbrt(P_max/(0.5*rho*np.pi*R**2*CP_max))
    omega_max = TSR_max*V0_rated/R * 30/np.pi

    print(f"V0_rated ={V0_rated:.3f} \nomega_max = {omega_max:.3f}")

    V0_list=np.arange(4,25,0.1)
    rotational_speed=TSR_max*V0_list/R * 30/np.pi
    for i,omega in enumerate(rotational_speed) :
        if omega>omega_max : rotational_speed[i]=omega_max
    fig = make_subplots()
    plot_line(fig,x=V0_list, y=rotational_speed, name='Rotational speed')
    Report=pd.read_csv('Report_Pitch.csv',sep=' ')
    plot_line(fig,x=Report['Wind_speed[m/s]'], y=Report['RPM'], name='Rotational speed (DTU Report)')
    plot_layout(fig,x_title="Wind Speed (m/s)",y_title="Rotational speed (rpm)",graph_title="Rotational speed")
    fig.show()
    
def question_4(list_V0,list_theta_p,list_radius,B,R,ac,eps,max_iteration,V0_rated,omega_max,P_max,TSR_max):
    
    list_TSR = compute_TSR_list(list_V0,V0_rated,TSR_max,omega_max,R)
    list_theta_p = rad(list_theta_p)    
    CP_dist,CT_dist,P_dist,T_dist=get_global_loads_over_TSR_and_pitch(list_TSR,list_theta_p,list_V0,list_radius,R,B,rho,ac,eps,max_iteration)
    
    Res = loads_as_function_of_V0(list_V0,list_theta_p,CP_dist,CT_dist,P_dist,T_dist,P_max)

def question_6(A,k,Cut_in_speed=4,Cut_out_speed=25,min_graph=0,max_graph=25):

    Res=pd.read_excel('Loads_as_function_of_V0.xlsx')
    wind_speed=list(Res['Wind speed'])
    power=list(Res['P'])
    f=np.zeros(len(wind_speed)-1)
    H=np.zeros(len(wind_speed)-1)
    AEO=np.zeros(len(wind_speed)-1)

    for i in range(len(wind_speed)-1) :
        if wind_speed[i]>=Cut_in_speed and wind_speed[i+1]<=Cut_out_speed:
            f[i]=np.exp(-(wind_speed[i]/A)**k)-np.exp(-(wind_speed[i+1]/A)**k)
            H[i]=8760*f[i]
            AEO[i]=0.5*(power[i]+power[i+1])*H[i]

    return(f,H,np.sum(AEO))

"""
Calcule & Display
"""

""" Q1 """ 
# V0=8
# x=np.arange(5,11,0.5)
# y=np.arange(-4,4,0.5)
# ac=1/3
# question_one(x,y,list_radius,V0,R,B,rho,ac,eps,max_iteration)


""" Q2 """
# r = 71.97
# TSR=8.01
# theta_p=0.09
# list_chord=np.arange(0,3.5,0.01)
# question_two(list_chord=list_chord,
#             TSR_max=TSR,
#             theta_p_max=theta_p,
#             r=r,R=R,
#             B=B,
#             ac=ac,eps=eps,max_iteration=max_iteration)


""" Q3 """
# CP_max =  0.469
# P_max = 10.64e6
# TSR_max=8.01
# question_3(TSR_max,P_max,CP_max,rho,R)

""" Q4 """
# P_max = 10.64e6
# TSR_max=8.01
# V0_rated =11.400 
# omega_max = 9.607
# list_V0=np.arange(4,26,1)
# list_theta_p=np.arange(-1,25.1,0.1)
# question_4(list_V0,list_theta_p,list_radius,B,R,ac,eps,max_iteration,V0_rated,omega_max*np.pi/30,P_max,TSR_max)

# x_names=['Wind speed']*5
# y_names=[['Pitch angle'],['T'],['P'],['CT','CP']]
# curve_names=[['Pitch angle'],['Thrust'],['Power'],['CT','CP']]

# x_titles=["Wind speed V0 (m/s)"]*5
# y_titles=["Pitch angle (°)","Thrust (N)","Power (W)","Coefficient"]
# graph_titles=["Pitch angle",'Thrust','Power','Coefficients']

# colors=[[None],[None],[None],['mediumslateblue','lightseagreen']]

# display_results('Loads_as_function_of_V0.xlsx',compare_file='Report_values.xlsx',
#                 x_names=x_names,y_names=y_names,curve_names=curve_names,
#                 x_titles=x_titles,y_titles=y_titles,graph_titles=graph_titles,
#                 compare_name='Report',colors=colors)


""" Q5 """

# TSR_max=8.01
# V0_rated=11.40
# omega_max=9.82*np.pi/30

# list_V0=np.array([5,9,11,20])
# list_TSR = compute_TSR_list(list_V0=list_V0,V0_rated=V0_rated,TSR_max=TSR_max,omega_max=omega_max,R=R)
# Res=pd.read_excel('Loads_as_function_of_V0.xlsx')
# Ashes=pd.read_excel('Ashes_values.xlsx')
# list_theta_p=list(Res['Pitch angle'].loc[Res['Wind speed'].isin(list_V0)])
# Val=pd.DataFrame()
# for i in range(len(list_V0)):
#     omega=omega_calc(list_TSR[i],list_V0[i],R)
#     list_result=get_local_loads_over_radius(list_radius,blade,B,list_V0[i],omega,rho,R,list_TSR[i],list_theta_p[i]*np.pi/180,ac,eps,max_iteration)
#     Val1=pd.DataFrame(list_result,columns=['Radius','pn','pt','T','P'])
#     Val1['Wind speed']=list_V0[i]
#     Val=pd.concat([Val,Val1])
# Val.to_excel('Loads_distribution_over_radius.xlsx')

# color=['mediumslateblue','lightseagreen','orange','deeppink']
# for load in ['pn','pt']:
#     fig = make_subplots()
#     for i,V0 in enumerate(list(Val['Wind speed'].unique())):
#         plot_line(fig,x=Val['Radius'].unique(),y=Val[load].loc[Val['Wind speed']==V0],name=str(V0)+" m/s",line_color=color[i])
#         plot_line(fig,x=Ashes['Radius'].unique()+2.8,y=Ashes[load].loc[Ashes['Wind speed']==V0],name=str(V0)+" m/s Ashes",dash='dot',line_color=color[i])
#     plot_layout(fig,x_title='Radius (m)',y_title=load+" (N/m)",graph_title=load+" over radius")
#     fig.show()