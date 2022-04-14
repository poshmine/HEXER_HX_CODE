'''
Date of last edit: August 24th, 2020
Author(s): Ryan McGuire*    Lane Carasik^
*Virginia CommonWealth University
*FAST Research Group
HEXER - base HEat EXchangER code
Equations for Reactor Heat Exchangers: 
Python Script for calculation of Shell-and Tube heat transfer coefficient using correctional factors
for Shell-and-Tube bundles based on models found in literature.
'''

'''
Revision Points
A_o_tb
N_rc
'''

## TEST COMMIT 
abefe=1

##Basic Imports and Function Representatives
import numpy as np
import HEXER_variables as v
import matplotlib.pyplot as plt

##Correctional Factor Calculations

#Baffle configuration correctional factor
D_ctl = v.D_otl-v.d_o
theta_ctl = 2*np.arccos((v.D_s-2*v.l_c)/D_ctl)
F_c = 1-(theta_ctl/np.pi)+(np.sin(theta_ctl)/np.pi)
J_c = (0.55+0.72*F_c)
#print(J_c)
#Bundle Leakage effects correctional factor
theta_b = 2*np.arccos(1-(2*v.l_c/v.D_s))
F_w = (theta_ctl/(2*np.pi))-(np.sin(theta_ctl)/(2*np.pi))
d_sb = v.D_s-v.D_baffle
d_tb = v.d_l-v.d_o
A_o_sb = np.pi*v.D_s*(d_sb/2)*(1-(theta_b/(2*np.pi)))
A_o_tb = (np.pi*v.d_o*d_tb*v.N_t*(1-F_w))/2
Case1 = v.p_t/v.d_o
if v.Tbl == 30 or 90:
    A_o_cr = (v.D_s-v.D_otl+2*(D_ctl/v.X_t)*(v.p_t-v.d_o))*v.L_bc
r_s = (A_o_sb)/(A_o_sb+A_o_tb)
r_lm = (A_o_sb+A_o_tb)/(A_o_cr)
J_l = (0.44*(1-r_s)+(1-0.44*(1-r_s))*np.e**(-2.2*r_lm))
#print(J_l)

#Bundle and pass partition bypass correctional factor
G_s = (v.m_s/A_o_cr)
A_o_bp = v.L_bc*(v.D_s-v.D_otl+0.5*v.N_p*v.w_p)
Re_s = round((G_s*v.d_o)/v.u_s)
N_ss_plus = v.N_ss/v.N_rcc
r_b = A_o_bp/A_o_cr
if Re_s<=100:
    C = 1.35
else:
    C = 1.25
if N_ss_plus >= 1/2:
    J_b = 1
else:
    J_b = (np.e**(-C*r_b*(1-(2*N_ss_plus)**(1/3))))
#print(J_b)

#Larger baffle spacing correctional factor
if v.Laminar_flow == 1:
    n = 1/3
else:
    n = 3/5
L_i_plus = v.L_bi/v.L_bc
L_o_plus = v.L_bo/v.L_bc
N_b = round(((v.L-v.L_bi-v.L_bo)/v.L_bc)+1)
J_s = ((N_b-1+L_i_plus**(1-n)+L_o_plus**(1-n))/(N_b-1+L_i_plus+L_o_plus))
#print(J_s)

#Adverse Temperature Gradient Buildup in Laminar Flow correctional factor
N_rcw = (0.8/v.X_t)*(v.l_c-(1/2)*(v.D_s-D_ctl))
N_rc = v.N_rcc+N_rcw
if Re_s >= 100:
    J_r = 1
elif Re_s <= 20:
    J_r = (10/N_rc)**0.18
else:
    J_r = 1+(((10/N_rc)**0.18)-1)/(20-100)*(Re_s-100)
#print(J_r)

##Heat Transfer Calculations - Core

#Perfect Heat Transfer Coefficent 
h_id = (((v.Nu_s*v.k)/v.d_o)*(1)**-0.14)
#print(h_id)

#Calculated Heating Coefficent
h_s = round((h_id*J_c*J_l*J_b*J_s*J_r),3)
h_o = h_s
#print("The Calculated Shell Side Heat Transfer Coefficient is " +str(h_s)+ " W/m^2 * K")

##Shell and Tube Pressure Drop

##Pressure Drop Calculations

#Tube-to-baffle leakage
if Re_s <= 100:
    D = 4.5
else:
    D = 3.7

if N_ss_plus >= 0.5:
    C_b = 1
else:
    C_b = np.e**(-D*r_b*(1-(2*N_ss_plus)**(1/3)))
#print(C_b)

#Bypass Flow
p = -0.15*(1+r_s)+0.8
C_l = np.e**(-1.33*(1+r_s)*r_lm**p)
#print(C_l)

#Differing Baffle Spacing from central section
if v.Laminar_flow == 0:    
    n_prime = 0.2
else:
    n_prime = 1.0
C_s = (v.L_bc/v.L_bo)**(2-n_prime)+(v.L_bc/v.L_bi)**(2-n_prime)
#print(C_s)

#Pressure Drop Variable Calculations
A_frt = (np.pi/4)*(v.d_o**2)*F_w*v.N_t
A_frw = (np.pi/4)*(v.D_s**2)*((theta_b/(2*np.pi))-((np.sin(theta_b))/(2*np.pi)))
A_o_w = A_frw-A_frt
g_c = 1
G_w = v.m_s/((A_o_cr*A_o_w)**0.5)
b = 6.59/(1+0.14*Re_s**0.52)
f_id = 3.5*(1.33*v.d_o/v.p_t)**b*Re_s**(-0.476)
Delta_p_b_id = (4*v.N_rcc*f_id*G_s**2)*(1)**(0.25)/(2*g_c*v.p_s)
Delta_p_w = N_b*(2+0.6*N_rcw)*((G_w**2)/(2*g_c*v.p_s))*C_l
Delta_p_cr = Delta_p_b_id*(N_b-1)*C_b*C_l
Delta_p_io = 2*Delta_p_b_id*(1+(N_rcw/v.N_rcc))*C_b*C_s

#Total Pressure drop
Delta_p_s = Delta_p_cr+Delta_p_w+Delta_p_io
#print(Delta_p_s)

##Energy Balance Equation
T_to = v.T_ti+v.E*(v.T_si-v.T_ti)
#q = m_s*c_ps*(T_to-T_ti)

##Overal heat Transfer Coefficient
#h_t = (Nu*k)/d_i
#h_i = h_t
#U_o = (1/(h_o)) + (1/(h_of)) + ((d_o*ln(d_o/d_i))/(2*k_w)) + (d_o/(h_if*d_i)) + (d_o/(h_i*d_i))

##Approximate Design Method (Sizing)
#A_s = (q)/(U_o*F*Delta_T_lm)

##Ploting Variables
'''
#Presure drop vs Interior Diameter
plt.scatter(v.D_s, Delta_p_s)
#plt.axis([0.31, 0.35,])
plt.title("Presure Drop vs Interior Diameter")
plt.xlabel("Inside Diameter (m)")
plt.ylabel("Total Pressure Drop")
plt.show()

#Heat transfer coefficent vs Interior Diameter
plt.scatter(v.D_s, h_id)
plt.title("Heat transfer Coefficent vs Interior Diameter")
plt.xlabel("Inside Diameter (m)")
plt.ylabel("Heat Transfer Coefficent")
plt.show()
'''
