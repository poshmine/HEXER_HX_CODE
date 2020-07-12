'''
Date of last edit: July 12th, 2020
Author(s): Ryan McGuire*    Lane Carasik^
*Virginia CommonWealth University
*FAST Research Group
SHEER - Shell-and-tube Heat ExchangER code
Equations for Shell-Side Reactor Heat Exchangers: 
Python Script for calculation of Shell-and Tube heat transfer coefficient using correctional factors
for Shell-and-Tube bundles based on models found in literature.
'''

'''
Revision Points
1) A_o_cr
2) A_o_tb
3) N_b
'''

##Basic Imports and Function Representatives
import numpy as np

##Measured Information used in calculations

##Correctional Factor Calculations

#Baffle configuration correctional factor
D_ctl = D_otl-d_o
theta_ctl = 2*np.acos((D_s-2*l_c)/D_ctl)
F_c = 1-(theta_ctl/np.pi)+(np.sin(theta_ctl)/np.pi)
J_c = (0.55+0.72*F_c)
#print(J_c)

#Bundle Leakage effects correctional factor
theta_b = 2*np.acos(1-(2*l_c/D_s))
F_w = (theta_ctl/(2*np.pi))-(np.sin(theta_ctl)/(2*np.pi))
d_sb = D_s-D_baffle
d_tb = d_l-d_o
A_o_sb = np.pi*D_s*(d_sb/2)*(1-(theta_b/(2*np.pi)))
A_o_tb = (np.pi*d_o*d_tb*N_t*(1-F_w))/2
A_o_cr = 0.03275
'''
if Bundle_layout in (30,90):
    A_o_cr = (D_s-D_otl+(D_ctl/X_t)*(X_t-d_o))*L_bc
elif Bundle_layout in (45,60):
    print("no1")
elif Bundle_layout in (30,90):
    print("no2")
else Bundle_layout in (30,90):
    print("no3")
'''
r_s = (A_o_sb)/(A_o_sb+A_o_tb)
r_lm = (A_o_sb+A_o_tb)/(A_o_cr)
J_l = (0.44*(1-r_s)+(1-0.44*(1-r_s))*np.e**(-2.2*r_lm))
#print(J_l)

#Bundle and pass partition bypass correctional factor
G_s = (m_s/A_o_cr)
A_o_bp = L_bc*(D_s-D_otl+0.5*N_p*w_p)
Re_s = round((G_s*d_o)/u_s)
N_ss_plus = N_ss/N_rcc
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
if Laminar_flow == 1:
    n = 1/3
else:
    n = 3/5
L_i_plus = L_bi/L_bc
L_o_plus = L_bo/L_bc
N_b = ((L-L_bi-L_bo)/L_bc)+1
J_s = ((N_b-1+L_i_plus**(1-n)+L_o_plus**(1-n))/(N_b-1+L_i_plus+L_o_plus))
#print(J_s)

#Adverse Temperature Gradient Buildup in Laminar Flow correctional factor
N_rcw = (0.8/X_t)*(l_c-(1/2)*(D_s-D_ctl))
N_rc = N_rcc+N_rcw
if Re_s >= 100:
    J_r = 1
elif Re_s <= 20:
    J_r = (10/N_rc)**0.18
else:
    J_r = 1+(((10/N_rc)**0.18)-1)/(20-100)*(Re_s-100)
#print(J_r)

##Heat Transfer Calculations - Core

#Perfect Heat Transfer Coefficent 
h_id = (((Nu_s*k)/d_o)*(1)**-0.14)
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
    C_b = e**(-D*r_b*(1-(2*N_ss_plus)**(1/3)))
#print(C_b)

#Bypass Flow
p = -0.15*(1+r_s)+0.8
C_l = e**(-1.33*(1+r_s)*r_lm**p)
#print(C_l)

#Differing Baffle Spacing from central section
if Laminar_flow == 0:    
    n_prime = 0.2
else:
    n_prime = 1.0
C_s = (L_bc/L_bo)**(2-n_prime)+(L_bc/L_bi)**(2-n_prime)
#print(C_s)

#Pressure Drop Variable Calculations
A_frt = (np.pi/4)*(d_o**2)*F_w*N_t
A_frw = (np.pi/4)*(D_s**2)*((theta_b/(2*np.pi))-((np.sin(theta_b))/(2*np.pi)))
A_o_w = A_frw-A_frt
g_c = 1
G_w = m_s/((A_o_cr*A_o_w)**0.5)
b = 6.59/(1+0.14*Re_s**0.52)
f_id = 3.5*(1.33*d_o/p_t)**b*Re_s**(-0.476)
Delta_p_b_id = (4*N_rcc*f_id*G_s**2)*(1)**(0.25)/(2*g_c*p_s)
Delta_p_w = N_b*(2+0.6*N_rcw)*((G_w**2)/(2*g_c*p_s))*C_l
Delta_p_cr = Delta_p_b_id*(N_b-1)*C_b*C_l
Delta_p_io = 2*Delta_p_b_id*(1+(N_rcw/N_rcc))*C_b*C_s

#Total Pressure drop
Delta_p_s = Delta_p_cr+Delta_p_w+Delta_p_io
#print(Delta_p_s)

##Energy Balance Equation
T_to = T_ti+E*(T_si-T_ti)
#q = m_s*c_ps*(T_to-T_ti)

##Overal heat Transfer Coefficient
h_t = (Nu*k)/d_i
h_i = h_t
#U_o = (1/(h_o)) + (1/(h_of)) + ((d_o*ln(d_o/d_i))/(2*k_w)) + (d_o/(h_if*d_i)) + (d_o/(h_i*d_i))

##Approximate Design Method (Sizing)
#A_s = (q)/(U_o*F*Delta_T_lm)
