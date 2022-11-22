from WATRO import WATRO
# from WATRO1 import WATRO1

"""Enter major ions concentrations in mol/l"""
Cl=971.0/35453.0; SO4=409.0/96000.0; Na=577.5/22989.769; Mg=83.0/24305.0
K=16.0/39098.3; Ca=167.0/40078.0; Sr=0.0; Br=0.0

"""Enter acid-base parameters""" 
feed_pH = 8.0 # Enter pH 
Bt_feed = 1.0 # Enter total boron (mg/l)
Alk_feed = 0.004789 #0.002992 # Enter Feed Alkalinity (eq/l)

"""Enter process operational conditions"""
P_feed = 9.0 #Enter Pressure (bars) 
t = 25.0 #Enter Temperature (celcius) 
u0 = 0.1 #Enter feed cross-flow velocity (m/s)
recovery = 82.0 #Enter Recovey Ratio (%)
pressure_drop = 2.0 #Enter total pressure drop (bars)

"""Enter Membrane Constants at 25C"""
Pw0 = 12.65e-7 #Enter water permeabiliy (if unavailable enter 0 - value will be derived from manufacturer data)
Ps0 = 9.404e-8 #Enter NaCl permeabiliy (if unavailable enter 0)
ks = 2.9404e-4 #Enter average mass transfer coefficient for charged solutes (if unavailable enter 0 - value will be derived from Sherwood correlations)
Pb0 = 2.78e-5 #Enter B(OH)3 permeabiliy (if unavailable enter 0)
kb = 2.62e-3 #Enter average mass transfer coefficient for unchraged solutes (if unavailable enter 0)
 
"""Enter Standard test conditions for estimating missing membrane constants"""
P_std = 41.0 #Enter standard pressure (bars)
NaCl_std = 32.0 #Enter standard NaCl concentration (g/l)
B_std = 5.0 #Enter standard B concentration (mg/l)
recovery_std = 50.0 #Enter recovery at standard conditions(%)
A = 7.9 #Enter Membrane surface area (m^2)
Qw = 4.7 #Enter Permeate flow at standard test conditions (m^3/d)
Rej_NaCl = 99.5 #Enter NaCl rejection at standard test conditions (%)
Rej_B = 83.0 # Enter B rejection at standard test conditions (%)
d_mil = 32.0 #enter feed spacer height (mil)

"""Run the program by pressing F5"""

"""The call for the function"""
(Jw,Cb,Cp,pH_b,pH_p,pH_m,Alkb,Alkm,
 Alkp,Btb,Btp,Btp_Accum_mgl,Ctb,Ctp)=WATRO(Cl,SO4,Na,Mg,K,Ca,Sr,Br,feed_pH,Bt_feed,Alk_feed,
                                           P_feed,t,u0,recovery,Pw0,Ps0,ks,Pb0,kb,P_std,
                                           NaCl_std,B_std,A,Qw,Rej_NaCl,Rej_B,d_mil,pressure_drop)

"""Output printing"""
print ('\n')
print ('\n')
print ('Accumulated permeate B (mg/l)')
print ('\n'.join(map(str, Btp_Accum_mgl)))
print ('\n')
print ('Brine pH')
print ('\n'.join(map(str, pH_b)))
print ('\n')
print ('Permeate pH')
print ('\n'.join(map(str, pH_p)))
print ('\n')
print ('\n')
print ('Brine_NaCl')
print ('\n'.join(map(str, Cb)))
print ('\n')
print ('Brine Alkalinity')
print ('\n'.join(map(str, Alkb)))
print ('\n')

 
"""Membrane Elements
HRLE440 = (41.0,31.0,92.0,99.8,8.0) #Units: A[m^2] Qw[m^3/d] Rej[%] Recovery[%]
XLE440 = (41.0,37.44,91.5,99.8,10.0)
ULE440 = (41.0,45.4,89.0,99.7,8.0)
BW30HR440 = (41.0, 48.0, 83.0, 99.7, 15.0)"""
