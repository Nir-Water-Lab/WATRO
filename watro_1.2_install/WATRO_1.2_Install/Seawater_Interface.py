
"""This is the original version of the WATRO program as published in
Environmental Science: Water Research and Technology (2015).
Updated version can be downloaded from https://dl.dropboxusercontent.com/u/64822322/WATRO.zip"""


import os
import numpy as np
from WATRO import WATRO

"""Enter major ions concentrations in mol/l"""
Cl=0.567125511; SO4=0.033763917; Na=0.497273368; Mg=0.052867719
K=0.011085392; Ca=0.010279573; Sr=0.0000776; Br=0.0

"""Enter acid-base parameters""" 
feed_pH = 8.0 # Enter pH 
Bt_feed = 5.0 # Enter total boron (mg/l)
Alk_feed = 0.002524
"""Enter process operational conditions"""
P_feed = 56.2 #Enter Pressure (bars) 
t = 25.0 #Enter Temperature (celcius) 
u0 = 0.17 #Enter feed cross-flow velocity (m/s)
recovery = 52.0 #Enter Recovey Ratio (%)
pressure_drop = 0.3 #Enter total pressure drop (bars)

"""Enter Membrane Constants at 25C. If unavailable enter 0 and it will be estimated by the software according to membrane manufacturer performance report"""
Pw0 = 5.793e-7 #1.084e-6 #Enter water permeabiliy (if unavailable enter 0 - value will be derived from manufacturer data)
Ps0 = 1.946e-8 #7.77e-8 #Enter NaCl permeabiliy (if unavailable enter 0)
ks = 2.32e-5 #7.73e-6 #Enter average mass transfer coefficient for charged solutes (if unavailable enter 0 - value will be derived from Sherwood correlations)
Pb0 = 1.875e-6 #2.78e-5 #Enter B(OH)3 permeabiliy (if unavailable enter 0)
kb = 1.0e-3 #1.62e-5 #Enter average mass transfer coefficient for unchraged solutes (if unavailable enter 0)
 
"""Enter manufacturer results from tandard test conditions for estimating missing membrane constants"""
P_std = 41.0 #Enter standard pressure (bars)
NaCl_std = 32.0 #Enter standard NaCl concentration (g/l)
B_std = 5.0 #Enter standard B concentration (mg/l)
recovery_std = 15.0 #Enter recovery at standard conditions(%)
A = 7.9 #Enter Membrane surface area (m^2)
Qw = 4.7 #Enter Permeate flow at standard test conditions (m^3/d)
Rej_NaCl = 99.5 #Enter NaCl rejection at standard test conditions (%)
Rej_B = 83.0 # Enter B rejection at standard test conditions (%)
d_mil = 28.0 #enter feed spacer height (mil)

"""Run the program by pressing F5"""

"""The call for the function"""
(Jw,Cb,Cp,pH_b,pH_p,pH_m,Alkb,Alkm,
 Alkp,Btb,Btp,Btp_Accum_mgl,Ctb,Ctp)=WATRO_no_salt(Cl,SO4,Na,Mg,K,Ca,Sr,Br,feed_pH,Bt_feed,Alk_feed,
                                           P_feed,t,u0,recovery,Pw0,Ps0,ks,Pb0,kb,P_std,
                                           NaCl_std,B_std,A,Qw,Rej_NaCl,Rej_B,d_mil,pressure_drop)

"""Output printing"""
print '\n Done\n Printing Results: \n'

print 'Accumulated permeate B (mg/l)'
print '\n'.join(map(str, Btp_Accum_mgl))
print '\n'
print 'Brine pH'
print '\n'.join(map(str, pH_b))
print '\n'
print 'Permeate pH'
print '\n'.join(map(str, pH_p))
print '\n'
print 'Permeate Alkalinty (eq/l)'
print '\n'.join(map(str, Alkp))
print '\n'
print 'Permeate flux (m/s)'
print '\n'.join(map(str, Jw))
print '\n'
print 'Permeate NaCl (mol/l)'
print '\n'.join(map(str, Cp))
print '\n'
print 'Brine_NaCl'
print '\n'.join(map(str, Cb))
print '\n'
print 'Brine Alkalinity'
print '\n'.join(map(str, Alkb))
print '\n'
print 'Momentary Permeate B (Mol/l)'
print '\n'.join(map(str, Btp))

 
"""Membrane Elements
HRLE440 = (41.0,31.0,92.0,99.8,8.0) #Units: A[m^2] Qw[m^3/d] Rej[%] Recovery[%]
XLE440 = (41.0,37.44,91.5,99.8,10.0)
ULE440 = (41.0,45.4,89.0,99.7,8.0)
BW30HR440 = (41.0, 48.0, 83.0, 99.7, 15.0)"""
