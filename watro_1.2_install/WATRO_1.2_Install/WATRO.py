def WATRO(Cl,SO4,Na,Mg,K,Ca,Sr,Br,feed_pH,Bt_feed,Alk_feed,P_feed,t,u0,recovery,Pw0,Ps0,ks,Pb0,kb,P_std,NaCl_std,B_std,A,Qw,Rej_NaCl,Rej_B,d_mil,pressure_drop):
    
    # Import standard library modules first.
    import os
    import sys
    # Then get third party modules.
    from win32com.client import Dispatch 
    #import matplotlib.pyplot as plt 
    import numpy as np 
    from math import exp,sqrt
    import scipy.optimize as optimize

    def selected_array(db_path, input_string):
        """Load database via COM and run input string."""
        dbase = Dispatch('IPhreeqcCOM.Object')
        dbase.LoadDatabase(db_path)
        dbase.RunString(input_string)
        return dbase.GetSelectedOutputArray()

    def phreecalc(input_string):
        """Get results from PHREEQC"""
        pitzer_result = selected_array('pitzer.dat', input_string)
        return pitzer_result

    def visco(T,S):
        """Calculate sewater viscosity based on Sharqawy et al. 2009"""
        S=S/1000
        a1 = 1.5700386464E-01;a2 = 6.4992620050E+01;a3 = -9.1296496657E+01;a4 = 4.2844324477E-05;
        mu_w = a4 + 1/(a1*(T+a2)**2+a3)
        a5 = 1.5409136040E+00;a6 = 1.9981117208E-02;a7 = -9.5203865864E-05
        a8 = 7.9739318223E+00;a9 = -7.5614568881E-02;a10 = 4.7237011074E-04
        A = a5 + a6*T + a7*T**2; B = a8 + a9*T + a10*T**2
        mu = mu_w*(1 + A*S + B*S**2)
        return mu

    def func(jw):
        """Calculate Jw using solution-diffusion and film layer models"""
        z=exp(jw/k[i])
        cm=((Ps+jw)*Cb[i]*z)/(jw+Ps*z)
        cp=Ps*Cb[i]*z/(jw+Ps*z)
        return Pw*(Pbar[i]-(PHI*cm-0.98*cp)*T*0.083145) - jw
      
    """numerical parameters"""
    step_num= int(recovery + 1)
    """Initialization of variables"""
    T=t+273.15; Ppa=P_feed*1e5; kphi=0; PHI=1.0; PHI_old=0
    r_f=recovery/100.0; r=np.linspace(0,r_f,step_num); dr=r[1]-r[0]; d = d_mil*2.54e-5; Pbar = np.zeros(len(r))    
    S=np.zeros(len(r));Cb=np.zeros(len(r));Ctb=np.zeros(len(r));Ctp=np.zeros(len(r));Alkb=np.zeros(len(r));Alkp=np.zeros(len(r))
    S0=(Cl*35.453+Na*22.98977+SO4*96.0626+Mg*24.305+Ca*40.078+K*39.098+Br*79.904+Sr*87.62+Bt_feed/1000);
    Cb[0]=(Cl+Na+SO4+Mg+Ca+K+Br+Sr+Bt_feed/10811); Cp=np.zeros(len(r)); C=np.zeros(len(r))
    Cm=np.zeros(len(r));k=np.zeros(len(r));Jw=np.zeros(len(r));CFb=np.zeros(len(r));CF=np.zeros(len(r));pH_b=np.zeros(len(r))
    pH_m=np.zeros(len(r)); pH_p=np.zeros(len(r));Btp=np.zeros(len(r));Btb=np.zeros(len(r));Btp_Accum=np.zeros(len(r));Btp_Accum_mgl=np.zeros(len(r))
    BOH3_b=np.zeros(len(r));BOH4_b=np.zeros(len(r));CO2_b=np.zeros(len(r));HCO3_b=np.zeros(len(r));CO3_b=np.zeros(len(r))
    Theta=np.zeros(len(r)); w_H_eff=np.zeros(len(r)); w_OH_eff=np.zeros(len(r))

    """Get constants from standard test conditions"""
    Bp=B_std*(1-Rej_B/100); NaClp=NaCl_std*(1-Rej_NaCl/100)
    Jw_avg = Qw/(A*24*3600)
    PHI_avg = 0.98; phi_boric = 0.95
    if NaCl_std>30.0:
        PHI_avg=0.922
        phi_boric= 0.9  #obtained from PHREEQC for the test conditions
    Pi= PHI_avg*2*(NaCl_std/58.443)*0.083145*T
    Pw_std= Jw_avg/(P_std-Pi)
    Ps_std= (Jw_avg*NaClp)/(NaCl_std - NaClp)
    Pb_std= (Jw_avg*Bp)/(phi_boric*B_std - Bp)
    
    """Temperature corrections for permeability constants"""
    if Pw0==0:
        Pw0=Pw_std
    if Ps0==0:
        Ps0=Ps_std
    if Pb0==0:
        Pb0=Pb_std

    Pw= Pw0*exp(0.0093*(T - 298.15))  #Taniguchi et al. 2001
    Pb= Pb0*exp(0.067*(T - 298.15))  #Hyung and Kim 2006
    Ps= Ps0*exp(0.0483*(T - 298.15))  #Taniguchi et al. 2001    
   
    """Pw=Pw0*exp(2640*(1/T-1/298.15))  #alternative ROSA equation"""
    Feed = """
            SOLUTION 1 seawater
            units     mol/l
            temp     %f
            pH       %f
            Cl       %e
            S(6)     %e
            Br       %e   
            Na       %e
            Mg       %e
            K        %e
            Ca       %e
            Sr       %e
            Alkalinity       %e
            B        %f mg/l
            USE solution 1
            USER_PUNCH
            -headings RHO osmotic_coefficient ALK Ct Bt
            -start           
            10 PUNCH ALK
            20 PUNCH TOT("C")
            30 PUNCH TOT("B")
            40 PUNCH RHO
             -end
            SELECTED_OUTPUT
            -reset          false
            -user_punch     true
            END"""%(t,feed_pH,Cl,SO4,Br,Na,Mg,K,Ca,Sr,Alk_feed,Bt_feed)
    sol_feed = phreecalc(Feed)
    rho= sol_feed[1][3]
    Alkb[0]=sol_feed[1][0] #/((1+S0/1000)/(rho/1000))
    Ctb[0]=sol_feed[1][1] #/((1+S0/1000)/(rho/1000))
    Btb[0]=sol_feed[1][2] #/((1+S0/1000)/(rho/1000))
    
    print ('Calculating water flux and salt concentration:'    )
    for i in range (len(r)):        
        """Water Flux and salt passage model"""
        Pbar[i]= P_feed - pressure_drop*(r[i]/r[len(r)-1])           
        PHI = 1.0
        """mass transfer coefficient"""
        k[i]= ks
        if ks==0:
            RHO_PHREE = """
                SOLUTION 1 seawater
                units     mol/l
                temp     %f
                pH       %f
                Cl       %e
                S(6)     %e
                Br       %e   
                Na       %e
                Mg       %e
                K        %e
                Ca       %e
                Sr       %e
                USE solution 1
                REACTION_PRESSURE 1
                %f
                USER_PUNCH
                -headings RHO osmotic_coefficient ALK Ct Bt
                -start
                10 PUNCH RHO
                20 PUNCH OSMOTIC
                 -end
                SELECTED_OUTPUT
                -reset          false
                -user_punch     true
                END"""%(t,7,Cl/(1-r[i]),SO4/(1-r[i]),Br/(1-r[i]),Na/(1-r[i]),Mg/(1-r[i]),K/(1-r[i]),Ca/(1-r[i]),Sr/(1-r[i]),Pbar[i])
            sol_rho = phreecalc(RHO_PHREE)
            rho = 1000*sol_rho[2][0]
            PHI = sol_rho[2][1]   
            S[i]=S0/(1-r[i])    # bulk salinity in kg/m^3
            visc=visco(t,S[i])   # (1.234*10**-6)*exp(0.00212*S[i]+1965/T) seawater viscocity in pa*s  from Sharkwy et al. 2009
            D_NaCl=(6.725*10**-6)*exp(0.1546*S[i]*10**-3 - 2513/T)   # Diffusivity  of NaCl in seawater in  m^2/s  from taniguchi et al 2001
            Sh=0.065*((rho*u0*(1-r[i])*2*d/visc)**0.875)*(visc/(rho*D_NaCl))**0.25  # sherwood number  from taniguchi et al 2001
            k[i] = Sh*D_NaCl/d   # mass transfer coefficient in m/sec                 
        
        """find Jw(i), Cm(i) and PHI(i)"""
        CF[i] = 1/(1-r[i])
        PHI_old =10
        while (abs(PHI-PHI_old)>0.001):
            
            PHI_old=PHI
            osmo_phree = """
                SOLUTION 1 mediterranean seawater
                units      mol/l
                temp       %f
                pH         %f
                Cl         %e 
                S(6)       %e 
                Br         %e  
                Na         %e 
                Mg         %e 
                K          %e 
                Ca         %e 
                Sr         %e              
                USE solution 1            
                USER_PUNCH
                -headings osmotic_coefficient
                -start
                10 PUNCH OSMOTIC
                20 PUNCH RHO
                 -end
                SELECTED_OUTPUT
                -reset                false
                -user_punch           true
                 END"""%(t,7,Cl*CF[i],SO4*CF[i],Br*CF[i],Na*CF[i],Mg*CF[i],K*CF[i],Ca*CF[i],Sr*CF[i])

            sol_osm=phreecalc(osmo_phree)
            PHI = sol_osm[1][0]; rho = sol_osm[1][1]        
            Jw[i] = optimize.bisect(func, 1e-8 ,1e-4 , xtol=1e-17, rtol=5e-15, maxiter=500)
            Cp[i] = (Cb[i]*Ps*exp(Jw[i]/k[i]))/(Jw[i]+Ps*exp(Jw[i]/k[i]))
            Cm[i] = Cp[i] +(Cb[i]-Cp[i])*exp(Jw[i]/k[i])
            CF[i] = Cm[i]/Cb[0]                   
            kphi=kphi+1
        
        if r[i]<recovery/100:
            Cb[i+1] = (Cb[i]*(1-r[i]) - dr*Cp[i])/(1-r[i+1])        
        CFb[i] = Cb[i]/Cb[0]        
    print ('Done \n\n'  )

     
    print ('Calculating acid-base dynamics: \n\n Recovery ratio count: \n')
    for i in range(len(r)-1):
        print(r[i])
                
        """Reactive transport model"""     
        bulk_speciation = """
            SOLUTION 1 seawater
            units     mol/kgw
            temp        %f
            pH          %f
            Cl          %e 
            S(6)        %e 
            Br          %e   
            Na          %e 
            Mg          %e 
            K           %e 
            Ca          %e 
            Sr          %e 
            C           %e 
            B           %e 
            Alkalinity    %e 
            USE solution 1
            REACTION_PRESSURE 1
            %f
            SELECTED_OUTPUT
            -reset    false
            -high_precision     true
            -ph       true
            -molalities      B(OH)4-  B(OH)3  HCO3-  CO2  CO3-2  OH-  H+  MgOH+  HSO4-  MgCO3
             END"""%(t,7.0,Cl*CFb[i],SO4/(1-r[i]),Br*CFb[i],Na*CFb[i],Mg/(1-r[i]),K*CFb[i],Ca/(1-r[i]),Sr/(1-r[i]),Ctb[i],Btb[i],Alkb[i],Pbar[i])
            
        sol=phreecalc(bulk_speciation)
        pH_b[i]=sol[2][0];BOH4_b[i]=sol[2][1];BOH3_b[i]=sol[2][2];HCO3_b[i]=sol[2][3];CO2_b[i]=sol[2][4]
        OH_b=sol[2][6]; H_b=sol[2][7]; MgOH_b=sol[2][8]; HSO4_b = sol[2][9]; MgCO3_b = sol[2][10]
        BOH4_b[i] = Btb[i] - BOH3_b[i]; CO3_b[i] = Ctb[i] - HCO3_b[i] - CO2_b[i]       
       
        if kb==0:
            kb=2*k[i]

        if i==0:            
            BOH3_p= (Pb*BOH3_b[0]*exp(Jw[i]/kb))/(Jw[0]+Pb*exp(Jw[i]/kb))
            BOH4_p= (Ps*BOH4_b[0]*exp(Jw[i]/k[i]))/(Jw[0]+Ps*exp(Jw[i]/k[i]))
            HCO3_p= (Ps*HCO3_b[0]*exp(Jw[i]/k[i]))/(Jw[0]+Ps*exp(Jw[i]/k[i]))
            CO2_p=  CO2_b[i]
            Btp[0]=BOH3_p+BOH4_p
            Ctp[0]=HCO3_p+CO2_p
            Btp_Accum[0]=Btp[0]
            
        
        OH_p = OH_b
        H_p = H_b
        MgOH_m=MgOH_b*exp(Jw[i]/k[i])
        HSO4_m=HSO4_b*exp(Jw[i]/k[i])      
        CO3_m=CO3_b[i]*exp(Jw[i]/k[i])
        MgCO3_m=MgCO3_b*exp(Jw[i]/k[i])
        CO2_m=CO2_b[i]        
        
        BOH3_p_old = 100; pH_m_old=100; pH_m[i]=pH_b[i]; Alkm= Alkb[i]*exp(Jw[i]/k[i]); Alkm_old=0 
        BOH3_p= (Pb*BOH3_b[i]*exp(Jw[i]/kb))/(Jw[i]+Pb*exp(Jw[i]/kb))
        kk=0        
        while(abs((pH_m[i]-pH_m_old)/pH_m[i])>0.0001)and(kk<50):            
            Alkm_old = Alkm
            pH_m_old = pH_m[i]
            BOH3_p_old = BOH3_p
            BOH3_m=BOH3_p+(BOH3_b[i]-BOH3_p)*exp(Jw[i]/kb)
            BOH4_m=BOH4_p+(BOH4_b[i]-BOH4_p)*exp(Jw[i]/k[i])
            HCO3_m=HCO3_p+(HCO3_b[i]-HCO3_p)*exp(Jw[i]/k[i])     
            OH_m = OH_p+(OH_b-OH_p)*exp(Jw[i]/(3.34*k[i]))
            H_m = H_p+(H_b-H_p)*exp(Jw[i]/(5.62*k[i]))            
                                 
            Btm= BOH3_m + BOH4_m
            
            Ctm= HCO3_m + CO2_m + CO3_m # + MgCO3_m
            Alkm= HCO3_m + 2*CO3_m + BOH4_m + OH_m - H_m + MgOH_m - HSO4_m # + 2*MgCO3_m          
            #Alkm= Alkb[i]*exp(Jw[i]/k[i])           
            film_speciation = """
                SOLUTION 1 mediterranean seawater
                units         mol/kgw
                temp          %f
                pH            %f
                Cl            %e 
                S(6)          %e 
                Na            %e 
                Mg            %e 
                K             %e 
                Ca            %e 
                Br            %e   
                Sr            %e
                C             %e 
                B             %e
                Alkalinity    %e 
                USE solution 1
                REACTION_PRESSURE 1
                %f
                SELECTED_OUTPUT
                -reset    false
                -high_precision       true
                -ph       true
                -molalities      B(OH)4-  B(OH)3  HCO3-  CO2  OH-  H+  MgOH+  HSO4- CO3-2 MgCO3
                 END"""%(t,7,CF[i]*Cl,CF[i]*SO4,CF[i]*Na,CF[i]*Mg,CF[i]*K,CF[i]*Ca,CF[i]*Br,CF[i]*Sr,Ctm,Btm,Alkm,Pbar[i])

            sol=phreecalc(film_speciation)
            pH_m[i]=sol[2][0];BOH4_m=sol[2][1];BOH3_m=sol[2][2];HCO3_m=sol[2][3];CO2_m=sol[2][4]
            OH_m = sol[2][5]; H_m = sol[2][6];  MgOH_m=sol[2][7]; HSO4_m = sol[2][8]
            CO3_m = Ctm - HCO3_m - CO2_m #sol[2][9]; MgCO3_m = sol[2][10]
            BOH4_m = Btm - BOH3_m
                        
            BOH3_p= (Pb*BOH3_m)/(Jw[i]+Pb)
            BOH4_p= (Ps*BOH4_m)/(Jw[i]+Ps)
            HCO3_p= (Ps*HCO3_m)/(Jw[i]+Ps)
            kk=kk+1            
                
            CO2_p= CO2_m
            
            k_Cb = 0.357/(1+exp(-52.63022629*(Cm[i]-0.12)))        
            Theta[i] = (1-k_Cb-0.05713078)/(1+exp(-1.72843187*(pH_m[i]-7))) + 0.05713078              
            w_H = 0.043; w_OH = 0.000033          
            w_H_eff[i] = w_H + (OH_m/H_m)*w_OH
            w_OH_eff[i] = w_OH + (H_m/OH_m)*w_H
            Rs = 1-Cp[i]/Cm[i]
            
            a= OH_m*w_OH_eff[i]
            b= w_OH_eff[i]*(1-Rs)**Theta[i]
            c=Jw[i]*(1-(1-Rs)**(1+Theta[i]))/(Rs*(1+Theta[i]))
            
            OH_p = a/(b+c)

            a2= H_m*w_H_eff[i]
            b2= w_H_eff[i]*(1-Rs)**(-Theta[i])
            c2=Jw[i]*(1-(1-Rs)**(1-Theta[i]))/(Rs*(1-Theta[i]))
                                     
            H_p = a2/(b2+c2)

            Btp[i]=BOH3_p+BOH4_p
            Ctp[i]=HCO3_p+CO2_p
            Alkp[i]=HCO3_p+BOH4_p +OH_p - H_p       
        
       
        permeate_speciation = """
            SOLUTION 1 permeate
            units         mol/kgw
            temp          %f
            pH            %f
            Na            %e 
            Cl            %e 
            C             %e 
            B             %e
            Alkalinity    %e
            USE solution 1
            SELECTED_OUTPUT
            -reset    false
            -ph       true
            -molalities      B(OH)4-  B(OH)3  HCO3-  CO2  OH-  H+ CO3-2            
             END"""%(t,7,Cp[i]/2,Cp[i]/2,Ctp[i],Btp[i],Alkp[i])
        
        sol=phreecalc(permeate_speciation)
        pH_p[i]=sol[1][0]; BOH3_p=sol[1][2]; HCO3_p=sol[1][3]; CO2_p=sol[1][4]; OH_p=sol[1][5]; H_p=sol[1][6]
        BOH4_p = Btp[i] - BOH3_p
        CO3_p = sol[1][7]
             
        Btb[i+1] = (Btb[i]*(1-r[i]) - dr*Btp[i])/(1-r[i+1])
        Ctb[i+1] = (Ctb[i]*(1-r[i]) - dr*Ctp[i])/(1-r[i+1])
        Alkb[i+1] = (Alkb[i]*(1-r[i]) - dr*Alkp[i])/(1-r[i+1])

    for t in range(len(r)):
        Btp_Accum[t] = np.average(Btp[0:t+1])
    Btp_Accum_mgl=10811*Btp_Accum

    i=i-1
    bulk_speciation = """
            SOLUTION 1 
            units     mol/kgw
            temp        %f
            pH          %f
            Cl          %e 
            S(6)        %e 
            Br          %e   
            Na          %e 
            Mg          %e 
            K           %e 
            Ca          %e 
            Sr          %e 
            C           %e 
            B           %e 
            Alkalinity    %e 
            USE solution 1
            REACTION_PRESSURE 1
            %f
            SELECTED_OUTPUT
            -reset    false
            -high_precision     true
            -ph       true
            -molalities      B(OH)4-  B(OH)3  HCO3-  CO2  CO3-2  OH-  H+  MgOH+  HSO4- MgCO3
            -totals               Ca
            -saturation_indices   Aragonite
            -equilibrium_phases   Aragonite
            EQUILIBRIUM_PHASES 1
                Aragonite 0 0
            END"""%(t,7,Cl*CFb[i],SO4/(1-r[i]),Br*CFb[i],Na*CFb[i],Mg/(1-r[i]),K*CFb[i],Ca/(1-r[i]),Sr/(1-r[i]),Ctb[i],Btb[i],Alkb[i],Pbar[i])
    
    

    return Jw,Cb,Cp,pH_b,pH_p,pH_m,Alkb,Alkm,Alkp,Btb,Btp,Btp_Accum_mgl,Ctb,Ctp
        
    

    
