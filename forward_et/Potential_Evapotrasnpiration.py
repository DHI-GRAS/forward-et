# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:20:22 2018

@author: gmgo
"""

def calPsycometricConstant(Tmean):#(psi)
    import gdal
    import os
    import numpy as np
    folder=(os.path.dirname(os.path.realpath(__file__)))
    ASL_Altitude_m=np.load(folder+'\\Data\\SRTM_Andalucia_1100m.npy')
    P=101.3*(((293-0.0065*ASL_Altitude_m)/293)**5.26)# (kPa)
    #Latent Heat of Vaporization (lambda)
    Lambda=2454-2.4*(Tmean-20)
    #Specific heat of moist air (cp
    cp=1013 #(kJ kg-1 ºC-1)
    #Ratio molecular weight of water vapour/dray air
    epsi=0.622;    
    #Psychrometric constant (psi)
    psi=(cp*P)/(epsi*Lambda)# (Pa)
    return psi,Lambda
   
def calSatVapPres(Tmean):
    es=10*(0.061121*2.718281828**(17.502*Tmean/(240.97+Tmean))*(1.0007+(3.46*10**(-8)*100)));
    return es
    
def calSlopVapPresCurv(Lambda,Tmean,es):
    #Slope vapour Pressure curve (s)
    s=(Lambda*18*1000*es)/(8.3144*(Tmean+273)**2)
    return s
def calPotET_CANOPY(s,psi,alphaPT,RnCanopy):
    import numpy as np
    #Results canopy potential evapotranspiration Etpc % Following Pristley Taylor (1972)(eq 14)    
    Etpc=alphaPT*(s/(s+psi))*RnCanopy#(Wm-2)
    Etpc[np.isinf(Etpc)]=np.nan
    psi[np.isinf(psi)]=np.nan
    s[np.isinf(s)]=np.nan
    return Etpc
    
def calPotET_SOIL(s,psi,alphaPT,Rn_SOIL): 
    import numpy as np
    G=0;
    Etps=alphaPT*(s/(s+psi))*(Rn_SOIL-G)#Wm-2)
    Etps[np.isnan(Rn_SOIL)]=np.nan
    return Etps