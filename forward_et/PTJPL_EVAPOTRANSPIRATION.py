# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 14:00:34 2018

@author: gmgo
"""

"""
This script calculates Evapotranspiration. 
This script is a continuation of the Net Radiation script and potential evapotranspiration. Therefore all comments to those code lines and indications
are removed to make it more easy to follow. It continues the same structure as the Net radiation script.
"""
import Functions as fn
import Radiation_Functions as rf
from Save_Image_As_Geotiff import getSceneMetadata, Save_array_tiff, SaveActEvapotrans_Image,SaveActSoilEva_Image,SaveActCanopyTransT_Image
import Potential_Evapotrasnpiration as PotET
import numpy as np
import Constraints

"""
The script starts getting a list of the files that are going to be proccesed. To do that we use the package glob
that has to be imported.
"""
import glob
InputsFolder='D:\\Tutorial 7 -GMGO\\1.Radiation\\Dictionaries_Andalucia\\'  #Specify the folder you have your dictionary data in.
FileList=glob.glob(InputsFolder+'*.pkl')

#Creates the output stack of images with the dimensions required for Andalusia
ActTranspiration=np.zeros([276,478,366])
ActEvaporation=np.zeros([276,478,366])
ActEvapotranspiration=np.zeros([276,478,366])
Fzang_all=np.zeros([366,276,478])
Fdrying_all=np.zeros([276,478,366])

i=0
for inputfile in FileList:
    Dictionary=fn.OpenDictionay(inputfile)
    SunHours,Latitude,Doy,Albedo=fn.GetShortWaveRadiationInputs(Dictionary)
    N,Ra,Rs,ws=rf.CalShortWaveIncomingRadiation(SunHours,Latitude,Doy)
    Sunset=(((ws)*(180/np.pi))/15)+12
    Sunrise=Sunset-N
    RSnet=rf.CalShortWaveNetRadiation(Rs, Albedo)
    LST,ObsTime,Emissivity,AirMax,AirMin=fn.GetLongWaveRadiationInputs(Dictionary)
    
    Rlw_in,Tair_MODIS_passtime=rf.CalLongWave_Incoming_Parton_Logan(LST,ObsTime,N,AirMax,AirMin)
    Rlw_out=rf.CalLongWave_Outgoing(LST,Emissivity)
    RLnet=Rlw_in-Rlw_out
    t=ObsTime-Sunrise
    J=fn.calcJparameter(N,t)
    RSnetInstant=RSnet/J
    RnDaily=(RSnetInstant+RLnet)*J
    LAI=fn.GetLAI(Dictionary)
    Tmean=fn.GetMeanTemperature(Dictionary)
    kRn=0.6
    Rn_SOIL,Rn_CANOPY=rf.SplitRn_SOIL_CANOPY(RnDaily,LAI,kRn) #W/m2
    psi,Lambda=PotET.calPsycometricConstant(Tmean)
    es=PotET.calSatVapPres(Tmean)
    s=PotET.calSlopVapPresCurv(Lambda,Tmean,es)
    alphaPT=1.26
    PotET_CANOPY=PotET.calPotET_CANOPY(s,psi,alphaPT,Rn_CANOPY)
    PotET_SOIL=PotET.calPotET_SOIL(s,psi,alphaPT,Rn_SOIL)   
    PotEVTR=PotET_CANOPY+PotET_SOIL
    """
    ****************************************
    CONSTRAINTS TO EVAPOTRANSPIRATION
    ****************************************
    In the previous part we have calculated the radiation components and the potential evapotranspiration. 
    One of the methods used to calculate ET using remote sensing data is the one proposed by Fisher and modified by Garcia et al.
    This method applies constraints to the potential ET and calculate the actual ET.
    In the next lines of code the functions needed to calculate these contraints are presented. Visit the papers from Fisher and Garcia
    to understand more on the nature of the constraints.
    
    """
    NDVI,FPAR=fn.GetConstraintsInputs(Dictionary)
    """
    BIOPHYSICAL CONSTRAINTS
    There are three constraints that need to be calculated:
    1. GREEN CANOPY FRACTION 
    2. PLANT TEMPERATURE CONSTRAINT
    3. PLANT MOISTURE CONSTRAINT
    """
    """GREEN CANOPY FRACTION  (fg)"""
    fg=Constraints.cal_fg(NDVI,FPAR)
    
    """PLANT TEMPERATURE CONSTRAINT (ft)"""
    ft=Constraints.cal_ft(Tmean)
    
    """PLANT MOISTURE CONSTRAINT (fm)"""
    fm=Constraints.calVegConstrains_fm(FPAR)
    #Eqs_all=Constraints.calEQ_all(Rn_SOIL,psi,s,count,Eqs_all)
    
    """SOIL MOISTURE CONSTRAIN"""
    Fzang_all,Fdrying_all=Constraints.calSOILMoistureConstrain_f_Zhang(i,Fzang_all,Fdrying_all)
        
    """
    ****************************************
    APPLY THE CONSTRAINTS TO THE POTENTIAL ET
    ****************************************
    """
    """ CONSTRAINTS TO TRASNPIRATION"""
    ET_DAILY_CANOPY=PotET_CANOPY*fm*ft*fg

    """ CONSTRAINTS TO EVAPORATION"""
    ET_DAILYSOIL=PotET_SOIL*Fdrying_all[:,:,i]

    """ TOTAL EVAPOTRASNPIRATION"""
    TOTAL_DAILY_ET_Constrained=ET_DAILYSOIL+ET_DAILY_CANOPY
    TOTAL_DAILY_ET_Constrained[np.isnan(TOTAL_DAILY_ET_Constrained)]=999.0
    TOTAL_DAILY_ET_Constrained[TOTAL_DAILY_ET_Constrained==999]=np.nan
    ActTranspiration[:,:,i]=ET_DAILY_CANOPY
    ActEvaporation[:,:,i]=ET_DAILYSOIL
    ActEvapotranspiration[:,:,i]=TOTAL_DAILY_ET_Constrained
    print(inputfile+' processed!')
    i=i+1
SaveActCanopyTransT_Image(ActTranspiration)
SaveActSoilEva_Image(ActEvaporation)
SaveActEvapotrans_Image(ActEvapotranspiration) #It creates and output in the output folder
    