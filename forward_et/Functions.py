import os
import pickle

import numpy as np


def calFIPAR(NDVI):
    m2 = 1.0
    b2 = -0.05
    fIPAR = m2 * NDVI + b2
    return fIPAR


def GetConstraintsInputs(datadict):
    # THis function gets only those inputs from the dictionary that are needed for
    # calculating constraints
    NDVI = datadict['NDVI_16']
    NDVI = np.float64(NDVI)
    NDVI[NDVI == 0] = np.nan
    NDVI = NDVI/10000
    fPAR = datadict['fPAR']
    fPAR = np.float64(fPAR)
    fPAR = fPAR / 100
    return NDVI, fPAR


def OpenDictionay(Filename):
    # THis function opens the dictionary
    datadict = pickle.load(open(Filename, 'rb'))
    return datadict


def GetShortWaveRadiationInputs(datadict):
    """Units:
        SunHours = hours*10
        Latitude = Degress
        Doy = Float64
    """
    # THis function gets only those inputs from the dictionary that are needed for
    # calculating the shortwave incoming radiation.
    SunHours = datadict['SunHours'][:]
    folder = (os.path.dirname(os.path.realpath(__file__)))
    Latitude = np.load(folder+'\\Data\\Latitude_Andalucia.npy')

    Doy = datadict['JulianDay']
    Albedo = datadict['Albedo']
    Albedo = Albedo*1.0/100.
    return SunHours, Latitude, Doy, Albedo


def GetLongWaveRadiationInputs(datadict):
    # THis function gets only those inputs from the dictionary that are needed for
    # calculating the shortwave incoming radiation.
    LST = datadict['LST_Composite'][:]  # land surface temperature
    ObsTime = datadict['LST_Composite_Time_day'][:]
    Emiss = datadict['emis_Composite'][:]
    Emiss = Emiss*0.02
    AirMax = (np.float64(datadict['DASEMON_Tmax'][:])/10)+273.15
    AirMin = (np.float64(datadict['DASEMON_Tmin'][:])/10)+273.15
    """Units:
        LST = Kelvin
        ObsTime = hours

    """
    return LST, ObsTime, Emiss, AirMax, AirMin


def calcJparameter(N, t):
    J = (2*N)/(np.pi*np.sin((np.pi*t)/N)*24)
    return J


def GetLAI(datadict):
    """Units:
        LAI = (m2/m2)
    """
    # THis function gets the leaf area index from the dataset. LAI units are scaled by 10
    # and need to be divided by that factr get then in the right units.
    LAI = datadict['LAI'][:]  # land surface temperature
    LAI = np.float64(LAI)/10.0
    return LAI


def GetMeanTemperature(datadict):
    Tmax = datadict['DASEMON_Tmax'][:]
    Tmin = datadict['DASEMON_Tmin'][:]
    Tmean = (np.float64(Tmax)+np.float64(Tmin))/20
    Tmean[Tmean == 0] = np.nan
    return Tmean
