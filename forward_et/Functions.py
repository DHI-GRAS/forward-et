# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 09:25:32 2018

@author: gmgo
"""


def calFIPAR(NDVI):
    m2 = 1.0
    b2 = -0.05
    fIPAR = m2 * NDVI + b2
    return fIPAR


def GetConstraintsInputs(Dictionary):
    # THis function gets only those inputs from the dictionary that are needed for
    # calculating constraints
    import numpy as np
    NDVI = Dictionary['NDVI_16']
    NDVI = np.float64(NDVI)
    NDVI[NDVI == 0] = np.nan
    NDVI = NDVI / 10000
    fPAR = Dictionary['fPAR']
    fPAR = np.float64(fPAR)
    fPAR = fPAR / 100

    return NDVI, fPAR


def OpenDictionay(Filename):
    # THis function opens the dictionary
    import pickle
    Dictionary = pickle.load(open(Filename, 'rb'))
    return Dictionary


def GetShortWaveRadiationInputs(Dictionary):
    # THis function gets only those inputs from the dictionary that are needed for
    # calculating the shortwave incoming radiation.
    import os
    import numpy as np
    SunHours = Dictionary['SunHours'][:]
    folder = (os.path.dirname(os.path.realpath(__file__)))
    Latitude = np.load(folder + '\\Data\\Latitude_Andalucia.npy')

    Doy = Dictionary['JulianDay']
    Albedo = Dictionary['Albedo']
    Albedo = Albedo * 1.0 / 100.
    """Units:
        SunHours=hours*10
        Latitude=Degress
        Doy=Float64
    """
    return SunHours, Latitude, Doy, Albedo


def GetLongWaveRadiationInputs(Dictionary):
    import numpy as np
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    # THis function gets only those inputs from the dictionary that are needed for
    # calculating the shortwave incoming radiation.
    LST = Dictionary['LST_Composite'][:]  # land surface temperature
    ObsTime = Dictionary['LST_Composite_Time_day'][:]
    Emiss = Dictionary['emis_Composite'][:]
    Emiss = Emiss * 0.02
    AirMax = (np.float64(Dictionary['DASEMON_Tmax'][:]) / 10) + 273.15
    AirMin = (np.float64(Dictionary['DASEMON_Tmin'][:]) / 10) + 273.15
    """Units:
        LST=Kelvin
        ObsTime=hours

    """

    return LST, ObsTime, Emiss, AirMax, AirMin


def calcJparameter(N, t):
    import numpy as np

    J = (2 * N) / (np.pi * np.sin((np.pi * t) / N) * 24)

    return J


def GetLAI(Dictionary):
    import numpy as np
    # THis function gets the leaf area index from the dataset. LAI units are scaled by 10
    # and need to be divided by that factr get then in the right units.
    LAI = Dictionary['LAI'][:]  # land surface temperature
    LAI = np.float64(LAI) / 10.0
    """Units:
        LAI=(m2/m2)
    """

    return LAI


def GetMeanTemperature(Dictionary):
    import numpy as np
    Tmax = Dictionary['DASEMON_Tmax'][:]
    Tmin = Dictionary['DASEMON_Tmin'][:]
    Tmean = (np.float64(Tmax) + np.float64(Tmin)) / 20
    Tmean[Tmean == 0] = np.nan
    return Tmean
