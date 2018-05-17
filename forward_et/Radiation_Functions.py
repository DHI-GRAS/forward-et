# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 15:08:16 2018

@author: gmgo
"""


def CalShortWaveIncomingRadiation(SunHours, Latitude, Doy):
    import numpy as np
    from matplotlib import pyplot as plt

    """------------------Parameters---------------
    factors for interception of solar radiation by aerosols and clouds
    IMPORTANT!!!! Please read the pages 41 to 51 from the FAO Irrigation and Drainage paper No. 56.
    You can get the pdf document in this website: http://academic.uprm.edu/abe/backup2/tomas/fao%2056.pdf
    In those pages there is a detail description of how to calculate the SWIn raditation and the equations that
    are implemented in the lines below.
    """
    SunHours = SunHours*1./10.  # It needs to be divided by 10 to get in the right units
    SunHours[SunHours == 0] = np.nan
    AS = 0.25
    BS = 0.50
    latRad = (Latitude*np.pi)/180
    delta = 0.409*np.sin(((Doy*2*np.pi)/365)-1.39)  # Declination in radians
    dr = 1+0.033*np.cos(Doy*((2*np.pi)/365))  # Eccentricity factor
    # ; %sunset hour angle radians
    ws = np.arccos(-np.tan(latRad)*np.tan(delta))
    N = (ws*24.)/np.pi  # day length in hours
    Ra = (((24*60)/np.pi)*0.0820)*dr*((ws*np.sin(latRad)*np.sin(delta)) +
                                      (np.cos(latRad)*np.cos(delta)*np.sin(ws)))  # Mjm2day-1
    Rs = Ra*(AS+(BS*SunHours/N))  # Shortwave Incoming radiation

    return N, Ra, Rs, ws


def CalShortWaveNetRadiation(Rs, Albedo):
    from matplotlib import pyplot as plt
    # RRsSwD is the previously calculated Incoming Shortwave radiation
    # Albedo is the albedo obtained from MODIS.
    RSnet = Rs*(1-Albedo)
    RSnet = RSnet*11.574  # Converts the units into W/m2

    return RSnet


def CalLongWave_Incoming_Parton_Logan(LST, ObsTime, N, AirMax, AirMin):
    import numpy as np
    import os
    # c= lag of the minimum temperature from the time of sunrise,(parameter following table 1, p.210)
    c1 = -0.17
    # d= d refers to a in Parton and Logan, (parameter following table 1, p.210)
    d1 = 1.86
    BB = 12-(N/2)+c1
    ObsTime = np.float64(ObsTime)
    ObsTime[ObsTime > 24] = np.nan
    ObsTime[ObsTime == 0] = np.nan

    BBD = ObsTime-BB  # BBD, number of hours after the minimum temperature occcurs until sunset

    Tair_MODIS_passtime = (AirMax-AirMin)*np.sin((np.pi*BBD)/(N+2*d1))+AirMin
    Tair_MODIS_passtime[np.isnan(ObsTime)] = np.nan
    # Calculating incoming longwave radiation
    Stef = 5.67e-008
    C = 0.261
    d = 7.77e-004
    EmissivityAir = 1-C*np.exp(-d*(Tair_MODIS_passtime-273.15)**2)
    Rlw_in = Stef*(EmissivityAir)*(Tair_MODIS_passtime**4)

    return Rlw_in, Tair_MODIS_passtime


def CalLongWave_Outgoing(LST, Emissivity):
    import numpy as np
    Stef = 5.67e-008
    LST = np.float64(LST)
    Rlw_out = Stef*(Emissivity)*(LST**4)

    return Rlw_out


def SplitRn_SOIL_CANOPY(Rn, LAI, kRn):
    import numpy as np
    Rn_SOIL = Rn*np.exp(-kRn-LAI)
    Rn_SOIL[np.isnan(Rn)] = np.nan
    Rn_Canopy = Rn-Rn_SOIL

    return Rn_SOIL, Rn_Canopy
