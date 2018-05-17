# -*- coding: utf-8 -*-
"""
Created on Wed Nov 08 09:19:12 2017

@author: gmgo
"""
#


def calcJparameter(N, t):
    import numpy as np
    J = (2*N)/(np.pi*np.sin((np.pi*t)/N)*24)

    return J


def calDailySolarIrradiance(Lat, lon, doy, In):
    import numpy as np

    # In-> is the insolation in hours obtained from Vicente Serrano Dataset
    """------------------Parameters---------------
    factors for interception of solar radiation by aerosols and clouds"""
    AS = 0.15
    BS = 0.30
#
#    AS=0.15
#    BS=0.10

    latRad = (Lat*np.pi)/180
#    lonRad=(lon*np.pi)/180
#    InRad=(In*np.pi)/180
    delta = 0.409*np.sin(((doy*2*np.pi)/365)-1.39)  # Declination in radians
    # delta=((-23.5*np.pi)/180)*np.cos(360/365*(doy+10))#Declination in radians
    dr = 1+0.033*np.cos(doy*((2*np.pi)/365))  # Eccentricity factor
    # ; %sunset hour angle radians
    ws = np.arccos(-np.tan(latRad)*np.tan(delta))
    N = (ws*24.)/np.pi  # day length in hours
    Ra = (((24*60)/np.pi)*0.0820)*dr*((ws*np.sin(latRad)*np.sin(delta)) +
                                      (np.cos(latRad)*np.cos(delta)*np.sin(ws)))  # Mjm2day-1
    Rs = Ra*(AS+(BS*In/N))
    # SunRiseHour

    return N, Ra, Rs, ws, In


def calculateRLDownInstIDsoandJackson(LST, Tmax, Tmin, latitude, longitude, year, step):
    import numpy as np
    LST = LST
    sigma = 5.67e-8  # Stefan Boltzman
    c1 = 0.261  # constant
    c2 = -0.17  # Lag time in minimum air temperature
    d = 7.77e-4  # constant
    # Assumed to be this time until we get the information regarding the observation time
    ObservationTime = 12.00
    # Converts the year and step from the netcdf file to gregorian calendar so it can be use in the function from noaa.
    latRad = (latitude*np.pi)/180
    delta = 0.4093*np.sin((step*2*np.pi/366)-1.405)  # Declination in radians
    # delta=((-23.5*np.pi)/180)*np.cos(360/366*(step+10))#Declination in radians
    dr = 1+0.033*np.cos(step*2*np.pi/366)  # Eccentricity factor
    # ; %sunset hour angle in radians
    ws = np.arccos(-np.tan(latRad)*np.tan(delta))
    N = (ws*24.)/np.pi  # day length in hours

    """
    Air temperature at sensor observation Parton & Logan (1981) (Eq. 3.2.7).
    Tmax--> Maximum pixel temperature of the day
    Tmin--> Minumum pixel temperature of the day
    m 
    N--> Elapsedtime between sunrise and sunset.Calculations are based on NOAA sunrise and sunset calculator
     using the pysunset scripts available in https://github.com/rconradharris/pysunset
     
    D
     """

    m = ObservationTime - (12-(N/2)+c2)
    Tair = ((Tmax-Tmin))*np.sin((np.pi*m)/(N+2*d))+Tmin
    Tair[Tair <= -1000] = np.nan
   # RLDInst=sigma*(((Tair)+273.15)**4)*(1-(c1*np.exp(-d*(Tair**2))))
    # RLDInst=sigma*(((Tair))**4)*(1-(c1*np.exp(-d*(Tair**2))))

    # Calculating incoming longwave radiation
    Stef = 5.67e-008
    C = 0.261
    d = 7.77e-004
    Tair_MODIS_passtime = Tair
    Rlw_in = Stef*(LST**4)*(1-(C*np.exp(-d*(Tair_MODIS_passtime**2))))

    return Rlw_in, Tair


def OutLongwaveRadiationInst(LST, Reflectance_b1, Reflectance_b2):
    """This function computes the longwave outcoming radiation at the overpass time
    of the sensor. It computes the NDVI and uses it to calculate the emissivity of the pixels
    """
    import numpy as np
    sigma = 5.67e-8  # Stefan Boltzman
    NDVI = (Reflectance_b2-Reflectance_b1)/(Reflectance_b2+Reflectance_b1)
    # Reescales the NDVI into emissivity assuming a linear relationship
    emissmax = 0.9999  # Vegetation emissivity
    emissmin = 0.9200  # Soil emissivity
    NDVI_max = 0.6  # Limit where is only vegetation
    NDVI_min = 0.2  # limint where is only soil
    emissivity = (emissmax-emissmin)/(NDVI_max-NDVI_min) * \
        (NDVI-NDVI_max)+emissmax
    emissivity[emissivity < emissmin] = emissmin
    emissivity[emissivity > emissmax] = emissmax

    emiss = 1.0094+0.047*np.log(NDVI)
    emiss[emiss > 0.986] = 0.986
    emiss[emiss < 0.914] = 0.914
    LongOut = -emissivity*sigma*LST**4
    LongOut2 = -emiss*sigma*LST**4

    return LongOut2


def shortWaveRadiationDown(RSwD, Albedo):
    RSnet = RSwD*(1-Albedo)
    RSout = RSwD-RSnet
    return RSout


def splitNetRadiationCanopySoil(Rn, LAI):
    kRn = 0.6  # García et al(2013)based on Impens&Lemeur (1969)
    # Rns is calculated as it follows in García et al (2013) which Fisher says is based in Beer (1852)
    RnS = Rn**(-kRn*LAI)  # (Wm-2)
    # Alternative: Rns can be calculated as follows in Norman et al (1995).Rns2=Rn_daily.*exp(0.9*log(1-fc))
    # NET RADIATION TO THE CANOPY
    RnC = Rn-RnS
    return RnC, RnS

# def calculateRLDownInstIDsoandJackson_Using_NOAA(Tmax,Tmin,latitude, longitude,year,step):
#    import noaa
#    import numpy as np
#    from jdcal import gcal2jd, jd2gcal
#    import datetime
#    sigma=5.67e-8 #Stefan Boltzman
#    c1=0.261 #constant
#    c2=-0.17 #Lag time in minimum air temperature
#    d=7.77e-4 #constant
#    ObservationTime=12.00 #Assumed to be this time until we get the information regarding the observation time
#    #Converts the year and step from the netcdf file to gregorian calendar so it can be use in the function from noaa.
#    jdcal=gcal2jd(year,01,01)
#    julday=jdcal[1]+step
#    date=jd2gcal(jdcal[0],julday)
#    date = datetime.datetime(date[0],date[1],date[2])
#    xdims=latitude.shape[0]
#    ydims=latitude.shape[1]
#    latitude=latitude.flatten()
#    longitude=longitude.flatten()
#    SunriseTime=np.zeros(latitude.shape)
#    SunsetTime=np.zeros(latitude.shape)
#    N_array=np.zeros(latitude.shape)
#
#    """
#    Air temperature at sensor observation Parton & Logan (1981) (Eq. 3.2.7).
#    Tmax--> Maximum pixel temperature of the day
#    Tmin--> Minumum pixel temperature of the day
#    m
#    N--> Elapsedtime between sunrise and sunset.Calculations are based on NOAA sunrise and sunset calculator
#     using the pysunset scripts available in https://github.com/rconradharris/pysunset
#
#    D
#     """
#    utc_offset=1
#
#    for i in xrange(0,np.int(latitude.shape[0])):
#        #for ii in xrange(0,np.int(latitude.shape[1])):
#
#        SunriseTime=noaa.get_sunrise(date, latitude[i], longitude[i], utc_offset)
#        SunsetTime=noaa.get_sunset(date, latitude[i], longitude[i], utc_offset)
#
#        N_temp=SunsetTime-SunriseTime
#        #N_temp=(SunriseTime)+(SunriseTime)-(SunriseTime)
#        N_array[i]=(N_temp.seconds)/3600.0
#    N_array=N_array.reshape(xdims,ydims)
#    m=ObservationTime - (12-(N_array/2)+c2)
#    Tair=(Tmax/10.-Tmin/10.)*np.sin((np.pi*m)/(N_array+2*d))+Tmin
#    Tair[Tair<=-1000]=np.nan
#    RLDInst=sigma*(((Tair)+273.15)**4)*(1-(c1*np.exp(-d*(Tair**2))))
#    return RLDInst,Tair,N_array
