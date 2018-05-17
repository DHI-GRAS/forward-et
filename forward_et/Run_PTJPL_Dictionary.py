# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 14:00:34 2018

@author: gmgo
"""

"""
This script calculates Evapotranspiration (ET) based on the PT-JPL algorithm. 
The input dataset to obtain ET was created on a combination of two different datasets, DASEMON and MODIS.
This script reads a daily input dictionary that contains all required inputs.

Code structure and ancillary functions:


"""
"""
We have created different ancillary scrips that contain functions to calculate the radiation and evapotrasnpiration.
To access those functions we need to import the modules.
import Functions
import Radiation_Functions
To make the name shorter in the code when we call the functions we can import the modules followed by "as" :
    
import Functions as fn
import Radiation_Functions as rf 

This way we can called the functions in the scripts by typing:

Example:
    N,Ra,Rs,ws=rf.CalShortWaveIncomingRadiation(SunHours,Latitude,Doy) instead of
    
    N,Ra,Rs,ws=Radiation_Functions.CalShortWaveIncomingRadiation(SunHours,Latitude,Doy)
"""
import Functions as fn
import Radiation_Functions as rf


"""
The script starts getting a list of the files that are going to be proccesed. To do that we use the package glob
that has to be imported.
"""
import glob
"""
In the next line we get a list of the pickle files that contains the inputs to calculate ET. It is only needed to indicate
the folder that contains all the files.
In the example here are in F:\Dictionaries
Paths can be expressed as:
    F:\\Dictionaries\\ or F:/Dictionaries/
The we just need to specify the file extention .pkl
And the line gets something like this.
FileList=glob.glob(InputsFolder+'*.pkl')
"""
InputsFolder = 'F:\\Dictionaries\\'
FileList = glob.glob(InputsFolder+'*.pkl')
for inputfile in FileList[75:76]:
    Dictionary = fn.OpenDictionay(inputfile)
    """
    STARTS CALCULATING THE SHORTWAVE COMPONENTS OF THE RADIATION
    In the first part of the code the inputs necesary to calculate the radiation components are read from the dictionary file
    """
    SunHours, Latitude, Doy, Albedo = fn.GetShortWaveRadiationInputs(
        Dictionary)

    """
    
    ****************************************
    SHORTWAVE INCOMING
    
    
    Please read the pages 41 to 51 from the FAO Irrigation and Drainage paper No. 56.
    You can get the pdf document in this website: http://academic.uprm.edu/abe/backup2/tomas/fao%2056.pdf
    In those pages there is a detail description of how to calculate the SWIn raditation and the equations that
    are implemented in the funtions below.
    ****************************************
    """

    N, Ra, Rs, ws = rf.CalShortWaveIncomingRadiation(SunHours, Latitude, Doy)

    """
    ****************************************
    SHORTWAVE NET
    ****************************************
    """

    RSnet = rf.CalShortWaveNetRadiation(Rs, Albedo)

    """
    ****************************************
    LONGWAVE INCOMING
    Visit the paper of Partorn and Logan from 1981 to understand how to convert land surface temperature to air temperature
    It is necesary to calculate the minimum and maximum airtemperature in order to get a good estimate of the air temperature
    at MODIS overpass time. We've done this step already for you, and are included as two numpy arrays called:
        YearMaxTemp.npy  and  YearMinTemp.npy
    Both arrays contain the temperature in Celsius degress. Therefore they must be converted to Kelvin.
    
    
    ****************************************
    """
    LST, ObsTime, Emissivity = fn.GetLongWaveRadiationInputs(Dictionary)

    Rlw_in, Tair_MODIS_passtime = rf.CalLongWave_Incoming_Parton_Logan(
        LST, ObsTime, N)
    """
    ****************************************
    LONGWAVE OUTGOING
    The longwave outgoing radiation uses the mean values of emissivity from chanels 31 and 32 from MODIS.
    ****************************************
    """
    Rlw_out = rf.CalLongWave_Outgoing(LST, Emissivity)
    """
    ****************************************
    LONGWAVE NET RADIATION
    ****************************************
    """
    RLnet = Rlw_in-Rlw_out

    """
    ******************************************************************************************************************************************
    ******************************************************************************************************************************************
    ******************************************************************************************************************************************
    NOTICE THAT LONGWAVE INCOMING RADIATION IS INSTANTANEOUS NOT DAILY AND DEREFORE IT HAS TO BE CONVERTED TO GET EVERYTHING IN THE SAME UNITS
    ******************************************************************************************************************************************
    ******************************************************************************************************************************************
    ******************************************************************************************************************************************    
    """

    """
    Calculates the J parameter to estimate daily/instantaneous from observations.
    """
    t = 12-(N/2)-1  # We assumed this formulation for sunrise time. Can be modified and improved
    J = fn.calcJparameter(N, t)
    """
    ****************************************
    LONGWAVE NET RADIATION DAILY
    ****************************************
    """
    RLnet_daily = RLnet*J

    """
    ****************************************
    CALCULATES NET RADIATION DAILY
    ****************************************
    """

    RnDaily = RSnet+RLnet_daily
