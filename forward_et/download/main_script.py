# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 11:10:21 2017

@author: gmgo
"""

""" This script calculates the ET using the JPL-PT model using MODIS data. It uses the aproach from M.Garcia et al. 2013"""


import GetMODISdata
import datetime 
import Prepare_input_data
import Functions_PT_JPL
from Reproject_Subset_Tiles2Reference import getSceneMetadata  
reference_Image='N:\\DASEMON\\DEM\\SRTM_DASEMON_1100m.tif'
#reference_Image='C:\\Procesos\\Andalucia\\AndaluciaRef.tif'
reference_metadata=getSceneMetadata(reference_Image)
year2process=2008
base = datetime.datetime(year2process,01,01)
date_list = [base + datetime.timedelta(days=x) for x in range(0, 366,1)]
date_list_8days = [base + datetime.timedelta(days=x) for x in range(0, 366,8)]
tiles='h17v04,h17v05,h18v04,h18v05'
Test=-1
Test2=-1
MainFolder='F:\\PT-JPL_MODIS\\' #Folder containing all inputs


InDataset,RhDataset,TmaxDataset,TminDataset,precip=Prepare_input_data.ReadDASEMON_inputs()          
#execfile("DonwloadERA_INTERIM4PTJPL_MODIS.py")

for date in date_list:
    print date
    try:
        Test,Test2=GetMODISdata.DownloadAllData(tiles,date,date_list_8days,Test,Test2,MainFolder)
        JulianDay,JulDay8=Functions_PT_JPL.JulDay_JulDay8fromDate(date)
        JulDay16=Functions_PT_JPL.JulDay_JulDay16fromDate(date)
        Dictionary=Prepare_input_data.MakeInputDictionary()
    
        Dictionary=Prepare_input_data.FillReflectanceVariables(MainFolder,Dictionary,JulianDay,year2process,reference_metadata)
        NDVI16=Prepare_input_data.FillNDVI_16_Days(MainFolder,Dictionary,JulDay16,year2process,reference_metadata)
        Dictionary=Prepare_input_data.FillLSTVariables(MainFolder,Dictionary,JulianDay,year2process,reference_metadata)
        Dictionary=Prepare_input_data.FillAlbedoInputs(MainFolder,Dictionary,JulDay8,year2process,reference_metadata)
        Dictionary=Prepare_input_data.FillAPAR_and_LAI_Inputs(MainFolder,Dictionary,JulDay8,year2process,reference_metadata)
        Dictionary=Prepare_input_data.fillECMWF_Variables(MainFolder,Dictionary,JulianDay,year2process)
 #       Dictionary=Prepare_input_data.FillDASEMON_inputs(InDataset,RhDataset,TmaxDataset,TminDataset,precip,Dictionary,JulianDay,reference_metadata)
        
        Dictionary ['Year']= year2process
        Dictionary ['JulianDay']= JulianDay
        import pickle
        DictionaryFolder='F:\\caca\\'
        JulDay=int(JulianDay)
        if int(JulDay) <100:
            if int(JulDay) <10:
                strJulDay='00'+str(JulDay)
            else:
                strJulDay='0'+str(JulDay)
        else:
            strJulDay=str(JulDay)
        Dictionary= Prepare_input_data.CleanDictionary4Students(Dictionary)
        Name_dictionary=DictionaryFolder+'DailyDataset_julianDay_'+strJulDay+'_Year_'+str(year2process)+'.pkl'
        f=open(Name_dictionary,'wb')
        pickle.dump(Dictionary,f)
        f.close()
        print (Name_dictionary+' has been created')
    except:
        year=str(date.year)
        month=date.month
        if month < 10:
            monthstr='0'+str(month)
        else:
            monthstr=str(month)
            
        day=(date.day)
        if day < 10:
            daystr='0'+str(day)
        else:
            daystr=str(day)
        day_full=year+'-'+monthstr+'-'+daystr
        JulianDay,JulDay8=Functions_PT_JPL.JulDay_JulDay8fromDate(date)
        day_full=day_full+'_Julday'+str(JulianDay)

