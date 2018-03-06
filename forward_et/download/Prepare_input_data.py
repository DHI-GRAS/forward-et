# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 09:44:29 2018

@author: gmgo
"""

def fillECMWF_Variables(MainFolder,Dictionary,JulDay,year):
    import gdal
    import numpy as np
    if JulDay<100:
        if JulDay<10:
            StrJul='00'+str(int(JulDay))
        else:
            StrJul='0'+str(int(JulDay))
    else:
        StrJul=str(int(JulDay))
        
    SSRD_file=MainFolder+'\\SSRD\\SSRD_'+str(year)+StrJul+'_ECMWF.tif'
    fid=gdal.Open(SSRD_file,gdal.GA_ReadOnly)
    SSRD=fid.GetRasterBand(1).ReadAsArray()
    Dictionary ['ECMWF_SSRD']= SSRD 
    
    TP_file=MainFolder+'\\TP\\TP_'+str(year)+StrJul+'_ECMWF.tif'
    fid=gdal.Open(TP_file,gdal.GA_ReadOnly)
    TP=fid.GetRasterBand(1).ReadAsArray()
    Dictionary ['ECMWF_Precipitation']= TP
    
    airTemp9_file=MainFolder+'\\2T\\2T_'+str(year)+StrJul+'_09_ECMWF.tif'
    fid=gdal.Open(airTemp9_file,gdal.GA_ReadOnly)
    airTemp9=fid.GetRasterBand(1).ReadAsArray()

    
    airTemp12_file=MainFolder+'\\2T\\2T_'+str(year)+StrJul+'_12_ECMWF.tif'
    fid=gdal.Open(airTemp12_file,gdal.GA_ReadOnly)
    airTemp12=fid.GetRasterBand(1).ReadAsArray()

    
    airTemp15_file=MainFolder+'\\2T\\2T_'+str(year)+StrJul+'_15_ECMWF.tif'
    fid=gdal.Open(airTemp15_file,gdal.GA_ReadOnly)
    airTemp15=fid.GetRasterBand(1).ReadAsArray()

    Airtemp=np.zeros([airTemp15.shape[0],airTemp15.shape[1]])
    Obstime=Dictionary ['LST_Composite_Time_day']
    Airtemp[Obstime<10.5]=airTemp9[Obstime<10.5]
    Airtemp[(Obstime>10.5) & (Obstime<13.5)]=airTemp12[(Obstime>10.5) & (Obstime<13.5)]
    Airtemp[(Obstime>13.5) & (Obstime<15)]=airTemp15[(Obstime>13.5) & (Obstime<15)]
    Dictionary ['ECMWF_Airtemp']= Airtemp
    
    return Dictionary
def cleanDirectory(folder):
    import os
    for root, dirs, files in os.walk(folder, topdown = False):
       for file in files:
#          print(os.path.join(root, file))
          os.remove(os.path.join(root, file))
          
    return
def MakeInputDictionary():
    Dictionary              =            {}
    Dictionary ['Year']= []
    Dictionary ['JulianDay']= []
    Dictionary ['Latitude']= [] 
    Dictionary ['Longitude']= [] 
    Dictionary ['SunHours']= [] 
    Dictionary ['Precipitation']= [] 
    Dictionary ['ECMWF_Precipitation']= [] 
    Dictionary ['ECMWF_SSRD']= [] 
    Dictionary ['ECMWF_Airtemp']= [] 
    Dictionary ['DASEMON_Tmin']= [] 
    Dictionary ['DASEMON_Tmax']= [] 
    Dictionary ['Reflectance_B1_Terra']= [] 
    Dictionary ['Reflectance_B2_Terra']= [] 
    Dictionary ['Reflectance_B1_Aqua']=  [] 
    Dictionary ['Reflectance_B2_Aqua']=  []
    Dictionary ['Reflectance_Aqua_Mask']=  []
    Dictionary ['Reflectance_Terra_Mask']=  []
    Dictionary ['NDVI']                 =  []
    Dictionary ['NDVI_16days']                 =  []
    Dictionary ['SAVI']                 =  []
    Dictionary ['fPAR']        =         []
    Dictionary ['LAI']        =         []
    Dictionary ['fPAR_mask']        =         []
    Dictionary ['emis_31_Aqua']   =      []
    Dictionary ['emis_32_Aqua']   =      []
    Dictionary ['emis_31_Terra']   =     []
    Dictionary ['emis_32_Terra']   =     []
    Dictionary ['emis_Composite']   =     []
    Dictionary ['LST_Terra']      =      []
    Dictionary ['LST_Aqua']      =       []
    Dictionary ['LST_Terra_Time_day'] =  []
    Dictionary ['LST_Aqua_Time_day'] =   []
    Dictionary ['LST_Aqua_Mask']=        []
    Dictionary ['LST_Terra_Mask']=       []
    Dictionary ['VZA_Terra_Mask']=       []
    Dictionary ['VZA_Terra_Mask']=       []
    Dictionary ['Albedo_WSA']        =   []  
    Dictionary ['Albedo_BSA']        =   []  
    Dictionary ['Albedo']           =   []  
    Dictionary ['Albedo_Mask']        =  []  
    Dictionary ['emis_31_Composite']   =  [] 
    Dictionary ['emis_32_Composite']   =     [] 
    Dictionary ['LST_Composite']      =      [] 
    Dictionary ['LST_Composite_Time_day'] =  [] 
    Dictionary ['VZA_Composite']=       [] 
    return Dictionary
    
def CleanDictionary4Students(Dictionary):
    del Dictionary ['Reflectance_B1_Terra']
    del Dictionary ['Reflectance_B2_Terra']
    del Dictionary ['Reflectance_B1_Aqua']
    del Dictionary ['Reflectance_B2_Aqua']
    del Dictionary ['Reflectance_Aqua_Mask']
    del Dictionary ['Reflectance_Terra_Mask']
    del Dictionary ['emis_31_Aqua']   
    del Dictionary ['emis_32_Aqua']   
    del Dictionary ['emis_31_Terra']   
    del Dictionary ['emis_32_Terra']  
    del Dictionary ['fPAR_mask']  
    del Dictionary ['Albedo_WSA']      
    del Dictionary ['Albedo_BSA'] 
    del Dictionary ['LST_Terra_Time_day'] 
    del Dictionary ['LST_Aqua_Time_day']  
    del Dictionary ['VZA_Terra']  
    del Dictionary ['VZA_Aqua']  
    del Dictionary ['VZA_Terra_Mask']  
    del Dictionary ['Albedo_Mask']  
    del Dictionary ['LST_Terra_Mask']  
    del Dictionary ['LST_Aqua_Mask'] 
    del Dictionary ['LST_Aqua']  
    del Dictionary ['LST_Terra'] 
    del Dictionary ['emis_31_Composite']  
    del Dictionary ['emis_32_Composite']  
    del Dictionary ['NDVI']  
    del Dictionary ['SAVI'] 
    del Dictionary ['VZA_Composite'] 
    del Dictionary ['Longitude'] 
    del Dictionary ['Latitude'] 
    del Dictionary ['NDVI_16days']
    del Dictionary ['Precipitation']
    return Dictionary
def FillNDVI_16_Days(MainFolder,Dictionary,JulDay,year,reference_metadata):
    import glob
    import gdal
    from pymasker import Masker
    import numpy as np
    import os
    from Reproject_Subset_Tiles2Reference import reprojectSubsettif  
    
    
    JulDay=int(JulDay)
    if int(JulDay) <100:
        if int(JulDay) <10:
            strJulDay='00'+str(JulDay)
        else:
            strJulDay='0'+str(JulDay)
    else:
        strJulDay=str(JulDay)

    Folder=MainFolder+'MOD_13\\'
    Folder_Subset=MainFolder+'MOD_13\\Subset\\'
    
    NDVI16=glob.glob(Folder+'MOD13A2.006'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( NDVI16[0], Folder_Subset,reference_metadata)
    NDVI16=glob.glob(Folder_Subset+'MOD13A2.006'+str(year)+strJulDay+'.mosaic.tif')
    fid_NDVI16=gdal.Open(NDVI16[0],gdal.GA_ReadOnly)
    StateQA_NDVI16=fid_NDVI16.GetRasterBand(2).ReadAsArray()
    NDVI16=fid_NDVI16.GetRasterBand(1).ReadAsArray()
#    NDVI16[StateQA_NDVI16>=2]=0 #Comemnted to see the resutls
#    NDVI16[StateQA_NDVI16<0]=0 #Comemnted to see the resutls
    Dictionary['NDVI_16']=NDVI16
       
    return Dictionary
def FillReflectanceVariables(MainFolder,Dictionary,JulDay,year,reference_metadata):
    import glob
    import gdal
    from pymasker import Masker
    import numpy as np
    import os
    from Reproject_Subset_Tiles2Reference import reprojectSubsettif  
    
    
    JulDay=int(JulDay)
    if int(JulDay) <100:
        if int(JulDay) <10:
            strJulDay='00'+str(JulDay)
        else:
            strJulDay='0'+str(JulDay)
    else:
        strJulDay=str(JulDay)

    Folder=MainFolder+'MOD_MYD_09\\'
    Folder_Subset=MainFolder+'MOD_MYD_09\\Subset\\'
    
    Terra=glob.glob(Folder+'MOD09GA.006'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( Terra[0], Folder_Subset,reference_metadata)
    Terra=glob.glob(Folder_Subset+'MOD09GA.006'+str(year)+strJulDay+'.mosaic.tif')
    fid_Terra=gdal.Open(Terra[0],gdal.GA_ReadOnly)
    StateQA_Terra=fid_Terra.GetRasterBand(1).ReadAsArray()
    masker=Masker(StateQA_Terra)
    CloudState_Terra=masker.get_mask(0,2,'00')
    CloudShadow_Terra=masker.get_mask(2,1,'0')
    LandWater_Terra=masker.get_mask(3,3,'001') #Land
    LandWater_Terra_inlandWater=masker.get_mask(3,3,'011') #Shallow inland
    
    
    SensorZenith_Terra=fid_Terra.GetRasterBand(2).ReadAsArray()
    SensorZenith_Terra=SensorZenith_Terra*1.0
    SensorZenith_Terra[SensorZenith_Terra>4500]=np.nan
    B1_Terra=fid_Terra.GetRasterBand(3).ReadAsArray()
    B2_Terra=fid_Terra.GetRasterBand(4).ReadAsArray()
    NDVI_Terra=(B2_Terra-B1_Terra)*1.00/(B2_Terra+B1_Terra)
    SAVI_Terra=1.5*((B2_Terra-B1_Terra)/(B2_Terra+B1_Terra+0.5))
    
    Aqua=glob.glob(Folder+'MYD09GA.006'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( Aqua[0], Folder_Subset,reference_metadata)
    Aqua=glob.glob(Folder_Subset+'MOD09GA.006'+str(year)+strJulDay+'.mosaic.tif')
    fid_Aqua=gdal.Open(Aqua[0],gdal.GA_ReadOnly)
    StateQA_Aqua=fid_Aqua.GetRasterBand(1).ReadAsArray()
    maskerAqua=Masker(StateQA_Aqua)
    CloudState_Aqua=maskerAqua.get_mask(0,2,'00')
    CloudShadow_Aqua=maskerAqua.get_mask(2,1,'0')
    LandWater_Aqua=masker.get_mask(3,3,'001')
    LandWater_Aqua_inlandWater=masker.get_mask(3,3,'011') #Shallow inland
    
    SensorZenith_Aqua=fid_Aqua.GetRasterBand(2).ReadAsArray()
    SensorZenith_Aqua=SensorZenith_Aqua*1.0
    SensorZenith_Aqua[SensorZenith_Aqua>4500]=np.nan
    B1_Aqua=fid_Aqua.GetRasterBand(3).ReadAsArray()
    B2_Aqua=fid_Aqua.GetRasterBand(4).ReadAsArray()
    NDVI_Aqua=(B2_Aqua-B1_Aqua)*1.00/(B2_Aqua+B1_Aqua)
    SAVI_Aqua=1.5*((B2_Aqua-B1_Aqua)/(B2_Aqua+B1_Aqua+0.5))
    
    NDVI_Terra[np.isnan(SensorZenith_Terra)]=np.nan
    NDVI_Aqua[np.isnan(SensorZenith_Aqua)]=np.nan
    NDVI=np.array([NDVI_Terra,NDVI_Aqua])
    NDVI=np.nanmean(NDVI,axis=0)
    
    SAVI_Terra[np.isnan(SensorZenith_Terra)]=np.nan
    SAVI_Aqua[np.isnan(SensorZenith_Aqua)]=np.nan
    SAVI=np.array([SAVI_Terra,SAVI_Aqua])
    SAVI=np.nanmean(SAVI,axis=0)
    SAVI[SAVI<-50000]=np.nan
    CloudState=CloudState_Aqua+CloudState_Terra
    CloudShadow=CloudShadow_Terra+CloudShadow_Aqua
    LandWater=LandWater_Aqua+LandWater_Terra+LandWater_Aqua_inlandWater+LandWater_Terra_inlandWater
    NDVI[CloudState==0]=np.nan
    NDVI[CloudShadow==0]=np.nan
    NDVI[LandWater==0]=np.nan
    NDVI[NDVI>1]=np.nan
    
    NDVI[np.isnan(NDVI)]=0
    NDVI=NDVI*1000
    NDVI=np.int16(NDVI)
    
    SAVI[np.isnan(SAVI)]=0
    SAVI=SAVI*1000
    SAVI=np.int16(SAVI)
    Dictionary['NDVI']=NDVI
    Dictionary['SAVI']=SAVI
    Dictionary ['Reflectance_B1_Terra']= B1_Terra 
    Dictionary ['Reflectance_B2_Terra']= B2_Terra 
    Dictionary ['Reflectance_B1_Aqua']=  B1_Aqua
    Dictionary ['Reflectance_B2_Aqua']=  B2_Aqua
    Dictionary ['Reflectance_Aqua_Mask']=  []
    Dictionary ['Reflectance_Terra_Mask']=  []
    return Dictionary
    
    
def FillLSTVariables(MainFolder,Dictionary,JulDay,year,reference_metadata):
    import glob
    import gdal
    from pymasker import Masker
    import numpy as np
    from Reproject_Subset_Tiles2Reference import reprojectSubsettif  
    
    JulDay=int(JulDay)
    if int(JulDay) <100:
        if int(JulDay) <10:
            strJulDay='00'+str(JulDay)
        else:
            strJulDay='0'+str(JulDay)
    else:
        strJulDay=str(JulDay)

    Folder=MainFolder+'MOD_MYD_11\\'
    Folder_Subset=Folder+'Subset\\'
    
    Terra=glob.glob(Folder+'MOD11A1.006'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( Terra[0], Folder_Subset,reference_metadata)
    Terra=glob.glob(Folder_Subset+'MOD11A1.006'+str(year)+strJulDay+'.mosaic.tif')
    fid_Terra=gdal.Open(Terra[0],gdal.GA_ReadOnly)
    LST_day_Terra=fid_Terra.GetRasterBand(1).ReadAsArray()
    QC_LST_Terra=fid_Terra.GetRasterBand(2).ReadAsArray()
    DayObsTime_Terra=fid_Terra.GetRasterBand(3).ReadAsArray()
    DayObsAngle_Terra=fid_Terra.GetRasterBand(4).ReadAsArray()
    DayObsAngle_Terra=np.int16(np.abs(DayObsAngle_Terra-65.0))
    Band31Emis_Terra=fid_Terra.GetRasterBand(9).ReadAsArray()
    Band32Emis_Terra=fid_Terra.GetRasterBand(10).ReadAsArray()
    
    masker=Masker(QC_LST_Terra)
    LST_Error_Terra=masker.get_mask(4,2,'00') #less 1 kelvin
    LST_Error_Terra2=masker.get_mask(4,2,'01') #less 1 kelvin
    LST_Error_Terra3=masker.get_mask(6,2,'10') #less 3 kelvin
    LST_Error_Terra=LST_Error_Terra+LST_Error_Terra2+LST_Error_Terra3
    LST_Error_Terra[LST_Error_Terra>=1]=1
    Aqua=glob.glob(Folder+'MYD11A1.006'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( Aqua[0], Folder_Subset,reference_metadata)
    Aqua=glob.glob(Folder_Subset+'MYD11A1.006'+str(year)+strJulDay+'.mosaic.tif')
    fid_Aqua=gdal.Open(Aqua[0],gdal.GA_ReadOnly)
    LST_day_Aqua=fid_Aqua.GetRasterBand(1).ReadAsArray()
    QC_LST_Aqua=fid_Aqua.GetRasterBand(2).ReadAsArray()
    DayObsTime_Aqua=fid_Aqua.GetRasterBand(3).ReadAsArray()
    DayObsAngle_Aqua=fid_Aqua.GetRasterBand(4).ReadAsArray()
    DayObsAngle_Aqua=np.int16(np.abs(DayObsAngle_Aqua-65.0))
    Band31Emis_Aqua=fid_Aqua.GetRasterBand(9).ReadAsArray()
    Band32Emis_Aqua=fid_Aqua.GetRasterBand(10).ReadAsArray()
    
    masker=Masker(QC_LST_Aqua)
    LST_Error_Aqua=masker.get_mask(4,2,'00') #less 1 kelvin
    LST_Error_Aqua2=masker.get_mask(6,2,'01') #less 3 kelvin
    LST_Error_Aqua3=masker.get_mask(6,2,'10') #less 3 kelvin
    LST_Error_Aqua=LST_Error_Aqua+LST_Error_Aqua2+LST_Error_Aqua3
    LST_Error_Aqua[LST_Error_Aqua>=1]=2
    LST_Mask=LST_Error_Terra+LST_Error_Aqua
    #Creates the composites of terra and aqua to run the PT-JPL 
    
    LST_Composite=np.zeros([LST_day_Aqua.shape[0],LST_day_Aqua.shape[1]])
    ObsTime_Composite=np.zeros([LST_day_Aqua.shape[0],LST_day_Aqua.shape[1]])
    ObsAngle_Composite=np.zeros([LST_day_Aqua.shape[0],LST_day_Aqua.shape[1]])
    Emiss31Composite=np.zeros([LST_day_Aqua.shape[0],LST_day_Aqua.shape[1]])
    Emiss32Composite=np.zeros([LST_day_Aqua.shape[0],LST_day_Aqua.shape[1]])
    
    for i in range(0,LST_day_Aqua.shape[0]):
        for ii in range(0,LST_day_Aqua.shape[1]):
            if LST_Mask[i,ii]==1:
                LST_Composite[i,ii]=LST_day_Terra[i,ii]
                ObsTime_Composite[i,ii]=DayObsTime_Terra[i,ii]
                ObsAngle_Composite[i,ii]=DayObsAngle_Terra[i,ii]
                Emiss31Composite[i,ii]=Band31Emis_Terra[i,ii]
                Emiss32Composite[i,ii]=Band32Emis_Terra[i,ii]
                
            if LST_Mask[i,ii]==2:
                LST_Composite[i,ii]=LST_day_Aqua[i,ii]
                ObsTime_Composite[i,ii]=DayObsTime_Aqua[i,ii]
                ObsAngle_Composite[i,ii]=DayObsAngle_Aqua[i,ii]
                Emiss31Composite[i,ii]=Band31Emis_Aqua[i,ii]
                Emiss32Composite[i,ii]=Band32Emis_Aqua[i,ii]
                
            if LST_Mask[i,ii]==3:
                if DayObsAngle_Aqua[i,ii]<DayObsTime_Terra[i,ii]:
                    LST_Composite[i,ii]=LST_day_Aqua[i,ii]
                    ObsTime_Composite[i,ii]=DayObsTime_Aqua[i,ii]
                    ObsAngle_Composite[i,ii]=DayObsAngle_Aqua[i,ii]
                    Emiss31Composite[i,ii]=Band31Emis_Aqua[i,ii]
                    Emiss32Composite[i,ii]=Band32Emis_Aqua[i,ii]
                else:
                    LST_Composite[i,ii]=LST_day_Terra[i,ii]
                    ObsTime_Composite[i,ii]=DayObsTime_Terra[i,ii]
                    ObsAngle_Composite[i,ii]=DayObsAngle_Terra[i,ii]
                    Emiss31Composite[i,ii]=Band31Emis_Terra[i,ii]
                    Emiss32Composite[i,ii]=Band32Emis_Terra[i,ii]
                    
    LST_Composite=LST_Composite*0.02 #in kelvin
    LST_Composite=np.int16(LST_Composite)  

    ObsTime_Composite=ObsTime_Composite*0.1
    ObsTime_Composite=np.int16(ObsTime_Composite)
    
    Emiss31Composite=Emiss31Composite*0.002
    Emiss32Composite=Emiss32Composite*0.002
    EmissComposite=((Emiss31Composite+Emiss32Composite)/2)*100
    EmissComposite  =np.int16(EmissComposite)
    Dictionary ['emis_31_Aqua']   =  Band31Emis_Aqua
    Dictionary ['emis_32_Aqua']   =  Band32Emis_Aqua
    Dictionary ['emis_31_Terra']   =     Band31Emis_Terra
    Dictionary ['emis_32_Terra']   =     Band32Emis_Terra
    Dictionary ['LST_Terra']      =      LST_day_Terra
    Dictionary ['LST_Aqua']      =       LST_day_Aqua
    Dictionary ['LST_Terra_Time_day'] =  DayObsTime_Terra
    Dictionary ['LST_Aqua_Time_day'] =   DayObsTime_Aqua
    Dictionary ['LST_Aqua_Mask']=        LST_Error_Aqua
    Dictionary ['LST_Terra_Mask']=       LST_Error_Terra
    Dictionary ['VZA_Terra']=       DayObsAngle_Aqua
    Dictionary ['VZA_Aqua']=       DayObsAngle_Terra
    
    Dictionary ['emis_31_Composite']   =  Emiss31Composite
    Dictionary ['emis_32_Composite']   =     Emiss32Composite
    Dictionary ['emis_Composite']   =     EmissComposite
    Dictionary ['LST_Composite']      =      LST_Composite
    Dictionary ['LST_Composite_Time_day'] =  ObsTime_Composite
    Dictionary ['VZA_Composite']=       ObsAngle_Composite
    
    
    return Dictionary

def FillAlbedoInputs(MainFolder,Dictionary,JulDay,year,reference_metadata):
    import glob
    import gdal
    import numpy as np
    from Reproject_Subset_Tiles2Reference import reprojectSubsettif  
    JulDay=int(JulDay)
    if int(JulDay) <100:
        if int(JulDay) <10:
            strJulDay='00'+str(JulDay)
        else:
            strJulDay='0'+str(JulDay)
    else:
        strJulDay=str(JulDay)

    Folder=MainFolder+'MCD_43\\'
    Folder_Subset=Folder+'Subset\\'
    
    AlbedoImage=glob.glob(Folder+'MCD43B3.005'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( AlbedoImage[0], Folder_Subset,reference_metadata)
    
    AlbedoImage=glob.glob(Folder_Subset+'MCD43B3.005'+str(year)+strJulDay+'.mosaic.tif')
    fid_AlbedoImage=gdal.Open(AlbedoImage[0],gdal.GA_ReadOnly)
    BSA=fid_AlbedoImage.GetRasterBand(1).ReadAsArray()
    WSA=fid_AlbedoImage.GetRasterBand(2).ReadAsArray()
    Albedo=(0.8*np.float64(BSA*0.001)+0.2*np.float64(WSA*0.001))
    
    AlbedoQuality=glob.glob(Folder+'MCD43B2.005'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( AlbedoQuality[0], Folder_Subset,reference_metadata)
    AlbedoQuality=glob.glob(Folder_Subset+'MCD43B2.005'+str(year)+strJulDay+'.mosaic.tif')
    fid_AlbedoQuality=gdal.Open(AlbedoQuality[0],gdal.GA_ReadOnly)
    Quality=fid_AlbedoQuality.GetRasterBand(1).ReadAsArray()
    Albedo[Quality!=0]=np.nan #Removed to see how it look
    Albedo[Albedo>1]=np.nan#Removed to see how it look
    Dictionary ['Albedo_WSA']        =   np.float64(WSA*0.001)  
    Dictionary ['Albedo_BSA']        =   np.float64(BSA*0.001) 
    Dictionary ['Albedo_Mask']        =  Quality # Zero are good data 
    Dictionary['Albedo']=   np.int16(Albedo*100)
    return Dictionary
    
def FillAPAR_and_LAI_Inputs(MainFolder,Dictionary,JulDay,year,reference_metadata):
    import glob
    import gdal
    import numpy as np
    from pymasker import Masker
    from Reproject_Subset_Tiles2Reference import reprojectSubsettif   
    JulDay=int(JulDay)
    if int(JulDay) <100:
        if int(JulDay) <10:
            strJulDay='00'+str(JulDay)
        else:
            strJulDay='0'+str(JulDay)
    else:
        strJulDay=str(JulDay)

    Folder=MainFolder+'MCD_15\\'
    Folder_Subset=Folder+'Subset\\'
    
    APAR_LAI=glob.glob(Folder+'MCD15A2H.006'+str(year)+strJulDay+'.mosaic.tif')
    reprojectSubsettif( APAR_LAI[0], Folder_Subset,reference_metadata)
    APAR_LAI=glob.glob(Folder_Subset+'MCD15A2H.006'+str(year)+strJulDay+'.mosaic.tif')
    fid_APARLAI=gdal.Open(APAR_LAI[0],gdal.GA_ReadOnly)
    LAI=fid_APARLAI.GetRasterBand(1).ReadAsArray()
    LAI=LAI
    APAR=fid_APARLAI.GetRasterBand(2).ReadAsArray()
    APAR=APAR
    QC=fid_APARLAI.GetRasterBand(3).ReadAsArray()
    masker=Masker(QC)
    MaskLAIAPAR=masker.get_mask(5,3,'000') #, Main (RT) method used, best result possible (no saturatio
    MaskLAIAPAR1=masker.get_mask(5,3,'001') #, Main (RT) method used with saturation. Good,very usable
    MaskLAIAPAR2=masker.get_mask(5,3,'010') #, Main (RT) method failed due to bad geometry, empirical algorithm used
    MaskLAIAPAR3=masker.get_mask(5,3,'011') #,  Main (RT) method failed due to problems other than geometry, empirical algorithm used
    MaskLAIAPAR=MaskLAIAPAR+MaskLAIAPAR1+MaskLAIAPAR2+MaskLAIAPAR3  
    MaskLAIAPAR[MaskLAIAPAR>=1]=1
#    LAI[MaskLAIAPAR!=1]=0  
#    APAR[MaskLAIAPAR!=1]=0 
    Dictionary ['fPAR']        =   np.int16( APAR)
    Dictionary ['LAI']        =   np.int16(LAI) 
    
    return Dictionary
def ReadDASEMON_inputs():  
    from netCDF4 import Dataset
    #Path where the files are contained
    folder='C:/Users/gmgo/Dropbox/Gorka/FORWARD/datasets/DESEMON/Vicente-Serrano/download_march7/'
    In_file                    =folder+'in.nc'
    Tmax_file                  =folder+'tmax.nc'
    Rh_file                    =folder+'rh.nc'
    Tmin_file                  =folder+'tmin.nc'
    precip                     =folder+'/good/pr_2008.nc'
    
    InDataset                    =Dataset(In_file, mode='r')
    RhDataset                    =Dataset(Rh_file, mode='r')
    TmaxDataset                  =Dataset(Tmax_file, mode='r')
    RhDataset                    =Dataset(Rh_file, mode='r')
    TminDataset                  =Dataset(Tmin_file, mode='r')
    precip                      =Dataset(precip, mode='r')
    return InDataset,RhDataset,TmaxDataset,TminDataset,precip
    
def deleteFile(path_to_file, attempts=0, timeout=100, sleep_int=2):
    import os
    from time import sleep
    if attempts < timeout and os.path.exists(path_to_file) and os.path.isfile(path_to_file): 
        try:
            os.remove(path_to_file)
            return True
        except:
            # perform an action
            sleep(sleep_int)
            deleteFile(path_to_file, attempts + 1)
    return False

def reprojectDASEMON(fileName, bandNum,Temp):
    import gdal
    from osgeo import gdal
    from osgeo.gdalnumeric import *
    import os
    import subprocess
    import re
    import gc
    import numpy
    from gdalconst import *
    import gc
    import subprocess
    from gdalconst import *
    reference_Image='N:\\DASEMON\\DEM\\SRTM_DASEMON_1100m.tif'
    ds = gdal.Open(reference_Image, gdal.GA_ReadOnly )
    dataOut=Temp
    # Constants
    C_srcProj = '+proj=longlat +a=6367470 +b=6367470 +no_defs' # projection of the ERA interim files
    C_gdal_translate = r'C:\WinPython\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\Lib\site-packages\osgeo\gdal_translate' 
    C_gdalwarp = r'C:\WinPython\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\Lib\site-packages\osgeo\gdalwarp'
    C_gdal_calc = r'C:\WinPython\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\Lib\site-packages\osgeo\gdal_calc.py' 

    dstProj = '+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs' 
    srcProj='+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs' 
    dstTr = '1000.00 -1000.00'
    dstTe = '100010.8359999999956926,3988484.6419999999925494 : 621010.8360000000102445,4288484.6419999999925494'
    #Write the out file  
    driver = gdal.GetDriverByName("GTiff") 
    outFileName = fileName
    dsOut = driver.Create(outFileName, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)  
    CopyDatasetInfo(ds,dsOut)  
    dsOut.GetRasterBand(1).WriteArray(Temp)  
    tiff_filename='C:\\Temporal_Trash\\temporal.tif'
    gdalWarpCmd = C_gdalwarp+' -wm 2000 -s_srs "'+srcProj+'" -t_srs "'+dstProj+'" -tr '+dstTr+' -te '+dstTe+' -r near -ot UInt32 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -overwrite '+outFileName+' '+tiff_filename
    proc=subprocess.Popen(gdalWarpCmd,shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.STDOUT, universal_newlines=False)
    for line in iter(proc.stdout.readline, ""):
        print line
    proc.wait()
       
    gc.collect()
    deleteFile(outFileName)
      
    #Close the datasets  
    
    ds = None  
      
    dsOut = None
    
    return 
def FillDASEMON_inputs(InDataset,RhDataset,TmaxDataset,TminDataset,precip,Dictionary,JulianDay,reference_Image):
    import numpy as np
    from pyproj import Proj, transform
    from Reproject_Subset_Tiles2Reference import reprojectSubsettif 
    
    
    
    outputDir='F:\\JPL-PT_MODIS_Andalucia\\Temp_Dasemon\\'
    fileName=outputDir+'temp.tif'
    
    i=int(JulianDay)-1           
    Temp=np.float64(InDataset.variables['pred'][i].data) 
    reprojectDASEMON(fileName, JulianDay,Temp)
    dataOut = Temp*1
    dataOut = dataOut.astype(int)
      
    Temp[Temp==-100000.000]=np.nan
    Dictionary ['SunHours'] =  np.int16(Temp *10) 
    Temp=None
    Temp=np.float64(precip.variables['pr'][i].data)  
    Temp[Temp==-100000.000]=np.nan     
    Dictionary ['Precipitation']=   Temp  
    Temp=None
    Temp=np.float64(TmaxDataset.variables['pred'][i].data) 
    Temp[Temp==-100000.000]=np.nan       
    Dictionary ['DASEMON_Tmax']=np.int16(Temp)
    Temp=None
    Temp=np.float64(TminDataset.variables['pred'][i].data) 
    Temp[Temp==-100000.000]=np.nan 
    Dictionary ['DASEMON_Tmin']=np.int16(Temp)  
    Temp=None  
    lon=InDataset.variables['lon'][:]
    lat=InDataset.variables['lat'][:]

    lonv, latv = np.meshgrid(lon, lat, sparse=False, indexing='xy')
    #Converts the utm coordinates to geographic
    inProj = Proj(init='epsg:32630')
    outProj = Proj(proj='latlong',datum='WGS84')
    Long,Lat = transform(inProj,outProj,lonv,latv)                    
    Dictionary ['Latitude']= Lat
    Dictionary ['Longitude']= Long 
    np.save('N:\\PythonScripts\\4Students\\Latitude.npy',Lat)
    return Dictionary
        