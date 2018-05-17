# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:28:35 2018

@author: gmgo
"""

from datetime import date, timedelta, datetime
import time
from osgeo import gdal
from osgeo.gdalnumeric import *
import os
import subprocess
import re
import gc
import numpy
from gdalconst import *
from ecmwfapi import ECMWFDataServer
from time import sleep
import os
import glob

# MODIS sinusoidal projection with extent of tile h18v03
dstProj = '+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs'
dstTr = '1100.00 -1100.00'
#dstTe = '-0.00000011548 5559752.5983332400 1111950.5196665400 6671703.1179999000'
dstTe = '-80950.0000000000000000 3979450.0000000000000000 1145550.0000000000000000 4896850.0000000000000000'


# Constants
# projection of the ERA interim files
C_srcProj = '+proj=longlat +a=6367470 +b=6367470 +no_defs'
C_gdal_translate = r'C:\WinPython\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\Lib\site-packages\osgeo\gdal_translate'
C_gdalwarp = r'C:\WinPython\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\Lib\site-packages\osgeo\gdalwarp'
C_gdal_calc = r'C:\WinPython\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\Lib\site-packages\osgeo\gdal_calc.py'

C_vrtTempFile = 'tmp.vrt'
C_vrtTempFile2 = 'tmp2.vrt'
MainFolder = 'F:\\PT-JPL_MODIS\\'  # Folder containing all inputs
# forecast air temp, wind dir and dew point temp
meteo = {
    'dataset': "interim",
    'stream': "oper",
              'step': "00/3/6/9/12",
              'number': "all",
              'levtype': "sfc",
              'date': "2012-01-01/to/2012-12-31",
              'time': "00/12",
              'origin': "all",
              'type': "fc",
              'param': "2t",
              'grid': "0.75/0.75",
              'class': "ei",
              'target': MainFolder + 'meteo_xxxx.grib'
}

Precipitation = {
    'dataset': "interim",
    'stream': "oper",
              'step': "00/3/6/9/12",
              'number': "all",
              'levtype': "sfc",
              'date': "2012-01-01/to/2012-12-31",
              'time': "00/12",
              'origin': "all",
              'type': "fc",
              'param': "tp",
              'grid': "0.75/0.75",
              'class': "ei",
              'target': MainFolder + 'meteo_xxxx.grib'
}

# forecast downward and total-clear-sky solar radiation
# https://software.ecmwf.int/wiki/pages/viewpage.action?pageId=56658233
solarRad = {
    'dataset': "interim",
    'stream': "oper",
              'step': "00/3/6/9/12",  # "3/6/9/12",
              'number': "all",
              'levtype': "sfc",
              'date': "2002-01-01/to/2002-12-31",
              'time': "00/12",
              'origin': "all",
              'type': "fc",
              'param': "169.128",  # 'param'   : "169.128/210.128/175.128",#Net Radiation176.128
              'grid': "0.75/0.75",
              'class': "ei",
              'target': MainFolder + 'solarRad_fc_xxxx.grib'
}


retrieveRequests = [Precipitation, solarRad, meteo]


def downloadEraData(request):
    # To run this example, you need an API key
    # available from https://api.ecmwf.int/v1/key/

    server = ECMWFDataServer()

    server.retrieve(request)


def deleteFile(path_to_file, attempts=0, timeout=100, sleep_int=2):
    if attempts < timeout and os.path.exists(path_to_file) and os.path.isfile(path_to_file):
        try:
            os.remove(path_to_file)
            return True
        except BaseException:
            # perform an action
            sleep(sleep_int)
            deleteFile(path_to_file, attempts + 1)
    return False


def multiplyBy100(fileName, bandNum):

    # Open the dataset
    ds = gdal.Open(fileName, GA_ReadOnly)
    band = ds.GetRasterBand(bandNum)
    data = BandReadAsArray(band)

    # The actual calculation
    dataOut = data * 100
    dataOut = dataOut.astype(int)

    # Write the out file
    driver = gdal.GetDriverByName("GTiff")
    outFileName = fileName + '_m100.tif'
    dsOut = driver.Create(outFileName, ds.RasterXSize,
                          ds.RasterYSize, 1, gdal.GDT_Int16)
    CopyDatasetInfo(ds, dsOut)
    dsOut.GetRasterBand(1).WriteArray(dataOut)
    #BandWriteArray(bandOut, dataOut)

    # Close the datasets
    data = None
    band = None
    ds = None
    dataOut = None
    dsOut = None

    return outFileName


def dailySSRD(fileName, bandNum):

    # Open the dataset
    ds = gdal.Open(fileName, GA_ReadOnly)
    band = ds.GetRasterBand(bandNum)
    data1 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 1)
    data2 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 2)
    data3 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 3)
    data4 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 4)
    data5 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 5)
    data6 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 6)
    data7 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 7)
    data8 = BandReadAsArray(band)
    SSRD = (data1 + data2 + data3 + data4 + data5 + data6 + data7 + data8) / 86400
    # The actual calculation
    dataOut = SSRD * 1
    dataOut = dataOut.astype(int)

    # Write the out file
    driver = gdal.GetDriverByName("GTiff")
    outFileName = fileName + '_m1.tif'
    dsOut = driver.Create(outFileName, ds.RasterXSize,
                          ds.RasterYSize, 1, gdal.GDT_Float32)
    CopyDatasetInfo(ds, dsOut)
    dsOut.GetRasterBand(1).WriteArray(dataOut)

    # Close the datasets
    data = None
    band = None
    ds = None
    dataOut = None
    dsOut = None

    return outFileName


def dailyRAin(fileName, bandNum):

    # Open the dataset
    ds = gdal.Open(fileName, GA_ReadOnly)
    band = ds.GetRasterBand(bandNum)
    data1 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 1)
    data2 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 2)
    data3 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 3)
    data4 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 4)
    data5 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 5)
    data6 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 6)
    data7 = BandReadAsArray(band)
    band = ds.GetRasterBand(bandNum + 7)
    data8 = BandReadAsArray(band)
    DPPT = (data1 + data2 + data3 + data4 + data5 + data6 + data7 + data8)
    # The actual calculation
    dataOut = (DPPT * 1000.)
    #dataOut = dataOut.astype(int)

    # Write the out file
    driver = gdal.GetDriverByName("GTiff")
    outFileName = fileName + '_m1.tif'
    dsOut = driver.Create(outFileName, ds.RasterXSize,
                          ds.RasterYSize, 1, gdal.GDT_Float32)
    CopyDatasetInfo(ds, dsOut)
    dsOut.GetRasterBand(1).WriteArray(dataOut)

    # Close the datasets
    data = None
    band = None
    ds = None
    dataOut = None
    dsOut = None

    return outFileName


def airtemptime09(fileName, bandNum):
    import numpy as np
    # Open the dataset
    ds = gdal.Open(fileName, GA_ReadOnly)
    band = ds.GetRasterBand(bandNum + 3)
    data_time09 = BandReadAsArray(band)
    # The actual calculation
    dataOut = data_time09 * 1.0
    # Write the out file
    driver = gdal.GetDriverByName("GTiff")
    outFileName = fileName + '_m1.tif'
    dsOut = driver.Create(outFileName, ds.RasterXSize,
                          ds.RasterYSize, 1, gdal.GDT_Float32)
    CopyDatasetInfo(ds, dsOut)
    dsOut.GetRasterBand(1).WriteArray(dataOut)
    # Close the datasets
    data = None
    band = None
    ds = None
    dataOut = None
    dsOut = None

    return outFileName


def airtemptime12(fileName, bandNum):
    import numpy as np
    # Open the dataset
    ds = gdal.Open(fileName, GA_ReadOnly)
    band = ds.GetRasterBand(bandNum + 4)
    data_time12 = BandReadAsArray(band)
    # The actual calculation
    dataOut = data_time12 * 1.0
    # Write the out file
    driver = gdal.GetDriverByName("GTiff")
    outFileName = fileName + '_m1.tif'
    dsOut = driver.Create(outFileName, ds.RasterXSize,
                          ds.RasterYSize, 1, gdal.GDT_Float32)
    CopyDatasetInfo(ds, dsOut)
    dsOut.GetRasterBand(1).WriteArray(dataOut)
    # Close the datasets
    data = None
    band = None
    ds = None
    dataOut = None
    dsOut = None

    return outFileName


def airtemptime15(fileName, bandNum):
    import numpy as np
    # Open the dataset
    ds = gdal.Open(fileName, GA_ReadOnly)
    band = ds.GetRasterBand(bandNum + 5)
    data_time15 = BandReadAsArray(band)
    # The actual calculation
    dataOut = data_time15 * 1.0
    # Write the out file
    driver = gdal.GetDriverByName("GTiff")
    outFileName = fileName + '_m1.tif'
    dsOut = driver.Create(outFileName, ds.RasterXSize,
                          ds.RasterYSize, 1, gdal.GDT_Float32)
    CopyDatasetInfo(ds, dsOut)
    dsOut.GetRasterBand(1).WriteArray(dataOut)
    # Close the datasets
    data = None
    band = None
    ds = None
    dataOut = None
    dsOut = None

    return outFileName


def gdal2GeoTiff_ECMWF_WGS84(filename):
    print("Translating to GeoTIFF...")
    tiff_filename_base = os.path.split(filename)[0] + os.sep
    tiff_filelist = []
    vrtTempFile = tiff_filename_base + C_vrtTempFile
    print vrtTempFile + 'cacacac'
    d = datetime(1970, 1, 1)

    # Read raster bands from file
    im = gdal.Open(filename, GA_ReadOnly)
    bandsNum = im.RasterCount
    srcProj = C_srcProj

    for i in range(1, bandsNum + 1, 8):
        print(str(i))

        # Get the band properties
        vrtTempFile = tiff_filename_base + str(i) + C_vrtTempFile
        band = im.GetRasterBand(i)
        forecastSeconds = int(re.findall(
            '\d+', band.GetMetadata()['GRIB_FORECAST_SECONDS'])[0])
        UTCtime_delta = int(re.findall(
            '\d+', band.GetMetadata()['GRIB_REF_TIME'])[0]) + forecastSeconds
        productName = band.GetMetadata()['GRIB_ELEMENT']
        band = None
        if not os.path.exists(os.path.join(tiff_filename_base, productName)):
            os.makedirs(os.path.join(tiff_filename_base, productName))
#        tiff_filename = os.path.join(tiff_filename_base, productName, productName+'_'+str((d + timedelta(seconds=UTCtime_delta)).year) + \
#                        str((d + timedelta(seconds=UTCtime_delta)).timetuple().tm_yday).zfill(3) + \
#                        str((d + timedelta(seconds=UTCtime_delta)).hour).zfill(2) + '_ECMWF.tif')
        tiff_filename = os.path.join(tiff_filename_base, productName, productName + '_' + str((d + timedelta(seconds=UTCtime_delta)).year) +
                                     str((d + timedelta(seconds=UTCtime_delta)).timetuple().tm_yday).zfill(3) + '_ECMWF.tif')

        print(str((d + timedelta(seconds=UTCtime_delta)).timetuple().tm_yday).zfill(3))
        print('Reprojecting')
        # Then subset and reproject to one MODIS tile using gdalwarp
        if productName == '2T':

            tiff_filename = os.path.join(tiff_filename_base, productName, productName + '_' + str((d + timedelta(seconds=UTCtime_delta)).year) +
                                         str((d + timedelta(seconds=UTCtime_delta)).timetuple().tm_yday).zfill(3) + '_09_ECMWF.tif')
            # Convert the band to integer to save space and save it as geotiff
            tempFile = airtemptime09(filename, i)

            gdalWarpCmd = C_gdalwarp + ' -wm 2000 -s_srs "' + srcProj + '" -t_srs "' + dstProj + '" -tr ' + dstTr + ' -te ' + dstTe + \
                ' -r near -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -overwrite ' + \
                tempFile + ' ' + tiff_filename
            proc = subprocess.Popen(gdalWarpCmd, shell=True, stdout=subprocess.PIPE,
                                    stdin=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=False)
            tiff_filename = None

            tiff_filename = os.path.join(tiff_filename_base, productName, productName + '_' + str((d + timedelta(seconds=UTCtime_delta)).year) +
                                         str((d + timedelta(seconds=UTCtime_delta)).timetuple().tm_yday).zfill(3) + '_12_ECMWF.tif')
            # Convert the band to integer to save space and save it as geotiff
            tempFile = airtemptime12(filename, i)

            gdalWarpCmd = C_gdalwarp + ' -wm 2000 -s_srs "' + srcProj + '" -t_srs "' + dstProj + '" -tr ' + dstTr + ' -te ' + dstTe + \
                ' -r near -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -overwrite ' + \
                tempFile + ' ' + tiff_filename
            proc = subprocess.Popen(gdalWarpCmd, shell=True, stdout=subprocess.PIPE,
                                    stdin=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=False)
            tiff_filename = None

            tiff_filename = os.path.join(tiff_filename_base, productName, productName + '_' + str((d + timedelta(seconds=UTCtime_delta)).year) +
                                         str((d + timedelta(seconds=UTCtime_delta)).timetuple().tm_yday).zfill(3) + '_15_ECMWF.tif')
            # Convert the band to integer to save space and save it as geotiff
            tempFile = airtemptime15(filename, i)

            gdalWarpCmd = C_gdalwarp + ' -wm 2000 -s_srs "' + srcProj + '" -t_srs "' + dstProj + '" -tr ' + dstTr + ' -te ' + dstTe + \
                ' -r near -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -overwrite ' + \
                tempFile + ' ' + tiff_filename
            proc = subprocess.Popen(gdalWarpCmd, shell=True, stdout=subprocess.PIPE,
                                    stdin=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=False)
            tiff_filename = None

        if productName == 'SSRD':
            # Test to correct
            tempFile = dailySSRD(filename, i)
            gdalWarpCmd = C_gdalwarp + ' -wm 2000 -s_srs "' + srcProj + '" -t_srs "' + dstProj + '" -tr ' + dstTr + ' -te ' + dstTe + \
                ' -r near -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -overwrite ' + \
                tempFile + ' ' + tiff_filename
            proc = subprocess.Popen(gdalWarpCmd, shell=True, stdout=subprocess.PIPE,
                                    stdin=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=False)
        if productName == 'TP':
            # Test to correct
            tempFile = dailyRAin(filename, i)
            gdalWarpCmd = C_gdalwarp + ' -wm 2000 -s_srs "' + srcProj + '" -t_srs "' + dstProj + '" -tr ' + dstTr + ' -te ' + dstTe + \
                ' -r near -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -overwrite ' + \
                tempFile + ' ' + tiff_filename
            proc = subprocess.Popen(gdalWarpCmd, shell=True, stdout=subprocess.PIPE,
                                    stdin=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=False)

        for line in iter(proc.stdout.readline, ""):
            print line
        proc.wait()

        gc.collect()
        deleteFile(tempFile)
        print(str((d + timedelta(seconds=UTCtime_delta)).timetuple().tm_yday).zfill(3))

        # tiff_filelist.append(tiff_filename)

    im = None
    return tiff_filelist


years = [2008]
if __name__ == "__main__":
    for request in retrieveRequests:
        for year in years:
            retrievalPeriods = [str(year) + '-01-01/to/' + str(year) + '-12-31']
            for period in retrievalPeriods:
                request['date'] = period
                request['target'] = request['target'][0:-9] + period[0:4] + '.grib'
                if not os.path.exists(request['target']):
                    downloadEraData(request)
                gdal2GeoTiff_ECMWF_WGS84(request['target'])
                files = glob.glob(MainFolder + '*.grib')
                for f in files:
                    os.remove(f)
