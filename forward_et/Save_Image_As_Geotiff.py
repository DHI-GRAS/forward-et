# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 09:52:41 2018

@author: gmgo
"""

import gdal
import numpy as np


def getSceneMetadata(reference_image):
    metadata = {'xsize': 0, 'ysize': 0, 'bands': 0, 'gt': [], 'proj': ''}
    # Get reprojection parameters from the first dataset
#    scene = '"'+scene+'":Grid:band1'
    inputImage = gdal.Open(reference_image, gdal.GA_ReadOnly)
    if inputImage:
        metadata['xsize'] = inputImage.RasterXSize
        metadata['ysize'] = inputImage.RasterYSize
        metadata['bands'] = inputImage.RasterCount
        metadata['gt'] = inputImage.GetGeoTransform()
        metadata['proj'] = inputImage.GetProjection()
        inputImage = None
    return metadata


def Save_array_tiff(array, outputfile, proj, geotransform):

    from osgeo import gdal, osr

    # Start the gdal driver for GeoTIFF
    driver = gdal.GetDriverByName("GTiff")

    # Write array
    shape = array.shape
    print(len(shape))

    if len(shape) > 2:
        ds = driver.Create(
            outputfile, shape[1], shape[0], shape[2], gdal.GDT_Float32)
        ds.SetProjection(proj)
        ds.SetGeoTransform(geotransform)

        for i in range(shape[2]):
            print(i)
            ds.GetRasterBand(i+1).WriteArray(array[:, :, i])

    else:
        ds = driver.Create(outputfile, shape[1], shape[0], 1, gdal.GDT_Float32)
        ds.SetProjection(proj)
        ds.SetGeoTransform(geotransform)
        ds.GetRasterBand(1).WriteArray(array)

    ds = None

    print('Saved ' + outputfile)

    return


def SaveRn_Image(Array):
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    ReferenceIMage = folder+'\\Data\\AndaluciaCoordinates'
    outputfile = folder + '\\OutPutImage\\Rn_Andalucia_All.tif'
    referenceMetadata = getSceneMetadata(ReferenceIMage)
    proj = referenceMetadata.get('proj')
    geotransform = referenceMetadata.get('gt')
    Save_array_tiff(Array, outputfile, proj, geotransform)
    return


def SaveCanopyTransT_Image(Array):
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    ReferenceIMage = folder+'\\Data\\AndaluciaCoordinates'
    outputfile = folder + '\\OutPutImage\\CanopyTrans_Andalucia_All.tif'
    referenceMetadata = getSceneMetadata(ReferenceIMage)
    proj = referenceMetadata.get('proj')
    geotransform = referenceMetadata.get('gt')
    Save_array_tiff(Array, outputfile, proj, geotransform)
    return


def SaveSoilEva_Image(Array):
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    ReferenceIMage = folder+'\\Data\\AndaluciaCoordinates'
    outputfile = folder + '\\OutPutImage\\SoilEva_Andalucia_All.tif'
    referenceMetadata = getSceneMetadata(ReferenceIMage)
    proj = referenceMetadata.get('proj')
    geotransform = referenceMetadata.get('gt')
    Save_array_tiff(Array, outputfile, proj, geotransform)
    return


def SavePotEvapotrans_Image(Array):
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    ReferenceIMage = folder+'\\Data\\AndaluciaCoordinates'
    outputfile = folder + '\\OutPutImage\\PotEvapotrans_Andalucia_All.tif'
    referenceMetadata = getSceneMetadata(ReferenceIMage)
    proj = referenceMetadata.get('proj')
    geotransform = referenceMetadata.get('gt')
    Save_array_tiff(Array, outputfile, proj, geotransform)
    return


def SaveActCanopyTransT_Image(Array):
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    ReferenceIMage = folder+'\\Data\\AndaluciaCoordinates'
    outputfile = folder + '\\OutPutImage\\ActCanopyTrans_Andalucia_All.tif'
    referenceMetadata = getSceneMetadata(ReferenceIMage)
    proj = referenceMetadata.get('proj')
    geotransform = referenceMetadata.get('gt')
    Save_array_tiff(Array, outputfile, proj, geotransform)
    return


def SaveActSoilEva_Image(Array):
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    ReferenceIMage = folder+'\\Data\\AndaluciaCoordinates'
    outputfile = folder + '\\OutPutImage\\ActSoilEva_Andalucia_All.tif'
    referenceMetadata = getSceneMetadata(ReferenceIMage)
    proj = referenceMetadata.get('proj')
    geotransform = referenceMetadata.get('gt')
    Save_array_tiff(Array, outputfile, proj, geotransform)
    return


def SaveActEvapotrans_Image(Array):
    import os
    folder = (os.path.dirname(os.path.realpath(__file__)))
    ReferenceIMage = folder+'\\Data\\AndaluciaCoordinates'
    outputfile = folder + '\\OutPutImage\\ActEvapotrans_Andalucia_All.tif'
    referenceMetadata = getSceneMetadata(ReferenceIMage)
    proj = referenceMetadata.get('proj')
    geotransform = referenceMetadata.get('gt')
    Save_array_tiff(Array, outputfile, proj, geotransform)
    return
