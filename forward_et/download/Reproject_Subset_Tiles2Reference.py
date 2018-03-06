# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 09:42:16 2018

@author: gmgo
"""

def getSceneMetadata(reference_image):
    import gdal
    metadata = {'xsize':0 ,'ysize': 0, 'bands':0, 'gt':[], 'proj':''}
    # Get reprojection parameters from the first dataset
#    scene = '"'+scene+'":Grid:band1'    
    inputImage = gdal.Open(reference_image,gdal.GA_ReadOnly)
    if inputImage:
        metadata['xsize'] = inputImage.RasterXSize
        metadata['ysize'] = inputImage.RasterYSize
        metadata['bands'] = inputImage.RasterCount
        metadata['gt'] = inputImage.GetGeoTransform()
        metadata['proj'] = inputImage.GetProjection()
        inputImage = None
    return metadata   
    
    
def reprojectSubsettif( inputimage, outputDir,referenceMetadata):
    import gdal
    import os
    from subprocess import Popen 
    import subprocess
    C_gdalwarp = 'C:\WinPython\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\Lib\site-packages\osgeo\gdalwarp'
    
    xsize = referenceMetadata['xsize']
    ysize = referenceMetadata['ysize']
    gt = referenceMetadata['gt']
    proj = referenceMetadata['proj']
    
    
    
    inputMetadata=getSceneMetadata(inputimage)
    outputFiles = []
    if inputimage:#modisScene:
        xsize = referenceMetadata['xsize']
        ysize = referenceMetadata['ysize']
        gt = referenceMetadata['gt']
        proj = referenceMetadata['proj']
        input_image = gdal.Open(inputimage, gdal.GA_ReadOnly)
        inputMetadata=getSceneMetadata(inputimage)
        
        gtInput = input_image.GetGeoTransform()
        modisImage = None
        
        pixSizeModis = [gt[1],gt[1]]
        
        UL = [gt[0], gt[3]]
        xsizeModis = xsize
        ysizeModis = ysize
        BR = [UL[0] + xsizeModis*pixSizeModis[0], UL[1] - ysizeModis*pixSizeModis[1]]
        file_name = os.path.basename(inputimage)
        
        outputImagePath=outputDir+file_name
        
        gdalCommand = C_gdalwarp+' -r near -t_srs '+str(proj)+' -te '+str(UL[0])+' '+str(BR[1])+' '+str(BR[0])+' '+str(UL[1])+' -tr '+str(pixSizeModis[0])+' '+str(pixSizeModis[1])+' -overwrite '+inputimage+' '+outputImagePath
        proc=subprocess.Popen(gdalCommand,shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.STDOUT, universal_newlines=False)
        for line in iter(proc.stdout.readline, ""):
            print line
        proc.wait()
        outputFiles.append(outputImagePath)

    return outputFiles  