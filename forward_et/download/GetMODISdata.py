# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 11:46:51 2017

@author: gmgo
"""


# Funcion para borrar todos los archivos intermedios y no crear mucha basura en el disco duro.
def delfilesinfolder(temporalfolder):
    import os
    filelist = [f for f in os.listdir(temporalfolder) if f.endswith(".hdf")]
    for f in filelist:
        remf = temporalfolder+f
        os.remove(remf)
    filelist = [f for f in os.listdir(temporalfolder) if f.endswith(".xml")]
    for f in filelist:
        remf = temporalfolder+f
        os.remove(remf)
    filelist = [f for f in os.listdir(temporalfolder) if f.endswith(".txt")]
    for f in filelist:
        remf = temporalfolder+f
        os.remove(remf)
    filelist = [f for f in os.listdir(temporalfolder) if f.endswith(".tif")]
    for f in filelist:
        remf = temporalfolder+f
        os.remove(remf)
    return
# Function to download the MOD product and make the mosaic


def downloadmosaic(outfolder, temporalfolder, product, tiles, password, user, date, delta):
    import pymodis
    import glob
    import os
    delfilesinfolder(temporalfolder)
    year = str(date.year)
    month = date.month
    if month < 10:
        monthstr = '0'+str(month)
    else:
        monthstr = str(month)
    day = (date.day)
    if day < 10:
        daystr = '0'+str(day)
    else:
        daystr = str(day)
    day = year+'-'+monthstr+'-'+daystr
    if (product == 'MOD09GA.006') or (product == 'MOD11A1.006'):
        path_download = 'MOLT'
        if product == 'MOD09GA.006':
            subset = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Bandas que queremos. 1 yes 0 No
        if product == 'MOD11A1.006':
            subset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
    if (product == 'MYD09GA.006')or (product == 'MYD11A1.006'):
        path_download = 'MOLA'
        if product == 'MYD09GA.006':
            subset = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Bandas que queremos. 1 yes 0 No
        if product == 'MYD11A1.006':
            subset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
    if product == 'MOD13A2.006':
        path_download = 'MOLT'
        if product == 'MOD13A2.006':
            subset = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    if product == 'MCD15A2H.006':
        path_download = 'MOTA'
        if product == 'MCD15A2H.006':
            subset = [1, 1, 1, 1, 1, 1]
    if product == 'MCD43B3.005':
        path_download = 'MOTA'
        if product == 'MCD43B3.005':
            subset = [0, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    if product == 'MCD43B2.005':
        path_download = 'MOTA'
        if product == 'MCD43B2.005':
            subset = [1, 0, 0, 0]
    modisDown = pymodis.downmodis.downModis(product=product, destinationFolder=temporalfolder,
                                            password=password, path=path_download, user=user, tiles=tiles, today=day, delta=delta)
    modisDown.connect()
    modisDown.downloadsAllDay()

    # create the list of files to use
    files = glob.glob(os.path.join(temporalfolder, product[0:7]+'*.hdf'))

    modisParse = pymodis. parsemodis.parseModis(files[0])
    # bounding box of the tiles
    modisParse.retBoundary()
    # Parse multiple data
    modisMultiParse = pymodis.parsemodis.parseModisMulti(files)
    dateout = files[0]
    pos_ini = dateout.index('.A')
    dateout = dateout[pos_ini+2:pos_ini+2+7]
    modisMultiParse.writexml(os.path.join(
        temporalfolder, 'modismultiparse.xml'))
    from pymodis.convertmodis_gdal import createMosaicGDAL
    output_tif = os.path.join(outfolder, product+dateout+'.mosaic.tif')

    mosaic = createMosaicGDAL(files, subset, 'GTiff')
    mosaic.run(output_tif)
    print('Mosaic for '+dateout+' generated!!!')
    return output_tif, day


def DownloadAllData(tiles, date, date_list_8days, Test, Test2, outfolder):
    print date
#    download_09Products(tiles,date,outfolder)
#    download_11Products(tiles,date,outfolder)
#
    Position = Comp8daydate(date)
    while Position != Test:

        #        download_15Products(tiles,date_list_8days[Position],outfolder)
        #        download_43Products(tiles,date_list_8days[Position],outfolder)
        Test = Position

    while Position != Test2:
        #        download_13Products(tiles,date_list_8days[Position],outfolder)
        Test2 = Position
    return Test, Test2


def download_09Products(tiles, date, outfolder):
    import Prepare_input_data
    mrtpath = 'C:\\MRT\\'
    outfolder = outfolder+'MOD_MYD_09\\'
    temporalfolder = outfolder+'Temporal\\'
    product_MOD = 'MOD09GA.006'
    product_MYD = 'MYD09GA.006'
    # tiles='h17v04,h17v05,h18v04,h18v05'
    password = 'Kalimero1980OG'
    user = 'gmendiguren'
    delta = 1
    Prepare_input_data.cleanDirectory(outfolder)
    Prepare_input_data.cleanDirectory(temporalfolder)
    downloadmosaic(outfolder, temporalfolder, product_MOD,
                   tiles, password, user, date, delta)
    downloadmosaic(outfolder, temporalfolder, product_MYD,
                   tiles, password, user, date, delta)
    return


def download_11Products(tiles, date, outfolder):
    import Prepare_input_data
    outfolder = outfolder+'MOD_MYD_11\\'
    temporalfolder = outfolder+'Temporal\\'
    product_MOD = 'MOD11A1.006'
    product_MYD = 'MYD11A1.006'
    # tiles='h17v04,h17v05,h18v04,h18v05'
    password = 'Kalimero1980OG'
    user = 'gmendiguren'
    delta = 1
    Prepare_input_data.cleanDirectory(outfolder)
    Prepare_input_data.cleanDirectory(temporalfolder)
    downloadmosaic(outfolder, temporalfolder, product_MOD,
                   tiles, password, user, date, delta)
    downloadmosaic(outfolder, temporalfolder, product_MYD,
                   tiles, password, user, date, delta)
    return


def download_13Products(tiles, date, outfolder):
    import Prepare_input_data
    outfolder = outfolder+'MOD_13\\'

    temporalfolder = outfolder+'Temporal\\'
    product_MOD = 'MOD13A2.006'

    # tiles='h17v04,h17v05,h18v04,h18v05'
    password = 'Kalimero1980OG'
    user = 'gmendiguren'
    delta = 1
    Prepare_input_data.cleanDirectory(outfolder)
    Prepare_input_data.cleanDirectory(temporalfolder)
    downloadmosaic(outfolder, temporalfolder, product_MOD,
                   tiles, password, user, date, delta)
    return


def download_15Products(tiles, date, outfolder):
    import Prepare_input_data
    outfolder = outfolder+'MCD_15\\'
    temporalfolder = outfolder+'Temporal\\'
    # temporalfolder='C:\\Temporal_Trash\\MCD18\\Temporal\\'
    product_MCD_15 = 'MCD15A2H.006'
    # tiles='h17v04,h17v05,h18v04,h18v05'
    password = 'Kalimero1980OG'
    user = 'gmendiguren'
    delta = 1
    Prepare_input_data.cleanDirectory(outfolder)
    Prepare_input_data.cleanDirectory(temporalfolder)
    downloadmosaic(outfolder, temporalfolder, product_MCD_15,
                   tiles, password, user, date, delta)
#    GetMODISdata.downloadmosaic(outfolder,temporalfolder,product_MYD,tiles,password,user,date,delta)
    return


def download_43Products(tiles, date, outfolder):
    import Prepare_input_data
    outfolder = outfolder+'\\MCD_43\\'
    temporalfolder = outfolder+'Temporal\\'
    # temporalfolder='C:\\Temporal_Trash\\MCD18\\Temporal\\'
    product_MCD_43 = 'MCD43B3.005'
    Quality_product_MCD_43 = 'MCD43B2.005'
    tiles = 'h17v04,h17v05,h18v04,h18v05'
    password = 'Kalimero1980OG'
    user = 'gmendiguren'
    delta = 1
    Prepare_input_data.cleanDirectory(outfolder)
    Prepare_input_data.cleanDirectory(temporalfolder)
    downloadmosaic(outfolder, temporalfolder, product_MCD_43,
                   tiles, password, user, date, delta)
    downloadmosaic(outfolder, temporalfolder,
                   Quality_product_MCD_43, tiles, password, user, date, delta)
#    GetMODISdata.downloadmosaic(outfolder,temporalfolder,product_MYD,tiles,password,user,date,delta)
    return


def Comp8daydate(date):
    from jdcal import gcal2jd

    year = str(date.year)
    month = date.month
    if month < 10:
        monthstr = '0'+str(month)
    else:
        monthstr = str(month)
    day = (date.day)
    if day < 10:
        daystr = '0'+str(day)
    else:
        daystr = str(day)
    day_full = year+'-'+monthstr+'-'+daystr

    Jul = gcal2jd(year, 1, 1)
    JulAct = gcal2jd(year, month, day)
    JulianDay = (JulAct[1]+1-Jul[1])
    JulDay8 = ('001', '009', '017', '025', '033', '041', '049', '057', '065', '073', '081',
               '089', '097', '105', '113', '121', '129', '137', '145', '153', '161', '169',
               '177', '185', '193', '201', '209', '217', '225', '233', '241', '249', '257',
               '265', '273', '281', '289', '297', '305', '313', '321', '329', '337', '345',
               '353', '361')
    JulDayIni = ('001', '009', '017', '025', '033', '041', '049', '057', '065', '073', '081',
                 '089', '097', '105', '113', '121', '129', '137', '145', '153', '161', '169',
                 '177', '185', '193', '201', '209', '217', '225', '233', '241', '249', '257',
                 '265', '273', '281', '289', '297', '305', '313', '321', '329', '337', '345',
                 '353')
    JulDayEnd = ('009', '017', '025', '033', '041', '049', '057', '065', '073', '081', '089',
                 '097', '105', '113', '121', '129', '137', '145', '153', '161', '169', '177',
                 '185', '193', '201', '209', '217', '225', '233', '241', '249', '257', '265',
                 '273', '281', '289', '297', '305', '313', '321', '329', '337', '345', '353',
                 '361')
    for i in range(0, len(JulDay8)):
        if JulianDay == int(JulDay8[i]):
            print i
            return i
        if (JulianDay > int(JulDayIni[i])) and (JulianDay < int(JulDayEnd[i])):
            print i
            return i
        if JulianDay >= 361:
            i = 45
            print i
            return i


def DownloadPJLPT_MODISDATA(DOY, YEAR, TILE, temporalfolder, outfolder):
    delfilesinfolder(temporalfolder)
    print('Starting to gather the MODIS data to run the JPL-PT model for ' +
          str(YEAR)+' '+str(DOY))
    mrtpath = 'C:\\MRT\\'  # directorio del Modis reprojection tool
    # test para comprobar que el MRT está instalado y en su sitio.
    pymodis.convertmodis.checkMRTpath(mrtpath)
    # la contraseña Teneis que registraros en esta web: https://urs.earthdata.nasa.gov/home
    password = 'Kalimero1980OG'
    user = 'gmendiguren'  # el usurario
    delta = 1
    year2process = YEAR
    tiles = 'h17v04,h17v05,h18v04,h18v05'
    base = datetime.datetime(year2process, 01, 01)
    temporalfolder = temporalfolder
    date_list = [base + datetime.timedelta(days=x)
                 for x in range(DOY-1, DOY, 1)]
    return
