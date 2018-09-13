//TODO: Q&A: the implementation is designed for producing daily AET within a certain year.
//TODO: Q&A: the GEE export output limitation
//TODO: Q&A: the UI enables extent selection (given coordinates or drawing extent), date selection
//TODO: Q&A: if we have time, I'm curious to know a bit more on the ET estimate and the f constrains

//**************************************************
// (P_0) Definition of inputs as global variables
//**************************************************

// Coordinates of points to build extent
var I_lon, I_lat;

// Extent
var I_clip;

//Study period
var I_startDate, I_endDate;


//**************************************************
// (P_1) Main ET module
//**************************************************
var getCollection = function(){
  //define variable
  var startDate = ee.Date(I_startDate);
  // var endDate = ee.Date(I_endDate);
  var endDate = ee.Date(I_endDate).advance(1, 'day');

  var endDateDOY = ee.Date(I_endDate);

  var requiredYear = startDate.get('year').toInt();


  //convert to daily frequency
  var startDoy = ee.Number.parse(startDate.getRelative('day', 'year')).add(1);
  var endDoy = ee.Number.parse(endDateDOY.getRelative('day', 'year')).add(1);
  var days = ee.List.sequence(startDoy,endDoy);

  var nod = endDate.difference(startDate, 'day');

  print ('Number of days/outputs', nod);


  var I_sr = ['state_1km', 'SensorZenith', 'sur_refl_b01', 'sur_refl_b02'];
  var I_te = ['LST_Day_1km', 'QC_Day', 'Day_view_time', 'Day_view_angle',
              'Emis_31', 'Emis_32'];
  var I_laifpar = ['Fpar', 'Lai', 'FparExtra_QC'];
  var I_ndvi = ['NDVI', 'SummaryQA'];
  var I_albedo = ['Albedo_BSA_shortwave', 'Albedo_WSA_shortwave',
                  'BRDF_Albedo_Band_Mandatory_Quality_shortwave'];
  var I_climate = ['AvgSurfT_inst', 'SWdown_f_tavg',
                  'LWdown_f_tavg'];
  //Gather GEE image collections
 var MODIS_SR = ee.ImageCollection("MODIS/006/MOD09GA")
        .filterDate(startDate, endDate).map(getDOY).map(reprojMODIS).select(I_sr).map(maskSR)
        .filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy))),
    MODIS_TE_TERRA = ee.ImageCollection("MODIS/006/MOD11A1")
        .filterDate(startDate, endDate).map(getDOY).map(reprojMODIS).select(I_te).map(renameTEterra).map(maskTEterra)
        .filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy))),
    MODIS_TE_AQUA = ee.ImageCollection("MODIS/006/MYD11A1")
        .filterDate(startDate, endDate).map(getDOY).map(reprojMODIS).select(I_te).map(renameTEaqua).map(maskTEaqua)
        .filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy))),
    MODIS_NDVI = ee.ImageCollection("MODIS/006/MOD13A1")
        .filterDate(startDate, endDate).map(getDOY).map(reprojMODIS).select(I_ndvi).map(maskNDVI)
        .filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy))),
    MODIS_LAIFPAR = ee.ImageCollection("MODIS/006/MCD15A3H")
        .filterDate(startDate, endDate).map(getDOY).map(reprojMODIS).select(I_laifpar).map(maskLAIFPAR)
        .filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy))),
    MODIS_BRDFA = ee.ImageCollection("MODIS/006/MCD43A3")
        .filterDate(startDate, endDate).map(getDOY).map(reprojMODIS).select(I_albedo).map(maskALBEDO)
        .filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy))),
    GLDAS = ee.ImageCollection("NASA/GLDAS/V021/NOAH/G025/T3H").filterDate(startDate, endDate).select(I_climate);

  //3-hour to daily
  var GLDAS_DAILY = ee.ImageCollection.fromImages(days.map(function(d) {
    var daily_mean = GLDAS.filter(ee.Filter.calendarRange({
      start: d,
      field: 'day_of_year'
    })).mean();
    var daily_Tmax = GLDAS.select('AvgSurfT_inst').filter(ee.Filter.calendarRange({
      start: d,
      field: 'day_of_year'
    })).max().rename('TMAX');
    var daily_Tmin = GLDAS.select('AvgSurfT_inst').filter(ee.Filter.calendarRange({
      start: d,
      field: 'day_of_year'
    })).min().rename('TMIN');
    return daily_mean.addBands(daily_Tmax).addBands(daily_Tmin).set('day', d);
  })).filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy)));
  print ('GLDAS', GLDAS_DAILY);
  //NDVI: 16-day to daily
  var NDVI_DAILY = make_DailyIC(16, MODIS_NDVI, 'NDVI', days);
  //LAI: 4-day to daily
  var LAI_DAILY = make_DailyIC(4, MODIS_LAIFPAR, 'Lai', days);
  //FPAR: 4-day to daily
  var FPAR_DAILY = make_DailyIC(4, MODIS_LAIFPAR, 'Fpar', days);

  // Define the join and filter
  var Join = ee.Join.inner();
  var FilterOnStartTime = ee.Filter.equals({
                                          'leftField': 'day',
                                          'rightField': 'day'
                                          });

  var TE_Joined = ee.ImageCollection(Join.apply(MODIS_TE_TERRA, MODIS_TE_AQUA, FilterOnStartTime)).map(merge_bands).map(maskTEcomposite);

  var NG_Joined = ee.ImageCollection(Join.apply(NDVI_DAILY, GLDAS_DAILY, FilterOnStartTime)).map(merge_bands);

  var LF_Joined = ee.ImageCollection(Join.apply(LAI_DAILY, FPAR_DAILY, FilterOnStartTime)).map(merge_bands);

  // Join the collections, passing entries through the filter
  var AL_Joined = ee.ImageCollection(Join.apply(MODIS_BRDFA, LF_Joined, FilterOnStartTime)).map(merge_bands);

  var ALTE_Joined = ee.ImageCollection(Join.apply(AL_Joined, TE_Joined, FilterOnStartTime)).map(merge_bands);
  var FinalDataset = ee.ImageCollection(Join.apply(NG_Joined, ALTE_Joined, FilterOnStartTime)).map(merge_bands);

  print ('Dictionary for data', FinalDataset)

  var PotEVTR = FinalDataset.map(get_VIs).map(calc_albedo).select('Albedo','NDVI','LAI', 'FPAR', 'AvgSurfT_inst', 'TMAX', 'TMIN', 'SWdown_f_tavg',
              'LWdown_f_tavg', 'EMIS_comp', 'LST_comp', 'DVT_comp', 'MASK_comp')
              .map(get_SWR_inputs)
              .map(calc_ShortWaveRadiation)
              .map(calc_LongWaveRadiation)
              .map(get_meantemp)
              .map(split_cannopy)
              .map(PET)
              .map(calc_Constrains);
  print ('PET related', PotEVTR)

  //soil moisture constraints
  var SMCcollection = calSOILMoistureConstrain_f_Zhang(requiredYear)
                      .filter(ee.Filter.and(ee.Filter.gte('day', startDoy), ee.Filter.lte('day', endDoy)));

  print ('SMC related', SMCcollection);
  return ee.ImageCollection(PotEVTR);
};

//**************************************************
// (P_2) Run the module with UI inputs and export
//**************************************************
var buildAndExportComposite = function() {
  // clear exceptions from previous run
  panelException.clear();

  //-------------------------------
  // read user Inputs from Graphical User Interface

  //Area of interest (rectangle)
  I_lon = parseFloat(lonBox.getValue());   // parseFloat converts string to number
  I_lat = parseFloat(latBox.getValue());

  var I_clip = ee.Geometry.Point(I_lon, I_lat)
                         .buffer(50000)
                         .bounds();

  I_startDate = I_startDateTexbox.getValue();
  I_endDate = I_endDateTexbox.getValue();

  var collection = getCollection();
  //define a test area
  var composite = collection.select('PET').mean().clip(I_clip);
  addCompositeToMapAndExport(composite);

};
var addCompositeToMapAndExport = function(composite){
  mapnew.layers().set(0, ui.Map.Layer(composite, {min: -5000, max: 0}, 'Daliy Mean PET'));

  // Simple export to Google Drive
  Export.image.toDrive({
    image: composite,
    description: 'FORWARD-ET_tester',
    folder: 'FORWARD-ET',
    scale: 1000,
    region: I_clip,
    maxPixels: 2e12,
    crs:'EPSG:4326'
  });

};

var buildAndPlotComposite = function() {
  // clear exceptions from previous run
  panelException.clear();

  //-------------------------------
  // read user Inputs from Graphical User Interface

  //Area of interest (rectangle)
  I_lon = parseFloat(lonBox.getValue());   // parseFloat converts string to number
  I_lat = parseFloat(latBox.getValue());

  I_startDate = I_startDateTexbox.getValue();
  I_endDate = I_endDateTexbox.getValue();

  var collection = getCollection();

  //plot ET collection
  // add a red pixel to the map where the user clicked or defined a coordinate
  var point = ee.Geometry.Point(I_lon, I_lat);
  var pixel = point.buffer(100).bounds();
  mapnew.layers().set(1, ui.Map.Layer(pixel, {color: 'FF0000'}, 'Plot pixel'));
  var chart = ui.Chart.image.series(collection.select('PET'), pixel, ee.Reducer.mean(), 100, 'day')
    .setOptions({
      title: 'FORWARD-PET',
      lineWidth: 1,
      pointSize: 3,
  });
  plotPanel.add(chart);
};

//**************************************************
// (P_3) Functions
//**************************************************
//produce daily frequency inputs by linear interpolation
//TODO: check the profile
function make_DailyIC(freq, collection, band, days){
  var IC_daily = ee.ImageCollection.fromImages(days.map(function(d) {
    freq = ee.Number(freq);
    var doy = ee.Number(d);
    //get the reminder btw doy and 16-day frequency
    var c = ee.Number(d).subtract(1).mod(freq);
    //get the divisible part btw doy and 16-day frequency
    var e = ee.Number(d).subtract(1).divide(freq).floor();
    //select MODIS image
    var MODIS_sdoy = e.multiply(freq).add(1);
    var MODIS_edoy = MODIS_sdoy.add(freq);
    // doy of the last MODIS image
    var MODIS_fdoy = ee.Number(ee.Image(collection.toList(collection.size()).get(-1)).get('day'));
    var simg = ee.Image(0);
    simg = collection.filter(ee.Filter.eq('day', MODIS_sdoy)).select(band).mean().toDouble();
    var eimg = ee.Image(0);
    eimg = collection.filter(ee.Filter.eq('day', MODIS_edoy)).select(band).mean().toDouble();
    // linear interpolate btw two adjacent images
    var condi = doy.lt(MODIS_fdoy);
    var interp_img = ee.Image(
        ee.Algorithms.If(
          condi,
          eimg.subtract(simg).divide(freq).multiply(c).add(simg),
          simg
        )
      );
    return interp_img.set('day', d);
  }));

  return IC_daily;
}
//create 'day' attribute in property: the calculation amongst different collections is control by 'day' attribute
function getDOY(image) {
  var I_doy = ee.Number.parse(image.date().getRelative('day', 'year')).add(1);
  return image.set('day', I_doy)
}
//add time band in collection as DOY
function getTIME(image) {
  var Doy = ee.Image.constant(ee.Number(image.get('day'))).long().rename('time');
  return image.addBands(Doy);
}

//reproject image to WGS84 projection
function reprojMODIS(image) {
  return image.reproject('EPSG:4326', null, 500).set('system:time_start',
  image.get('system:time_start')).set('day',
  image.get('day'));
}

//This function merges bands from multiply collections
function merge_bands(element) {
  return ee.Image.cat(element.get('primary'), element.get('secondary'));
}

// helper function to extract the QA bits
function getQABits(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    // Return a single band image of the extracted QA bits, giving the band
    // a new name.
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
}

//function for surface reflectance data masking
//TODO: check QC bits
function maskSR(image) {
  // Select the QA band.
  var QA = image.select('state_1km');
  // Get the internal_cloud_algorithm_flag bit.
  var cloud = getQABits(QA, 0, 1, 'cloud');
  var cloudshadow = getQABits(QA, 2, 2, 'cloudshadow');
  // Return an image masking out cloudy areas.
  return image.updateMask(cloud.eq(0).and(cloudshadow.eq(0)));
}

//function for LAI/FPAR data masking
//TODO: check QC bits
function maskLAIFPAR(image) {
  // Select the QA band.
  var QA = image.select('FparExtra_QC');
  var internalQuality = getQABits(QA, 4, 6, 'internal_quality_flag')
  return image.updateMask(internalQuality.eq(0));
  // // Get the internal_cloud_algorithm_flag bit.
  // var snowice = getQABits(QA, 2, 2, 'snow');
  // var aerosol = getQABits(QA, 3, 3, 'aerosol');
  // var cirrus = getQABits(QA, 4, 4, 'cirrus');
  // var cloud = getQABits(QA, 5, 5, 'cloud');
  // var cloudshadow = getQABits(QA, 6, 6, 'cloudshadow');
  // // Return an image masking out cloudy areas.
  // return image.updateMask(cloud.eq(0).and(cloudshadow.eq(0)).and(cirrus.eq(0)));
}

//function for NDVI data masking
//TODO: check QC bits
function maskNDVI(image) {
  // Select the QA band.
  var QA = image.select('SummaryQA');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 0, 1, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.lte(1));
}

//function for Albedo data masking
//TODO: check QC bits
function maskALBEDO(image) {
  // Select the QA band.
  var QA = image.select('BRDF_Albedo_Band_Mandatory_Quality_shortwave');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 0, 0, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.eq(0));
}

//function for Temperature and Emissivity data masking
//Note: Temperature and Emissivity mask function required to create a composite from Terra and Aqua products
//functions include rename the bands and mask for each production and make the composite with rules

//rename terra bands
function renameTEterra(image) {
  var bn = ['LST_Day_1km', 'QC_Day', 'Day_view_time', 'Day_view_angle',
            'Emis_31', 'Emis_32'];
  var vertLabels = [];
  for (var i = 0; i <= 5; i++) {
      var pattern = bn[i] + '_terra';
      vertLabels.push(pattern);
  }
  var renamed = image.select(
    bn, // old names
    vertLabels // new names
  );
  return renamed;
}

//rename aqua bands
function renameTEaqua(image) {
  var bn = ['LST_Day_1km', 'QC_Day', 'Day_view_time', 'Day_view_angle',
            'Emis_31', 'Emis_32'];
  var vertLabels = [];
  for (var i = 0; i <= 5; i++) {
      var pattern = bn[i] + '_aqua';
      vertLabels.push(pattern);
  }
  var renamed = image.select(
    bn, // old names
    vertLabels // new names
  );
  return renamed;
}

//mask terra bands
//TODO: check QC bits
function maskTEterra(image) {
  // Select the QA band.
  var QA = image.select('QC_Day_terra');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 6, 7, 'internal_quality_flag_terra');
  // Return an image masking out cloudy areas.
  return image.addBands(internalQuality.neq(3).remap([1], [1], 0).rename('internal_quality_flag_terra'));
}

//mask aqua bands
//TODO: check QC bits
function maskTEaqua(image) {
  // Select the QA band.
  var QA = image.select('QC_Day_aqua');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 6, 7, 'internal_quality_flag_aqua');
  // Return an image masking out cloudy areas.
  return image.addBands(internalQuality.neq(3).remap([1], [2], 0).rename('internal_quality_flag_aqua'));
}

//Create composite from Terra and Aqua
//TODO: check output
function maskTEcomposite(image) {
  var LST_mask = image.select('internal_quality_flag_terra').add(image.select('internal_quality_flag_aqua'))
                  .rename('MASK_comp');
  var LST_composite = ee.Image(LST_mask.eq(1).multiply(image.select('LST_Day_1km_terra')))
                      .add(ee.Image(LST_mask.eq(2).multiply(image.select('LST_Day_1km_aqua'))))
                      .add(ee.Image(LST_mask.eq(3)).multiply(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).lt(0).multiply(image.select('LST_Day_1km_aqua')))
                      .add(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).gte(0).multiply(image.select('LST_Day_1km_terra'))))))
                      .multiply(0.02).rename('LST_comp');
  var DVT_composite = ee.Image(LST_mask.eq(1).multiply(image.select('Day_view_time_terra')))
                      .add(ee.Image(LST_mask.eq(2).multiply(image.select('Day_view_time_aqua'))))
                      .add(ee.Image(LST_mask.eq(3)).multiply(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).lt(0).multiply(image.select('Day_view_time_aqua')))
                      .add(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).gte(0).multiply(image.select('Day_view_time_terra'))))))
                      .multiply(0.1).rename('DVT_comp');
  var EMIS31_composite = ee.Image(LST_mask.eq(1).multiply(image.select('Emis_31_terra')))
                      .add(ee.Image(LST_mask.eq(2).multiply(image.select('Emis_31_aqua'))))
                      .add(ee.Image(LST_mask.eq(3)).multiply(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).lt(0).multiply(image.select('Emis_31_aqua')))
                      .add(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).gte(0).multiply(image.select('Emis_31_terra'))))))
                      .multiply(0.002).rename('EMIS31_comp');
  var EMIS32_composite = ee.Image(LST_mask.eq(1).multiply(image.select('Emis_32_terra')))
                      .add(ee.Image(LST_mask.eq(2).multiply(image.select('Emis_32_aqua'))))
                      .add(ee.Image(LST_mask.eq(3)).multiply(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).lt(0).multiply(image.select('Emis_32_aqua')))
                      .add(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_angle_terra'))).gte(0).multiply(image.select('Emis_32_terra'))))))
                      .multiply(0.002).rename('EMIS32_comp');
  var EMIS_composite = EMIS31_composite.add(EMIS32_composite).divide(2).multiply(10).rename('EMIS_comp');
  return image.addBands(LST_mask).addBands(LST_composite)
              .addBands(DVT_composite.updateMask(DVT_composite.neq(0).and(DVT_composite.lte(24)))).addBands(EMIS_composite);
}

//fn.GetShortWaveRadiationInputs
function get_SWR_inputs(image) {
  // TODO: sunhours and latitude
  // """Units:
  //     Sunhours=hours*10
  //     Latitude=Degress
  //     Doy=Float64
  // """
  var Sunhours = ee.Image(10).divide(10).rename('Sunhours');
  var Latitude = ee.Image(1).rename('Latitude');
  var Doy = ee.Image(ee.Number(image.get('day'))).rename('Doy');
  return image.addBands(Sunhours.updateMask(Sunhours.neq(0))).addBands(Latitude).addBands(Doy);
}




//rescale the vegetation indexes
function get_VIs(image) {
  var LAI = image.select('Lai').rename('LAI').divide(10);
  var FPAR = image.select('Fpar').rename('FPAR').divide(100);
  var NDVI = image.select('NDVI').divide(10000);
  return image.addBands(LAI).addBands(FPAR).addBands(NDVI);
}


//fn.GetMeanTemperature
//TODO: GLDAS temperature unit
function get_meantemp(image) {
  var AirMax = image.select('TMAX');
  var AirMin = image.select('TMIN');
  var Tmean = (AirMax.add(AirMin)).divide(20).rename('Tmean');
  return image.addBands(Tmean);
}

// Prepare_input_data.py-FillAlbedoInputs
function calc_albedo(image) {
  var BWSA  = image.expression(
      ' (0.8 * BSA * 0.001 + 0.2 * WSA * 0.001) * 100', {
        'BSA': image.select('Albedo_BSA_shortwave'),
        'WSA': image.select('Albedo_WSA_shortwave')
  }).rename('Albedo');
  return image.addBands(BWSA).updateMask(BWSA.lte(100));
}

// rf.CalShortWaveIncomingandNetRadiation
function calc_ShortWaveRadiation(image) {
  var Albedo = image.select('Albedo');
  var Sunhours = image.select('Sunhours');
  var Latitude = image.select('Latitude');
  var Doy = image.select('Doy');

  var fimage = image.select(['Albedo', 'Sunhours', 'Latitude', 'Doy'])

  var AS = ee.Image(0.25);
  var BS = ee.Image(0.50);
  //latRad = (Latitude * np.pi) / 180
  var latRad = ee.Image(Latitude.multiply(Math.PI).divide(180)).rename('latRad');
  var fimage = fimage.addBands(latRad);
  // Declination in radians: delta = 0.409 * np.sin(((Doy * 2 * np.pi) / 365) - 1.39)
  var delta = fimage.expression(
      '0.409 * sin(((Doy * 2 * pi) / 365) - 1.39)', {
        'Doy': fimage.select('Doy'),
        'pi': Math.PI
  }).rename('delta');
  var fimage = fimage.addBands(delta);
  // Eccentricity factor: dr = 1 + 0.033 * np.cos(Doy * ((2 * np.pi) / 365))
  var dr = fimage.expression(
      '1 + 0.033 * cos(Doy * ((2 * pi) / 365)) ', {
        'Doy': fimage.select('Doy'),
        'pi': Math.PI
  }).rename('dr');
  var fimage = fimage.addBands(dr);
  // Sunset hour angle radians ws: = np.arccos(-np.tan(latRad) * np.tan(delta))
  var ws = fimage.expression(
      'acos(-tan(latRad) * tan(delta))', {
        'latRad': fimage.select('latRad'),
        'delta': fimage.select('delta'),
        'pi': Math.PI
  }).rename('ws');
  var fimage = fimage.addBands(ws);
  // Day length in hours: N = (ws * 24.) / np.pi
  var N = ws.multiply(24).divide(Math.PI).rename('N');
  var fimage = fimage.addBands(N);
  // Mjm2day-1: Ra = (((24 * 60) / np.pi) * 0.0820) *
  // dr * ((ws * np.sin(latRad) * np.sin(delta)) +
  // (np.cos(latRad) * np.cos(delta) * np.sin(ws)))
  var Ra = fimage.expression(
      '(((24 * 60) / pi) * 0.0820) * dr * ((ws * sin(latRad) * sin(delta)) + (cos(latRad) * cos(delta) * sin(ws)))', {
        'latRad': fimage.select('latRad'),
        'delta': fimage.select('delta'),
        'dr': fimage.select('dr'),
        'ws': fimage.select('ws'),
        'pi': Math.PI
  }).rename('Ra');
  var fimage = fimage.addBands(Ra);
  // Shortwave Incoming radiation: Rs = Ra * (AS + (BS * SunHours / N))
  var Rs = Ra.multiply(AS.add(BS.multiply(Sunhours).divide(N))).rename('Rs');
  //RSnet = Rs * (1 - Albedo)
  var RSnet = Rs.multiply(Albedo.multiply(-1).add(1)).multiply(11.574).rename('RSnet');
  //Sunset = Sunrise + N = (((ws) * (180 / np.pi)) / 15) + 12
  var Sunrise = ws.multiply(180 / Math.PI).divide(15).add(12).subtract(N).rename('Sunrise');

  return image.addBands(N).addBands(Rs).addBands(RSnet).addBands(Sunrise);
}

//rf.CalLongWave_Incoming_Parton_Logan
function calc_LongWaveRadiation(image) {
  var fimage = image.select(['LST_comp', 'DVT_comp', 'EMIS_comp',
                              'N', 'TMAX', 'TMIN', 'Sunrise', 'RSnet'])
  //c refers to lag of the minimum temperature from the time of sunrise,(parameter
  // following table 1, p.210)
  var c1 = ee.Image(-0.17);
  //d refers to a in Parton and Logan, (parameter following table 1, p.210)
  var d1 = ee.Image(1.86).rename('d1');
  var fimage = fimage.addBands(d1);
  //BB = 12 - (N / 2) + c1
  var BB = ee.Image(12).subtract(fimage.select('N').divide(2)).add(c1);
  //BBD, number of hours after the minimum temperature occcurs until sunset
  var BBD = fimage.select('DVT_comp').subtract(BB).rename('BBD');
  var fimage = fimage.addBands(BBD);
  //Tair_MODIS_passtime = (AirMax - AirMin) * np.sin((np.pi * BBD) / (N + 2 * d1)) + AirMin
  var Tair_MODIS_passtime = fimage.expression(
        '(AirMax - AirMin) * sin((pi * BBD) / (N + 2 * d1)) + AirMin ', {
        'AirMax': fimage.select('TMAX'),
        'AirMin': fimage.select('TMIN'),
        'BBD': fimage.select('BBD'),
        'N': fimage.select('N'),
        'd1': fimage.select('d1'),
        'pi': Math.PI
  }).rename('Tair_MODIS_passtime');
  var fimage = fimage.addBands(Tair_MODIS_passtime);
  //Calculating incoming longwave radiation
  var Stef = ee.Image(5.67e-008).rename('Stef');
  var C = ee.Image(0.261).rename('C');
  var d = ee.Image(7.77e-004).rename('d');
  var fimage = fimage.addBands(Stef).addBands(C).addBands(d);
  //EmissivityAir = 1 - C * np.exp(-d * (Tair_MODIS_passtime - 273.15)**2)
  var EmissivityAir = fimage.expression(
        '1 - C * exp( -d * (Tair_MODIS_passtime - 273.15)**2)', {
        'C': fimage.select('C'),
        'd': fimage.select('d'),
        'Tair_MODIS_passtime': fimage.select('Tair_MODIS_passtime')
  }).rename('EmissivityAir');
  var fimage = fimage.addBands(EmissivityAir);
  //Rlw_in = Stef * (EmissivityAir) * (Tair_MODIS_passtime**4)
  var Rlw_in = fimage.expression(
        'Stef * (EmissivityAir) * (Tair_MODIS_passtime**4)', {
        'Stef': fimage.select('Stef'),
        'EmissivityAir': fimage.select('EmissivityAir'),
        'Tair_MODIS_passtime': fimage.select('Tair_MODIS_passtime')
  }).rename('Rlw_in');

  //CalLongWave_Outgoing
  //Rlw_out = Stef * (Emissivity) * (LST**4)
  var Rlw_out = fimage.expression(
        'Stef * (Emissivity) * (LST**4)', {
        'Stef': fimage.select('Stef'),
        'Emissivity': fimage.select('EMIS_comp'),
        'LST': fimage.select('LST_comp')
  }).rename('Rlw_out');
  //PTJPL_EVAPOTRANSPIRATION.py line 51-52
  var RLnet = Rlw_in.subtract(Rlw_out);
  var t = fimage.select('DVT_comp').subtract(fimage.select('Sunrise')).rename('t');
  var fimage = fimage.addBands(t);
  //fn.calcJparameter
  //J = (2 * N) / (np.pi * np.sin((np.pi * t) / N) * 24)
  var J = fimage.expression(
        '(2 * N) / (pi * sin((pi * t) / N) * 24)', {
        'N': fimage.select('N'),
        't': fimage.select('t'),
        'pi': Math.PI
  });
  //PTJPL_EVAPOTRANSPIRATION.py line 54-55
  var RSnetInstant = fimage.select('RSnet').divide(J);
  var RnDaily = (RSnetInstant.add(RLnet)).multiply(J).rename('RnDaily');
  return image.addBands(RnDaily).addBands(Tair_MODIS_passtime);
}
//TODO: clarify the GLDAS input
function calc_NetRadiation(image) {
  var fimage = image.select(['LST_comp', 'DVT_comp', 'EMIS_comp',
                              'LWdown_f_tavg', 'Swnet_tavg'])
  //CalLongWave_Outgoing
  //Rlw_out = Stef * (Emissivity) * (LST**4)
  var Stef = ee.Image(5.67e-008).rename('Stef');
  var Rlw_out = fimage.expression(
        'Stef * (Emissivity) * (LST**4)', {
        'Stef': fimage.select('Stef'),
        'Emissivity': fimage.select('EMIS_comp'),
        'LST': fimage.select('LST_comp')
  }).rename('Rlw_out');
  var RLnet = fimage.select('LWdown_f_tavg').subtract(Rlw_out);
  //PTJPL_EVAPOTRANSPIRATION.py line 54-55
  var RSnetInstant = fimage.select('Swnet_tavg').divide(J);
  var RnDaily = (RSnetInstant.add(RLnet)).multiply(J).rename('RnDaily');
  return image.addBands(RnDaily);
}


//rf.SplitRn_SOIL_CANOPY
function split_cannopy(image) {
  var kRn = 0.6;
  var LAI = image.select('LAI');
  var Rn = image.select('RnDaily');
  var Rn_SOIL = Rn.multiply((ee.Image(kRn).add(LAI).multiply(-1)).exp()).rename('Rn_SOIL');
  var Rn_Canopy = Rn.subtract(Rn_SOIL).rename('Rn_Canopy');

  return image.addBands(Rn_SOIL).addBands(Rn_Canopy);
}

//This function gives PET following PTJPL_EVAPOTRANSPIRATION.py line 60-66
function PET(image) {
  //TODO: Tmean is GLDAS daily mean temperature
  var Tmean = image.select('AvgSurfT_inst');
  var RnCanopy = image.select('Rn_Canopy');
  var Rn_SOIL = image.select('Rn_SOIL');
  //TODO: DEM check
  var Altitude = ee.Image('CGIAR/SRTM90_V4');

  //PotET.calPsycometricConstant
  //P = 101.3 * (((293 - 0.0065 * ASL_Altitude_m) / 293)**5.26)
  //unit is (kPa)
  var P = image.expression(
        '101.3 * (((293 - 0.0065 * ASL_Altitude_m) / 293)**5.26)', {
        'ASL_Altitude_m': Altitude
  });
  // Latent Heat of Vaporization (lambda)
  //Lambda = 2454 - 2.4 * (Tmean - 20)
  var Lambda = ee.Image(2454).subtract(ee.Image(2.4).multiply(Tmean.subtract(20)));
  // Specific heat of moist air (cp)
  // Unit (kJ kg-1 ºC-1)
  var cp = 1013;
  //Ratio molecular weight of water vapour/dray air
  var epsi = 0.622;
  //Psychrometric constant (psi)
  //Unit (Pa)
  var psi = (P.multiply(cp)).divide(Lambda.multiply(epsi));

  //PotET.calSatVapPres
  //es = 10 * (0.061121 * 2.718281828**(17.502 * Tmean / (240.97 + Tmean))
  //               * (1.0007 + (3.46 * 10**(-8) * 100)))
  var es = image.expression(
        '10 * (0.061121 * 2.718281828**(17.502 * Tmean / (240.97 + Tmean)) * (1.0007 + (3.46 * 10**(-8) * 100)))', {
        'Tmean': Tmean
  });

  //PotET.calSlopVapPresCurv
  //Slope vapour Pressure curve (s)
  //s = (Lambda * 18 * 1000 * es) / (8.3144 * (Tmean + 273)**2)
  var s = image.expression(
        '(Lambda * 18 * 1000 * es) / (8.3144 * (Tmean + 273)**2)', {
        'Lambda': Lambda,
        'es': es,
        'Tmean': Tmean
  });

  //PotET.calPotET_CANOPY
  //Results canopy potential evapotranspiration Etpc % Following Pristley Taylor (1972)(eq 14)
  var alphaPT = ee.Image(1.26);
  //Etpc = alphaPT * (s / (s + psi)) * RnCanopy
  //Unit (Wm-2)
  //PotET_CANOPY
  var Etpc = image.expression(
        'alphaPT * (s / (s + psi)) * RnCanopy', {
        'alphaPT': alphaPT,
        's': s,
        'psi': psi,
        'RnCanopy': RnCanopy
  }).rename('PotET_CANOPY');

  //PotET.calPotET_SOIL
  var G = 0;
  //Etps = alphaPT * (s / (s + psi)) * (Rn_SOIL - G)
  //PotET_SOIL
  var Etps = ee.Image(alphaPT).multiply(s.divide(s.add(psi))).multiply(Rn_SOIL.subtract(G)).rename('PotET_SOIL');

  //PET
  var PotEVTR = Etpc.add(Etps).rename('PET');

  return image.addBands(PotEVTR).addBands(Etpc).addBands(Etps);
}
// calculate the constrains for temperature, vegetation, plant moisture
function calc_Constrains(image) {
  //temperature constraints: ft
  var Topt = ee.Image(25);
  var ft = image.expression(
        '1.1814 * (1 + exp(0.3 * (Topt - 10 - Tmean)))**-1', {
        'Topt': Topt,
        'Tmean': image.select('AvgSurfT_inst')
  }).rename('ft');

  //vegetation constraints: fg
  var fIPAR = ee.Image(1).multiply(image.select('NDVI')).add(ee.Image(-0.05));
  var fg = image.select('FPAR').divide(fIPAR).rename('fg');

  //plant moisture constrains: fm
  //TODO: fPARmax CHANGE LATER
  var fPARmax = ee.Image(1);
  var fm = image.select('FPAR').divide(fPARmax).rename('fm');

  return image.addBands(ft).addBands(fg).addBands(fm);
}
//get images indicating whether the daily rainfall is effective
function effRainyday(image) {
  var Rainyday = image.gte(0.5).toInt().rename('effectRain');
  return image.addBands(Rainyday);

}
//get a collection of images indicating the count number of cumulative rainfall
//TODO: check correctness
function drylength(current_one, previous_all){
  var previous = ee.Image(ee.List(previous_all).get(-1));  // cast to establish the client-side type
  var sum = current_one.add(previous).multiply(current_one).set('day', current_one.get('day'));

  return ee.List(previous_all).add(sum);
}


//get a collection of images indicating the daily fzhang and fdrying constrains
//TODO: GO through the Python code together
function calSOILMoistureConstrain_f_Zhang(requiredYear){
  var sYear = requiredYear;
  var sDay = '-01-01';
  var eDay = '-12-31';
  var sDate = ee.Date(ee.String(ee.Number(sYear).toInt()).cat(sDay));
  var eDate = ee.Date(ee.String(ee.Number(sYear).toInt()).cat(eDay));

  var sDoy = ee.Number.parse(sDate.getRelative('day', 'year')).add(1);
  var eDoy = ee.Number.parse(eDate.getRelative('day', 'year')).add(1);
  print (eDoy, 'endDOY');

  //precipitation
  var CHIRPS_DAILY = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate(sDate, eDate).map(getDOY);
  //rainy day
  var CHIRPS_EFFECT = CHIRPS_DAILY.map(effRainyday).select('effectRain');
  var time0 = CHIRPS_EFFECT.first().get('day');
  var RDcollection = ee.ImageCollection(ee.List(CHIRPS_EFFECT.iterate(drylength, ee.List([
    CHIRPS_EFFECT.select('effectRain').max().set('day', time0)
      .remap([0,1], [0,0], 0, "effectRain")
      .rename('effectRain').cast({'effectRain': 'long'})
  ])))).map(getTIME);
  //soil equilibrium evaporation rate
  //TODO: change later
  var EQcollection = CHIRPS_DAILY;

  var days = ee.List.sequence(sDoy,eDoy);
  var byDay = ee.ImageCollection.fromImages(
        days.map(function (d) {
          var sumPrec = ee.Image(0);
          var sumEqs = ee.Image(0);
          var maxRainyday = RDcollection.filter(ee.Filter.lte('day', d)).max();
          var Delta_t = ee.Image(d-1).subtract(maxRainyday);
          var alpha_f = ee.Image(0.5);
          //****Get the position of the maximum value***//
          // turn image collection into an array
          var array = RDcollection.filter(ee.Filter.lte('day', d)).toArray();

          // sort array by the first band, keeping other bands
          var axes = { image:0, band:1 };
          var sort = array.arraySlice(axes.band, 0, 1);  // select bands from index 0 (inclusive) to 1 (exclusive)
          var sorted = array.arraySort(sort);

          // take the first image only
          var length = sorted.arrayLength(axes.image);
          var values = sorted.arraySlice(axes.image, length.subtract(1), length);

          // convert back to an image
          var max = values.arrayProject([axes.band]).arrayFlatten([['effectRain', 'time']]);
          var MaxPos = max.select(1);
          //****end****//

          //different with Python since d starts from 1
          if (d <= 15){
            sumPrec = CHIRPS_DAILY.filter(ee.Filter.and(ee.Filter.gte('day', 0), ee.Filter.lte('day', d))).sum();
            sumEqs = EQcollection.filter(ee.Filter.and(ee.Filter.gte('day', 0), ee.Filter.lte('day', d))).sum();
          }
          else if (d > 15){
            sumPrec = CHIRPS_DAILY.filter(ee.Filter.and(ee.Filter.gte('day', d-16), ee.Filter.lte('day', d))).sum();
            sumEqs = EQcollection.filter(ee.Filter.and(ee.Filter.gte('day', d-16), ee.Filter.lte('day', d))).sum();
          }
          var Fzang = sumPrec.divide(sumEqs).rename('precipitation');
          var Fzang_all = Fzang.expression(
              "(b('precipitation') <= 0) ? 0" +
                ": (b('precipitation') > 1) ? 1" +
                  ": b('precipitation')"
          ).rename('Fzang_all');

    //TODO: fdrying to be implemented



          return Fzang_all.toDouble().set('day', d);
  }));



  return byDay;


}
//apply the constrains to PET
function apply_contrains_ET(image) {
  // """ CONSTRAINTS TO TRASNPIRATION"""
  var ET_DAILY_CANOPY = image.expression(
        'PotET_CANOPY * fm * ft * fg', {
        'PotET_CANOPY': image.select('PotET_CANOPY'),
        'fm': image.select('fm'),
        'ft': image.select('ft'),
        'fg': image.select('fg')
  });

  // """ CONSTRAINTS TO EVAPORATION"""
  var ET_DAILYSOIL = image.select('PotET_SOIL').multiply(image.select('Fdrying_all'));

  // """ TOTAL EVAPOTRASNPIRATION"""
  var TOTAL_DAILY_ET_Constrained = ET_DAILYSOIL.add(ET_DAILY_CANOPY).rename('AET');

  return TOTAL_DAILY_ET_Constrained;
}

function merge_bands(element) {
  return ee.Image.cat(element.get('primary'), element.get('secondary'));
}


/**
 * returns 1 if I_startDate is later in the year than I_endDate, 0 otherwise.
 */
function testIfStartDateLater(){
  var eeStartDate = ee.Date(I_startDate); // compares wheter I_startDate is later than I_endDate. Only MM-DD is compared, because I_startYear is just a placeholder
  var eeEndDate = ee.Date(I_endDate);

  return eeStartDate.difference(eeEndDate, 'day').gt(0);
}

//**************************************************
// (P_4) UIs
//**************************************************
//https://developers.google.com/earth-engine/ui_panels about the function of pannels

// Title
var toolTitle = ui.Label('Forward-ET Tool');
toolTitle.style().set('fontWeight', 'bold');
toolTitle.style().set('color', 'purple');
toolTitle.style().set({
  fontSize: '22px',
});

// Steps
var toolFilter = ui.Label('1) Select Filters', {fontWeight: 'bold'});
var toolRun = ui.Label('2) Run Model', {fontWeight: 'bold'});

// Textboxes
var I_startDateTexbox = ui.Textbox({
  placeholder: 'Start date',
  style: {width: '95px'}
}).setValue('2004-01-01');   // sets default value


var I_endDateTexbox = ui.Textbox({
  placeholder: 'End date',
  style: {width: '95px'}
}).setValue('2004-01-31');

// coordinate panel
var coordSectionLabel = ui.Label('Define AOI Centroid Coordinates', {color: 'gray'});

var latLabel = ui.Label('Latitude:');
var latBox = ui.Textbox({value:43.7929});
latBox.style().set('stretch', 'horizontal');

var lonLabel = ui.Label('Longitude:');
var lonBox = ui.Textbox({value:-122.8848});
lonBox.style().set('stretch', 'horizontal');

var latLonPanel = ui.Panel(
  [
    coordSectionLabel,
    ui.Panel([lonLabel, lonBox, latLabel, latBox],ui.Panel.Layout.Flow('horizontal'))
  ],
  null,
  {stretch: 'horizontal'}
);


// Run Button
var runButton = ui.Button('Run FORWARD-ET Model');
runButton.onClick(buildAndExportComposite);


//***Add Date textbox in displayed panel***//
var panelDate = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {width: '240px', margin: '-10px 0 0 0'}});

panelDate.add(I_startDateTexbox);
panelDate.add(I_endDateTexbox);

//***Add Run button in displayed panel***//

var panelRunButton = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {width: '240px', height: '120px', margin: '-10px 0 0 0'}});
panelRunButton.add(runButton);


// create pannel in which buttons are to be added
var panel = ui.Panel({layout: ui.Panel.Layout.flow('vertical'), style: {width: '300px', padding: '8px'}});

// add contents
panel.add(toolTitle);
panel.add(toolFilter);
panel.add(latLonPanel);

panel.add(ui.Label("Start/End date (YYYY-MM-DD)", {color: 'gray'}));
panel.add(panelDate);

panel.style({position: "top-left"});

panel.add(toolRun);
panel.add(panelRunButton);

// define panel for no image exception
var panelException = ui.Panel({style: {width: '240px'}});
panelException.setLayout(ui.Panel.Layout.flow());


// plot panel
var plotsPanelLabel = ui.Label('Time Series Plots', {fontWeight: 'bold', stretch: 'horizontal'});
plotsPanelLabel.style().set('color', 'purple');
plotsPanelLabel.style().set({
  fontSize: '22px',
});
var plotPanel = ui.Panel(null, null, {stretch: 'horizontal'});
var plotPanelParent = ui.Panel([plotsPanelLabel, plotPanel], null, {width: '300px'});

// map panel
var mapnew = ui.Map();
mapnew.style().set({cursor:'crosshair'});
mapnew.setOptions('HYBRID');
mapnew.setCenter(20, 5, 9);
mapnew.setZoom(2);
var processingLabel = ui.Label('Processing, please wait...', {shown:false, position: 'top-center'});
mapnew.add(processingLabel);


// plot time series for clicked point on map
mapnew.onClick(function(coords) {
  var x = coords.lon;
  var y = coords.lat;
  lonBox.setValue(x);
  latBox.setValue(y);
  buildAndPlotComposite();
  mapnew.setCenter(x, y, 16);
});

ui.root.clear();
ui.root.add(plotPanelParent);
ui.root.add(mapnew);
ui.root.add(panel);
