//define time variable
var startYear = 2003;
var endYear = 2003;
var startDay = '-01-01';
var endDay = '-01-02';
var startDate = ee.Date(ee.String(ee.Number(startYear).toInt()).cat(startDay));
var endDate = ee.Date(ee.String(ee.Number(endYear).toInt()).cat(endDay));
//define a test area
var framed = ee.FeatureCollection('ft:1uZ5LdtjfS6fcrCpiKJw9vko91P6do30nFtt2FNp2', 'geometry');
//Map.addLayer(framed, {color: '800080'},'frame');
Map.centerObject(framed,4);

//define the required bands for building the model
var I_sr = ['state_1km', 'SensorZenith', 'sur_refl_b01', 'sur_refl_b02']
var I_te = ['LST_Day_1km', 'QC_Day', 'Day_view_time', 'Day_view_angle',
            'Emis_31', 'Emis_32']
var I_laifpar = ['Fpar', 'Lai', 'FparExtra_QC']
var I_ndvi = ['NDVI', 'SummaryQA']
var I_albedo = ['Albedo_BSA_shortwave', 'Albedo_WSA_shortwave',
                'BRDF_Albedo_Band_Mandatory_Quality_shortwave']

var MODIS_SR = ee.ImageCollection("MODIS/006/MOD09GA")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_sr).map(maskSR),
    MODIS_TE_TERRA = ee.ImageCollection("MODIS/006/MOD11A1")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_te).map(renameTEterra).map(maskTEterra),
    MODIS_TE_AQUA = ee.ImageCollection("MODIS/006/MYD11A1")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_te).map(renameTEaqua).map(maskTEaqua),
    MODIS_NDVI = ee.ImageCollection("MODIS/006/MOD13A1")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_ndvi).map(maskNDVI),
    MODIS_LAIFPAR = ee.ImageCollection("MODIS/006/MCD15A3H")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_laifpar).map(maskLAIFPAR),
    MODIS_BRDFA = ee.ImageCollection("MODIS/006/MCD43A3")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_albedo).map(maskALBEDO);

//This function reprojects MODIS Sinusoidal projection to WGS84
function reprojMODIS(image) {
  return image.reproject('EPSG:4326', null, 500).set('system:time_start',
  image.get('system:time_start'));
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

//**************************************************
// Loads of functions to keep desired pixels for MODIS products.
//TODO: check QC bits
//**************************************************

//function for surface reflectance data masking
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
function maskNDVI(image) {
  // Select the QA band.
  var QA = image.select('SummaryQA');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 0, 1, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.lte(1));
}

//function for Albedo data masking
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
function maskTEterra(image) {
  // Select the QA band.
  var QA = image.select('QC_Day_terra');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 6, 7, 'internal_quality_flag_terra');
  // Return an image masking out cloudy areas.
  return image.addBands(internalQuality.neq(3).remap([1], [1], 0).rename('internal_quality_flag_terra'));
}

//mask aqua bands
function maskTEaqua(image) {
  // Select the QA band.
  var QA = image.select('QC_Day_aqua');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 6, 7, 'internal_quality_flag_aqua');
  // Return an image masking out cloudy areas.
  return image.addBands(internalQuality.neq(3).remap([1], [2], 0).rename('internal_quality_flag_aqua'));
}

//Create composite from Terra and Aqua
//TODO: revise after meeting
function maskTEcomposite(image) {
  var LST_mask = image.select('internal_quality_flag_terra').add(image.select('internal_quality_flag_aqua'))
                  .rename('MASK_comp');
  var LST_composite = ee.Image(LST_mask.eq(1).multiply(image.select('LST_Day_1km_terra')))
                      .add(ee.Image(LST_mask.eq(2).multiply(image.select('LST_Day_1km_aqua'))))
                      .add(ee.Image(LST_mask.eq(3)).multiply(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_time_terra'))).lt(0).multiply(image.select('LST_Day_1km_aqua')))
                      .add(ee.Image((image.select('Day_view_angle_aqua').subtract(image.select('Day_view_time_terra'))).gte(0).multiply(image.select('LST_Day_1km_terra'))))))
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
  return image.addBands(LST_mask).addBands(LST_composite).addBands(DVT_composite).addBands(EMIS_composite);
}

//**************************************************
// Loads of functions to keep desired pixels for MODIS products.
//TODO: cross-check equation validity
//**************************************************
//fn.GetShortWaveRadiationInputs
function get_SWR_inputs(image) {
  // TODO: temporal numeric constant: corresponding inputs will be defined later
  // """Units:
  //     Sunhours=hours*10
  //     Latitude=Degress
  //     Doy=Float64
  // """
  var Sunhours = ee.Image(10).divide(10).rename('Sunhours');
  var Latitude = ee.Image(1).rename('Latitude');
  var Doy = ee.Image(1).rename('Doy');
  return image.addBands(Sunhours.updateMask(Sunhours.neq(0))).addBands(Latitude).addBands(Doy);
}

//fn.GetLongWaveRadiationInputs
//TODO: add DASEMON_Tmax, DASEMON_Tmin
function get_LWR_inputs(image) {
  var DASEMON_Tmax = ee.Image(1).rename('TMAX');
  var DASEMON_Tmin = ee.Image(1).rename('TMIN');
  return image.addBands(DASEMON_Tmax).addBands(DASEMON_Tmin);
}

//fn.GetLAI
function get_LAI(image) {
  var LAI = image.select('Lai').rename('LAI').divide(10);
  return image.addBands(LAI);
}

//fn.GetMeanTemperature
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
  var LST = image.select('LST_comp');
  var ObsTime = image.select('DVT_comp');
  var Emissivity = image.select('EMIS_comp');
  var N = image.select('N');
  var AirMax = image.select('TMAX');
  var AirMin = image.select('TMIN');

  var Sunrise = image.select('Sunrise');
  var RSnet = image.select('RSnet');

  var c1 = -0.17;
  var d1 = 1.86;
  var BB = ee.Image(12).subtract(N.divide(2)).add(c1);
  var BBD = ObsTime.subtract(BB);
  var Tair_passtime = (AirMax.subtract(AirMin))
                      .multiply((ee.Image(Math.PI).multiply(BBD)
                      .divide((ee.Image(d1).multiply(2).add(N)))).sin())
                      .add(AirMin).rename('Tair_passtime');
  var Stef = 5.67e-008;
  var C = 0.261;
  var d = 7.77e-004;
  var EmissivityAir = (Tair_passtime.subtract(273.15)).pow(ee.Image.constant(2))
                      .multiply(ee.Image(d).multiply(-1)).exp().multiply(ee.Image(C).multiply(-1)).add(1);
  var Rlw_in = Tair_passtime.pow(ee.Image.constant(4)).multiply(EmissivityAir).multiply(Stef).rename('Rlw_in');
  var Rlw_out = LST.pow(4).multiply(Emissivity).multiply(Stef).rename('Rlw_out');

  var RLnet = Rlw_in.subtract(Rlw_out);

  var t = ObsTime.subtract(Sunrise).rename('t');

  var J = N.multiply(2).divide(t.multiply(Math.PI).divide(N).cos().multiply(24));
  var RSnetInstant = RSnet.divide(J);
  var RnDaily = (RSnetInstant.add(RLnet)).multiply(J).rename('RnDaily');
  return image.addBands(RnDaily).addBands(Tair_passtime);
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

function PET(image) {
  var Tmean = image.select('Tmean');
  var RnCanopy = image.select('Rn_Canopy');
  var Rn_SOIL = image.select('Rn_SOIL');
  //DEM to be changed later
  var Altitude = ee.Image(1);

  //calPsycometricConstant
  var P = ee.Image(101.3).multiply(((ee.Image(293).subtract(ee.Image(0.0065).multiply(Altitude)))
          .divide(293)).pow(ee.Image.constant(5.26)));
  var Lambda = ee.Image(2454).subtract(ee.Image(2.4).multiply(Tmean.subtract(20)));
  var cp = 1013;
  var epsi = 0.622;
  var psi = (P.multiply(cp)).divide(Lambda.multiply(epsi));

  //calSatVapPres
  //es = 10 * (0.061121 * 2.718281828**(17.502 * Tmean / (240.97 + Tmean))
  //               * (1.0007 + (3.46 * 10**(-8) * 100)))
  var es = ee.Image(10).multiply(ee.Image(0.061121)
            .multiply(ee.Image(2.718281828).pow(Tmean.multiply(17.502).divide(Tmean.add(240.97))))
            .multiply((ee.Image(3.46).multiply(ee.Image(10).pow(-8)).multiply(100)).add(1.0007)))

  //calSlopVapPresCurv
  var s = (Lambda.multiply(18).multiply(1000).multiply(es))
          .divide(ee.Image(8.3144).multiply(Tmean.add(273).pow(2)));

  //calPotET_CANOPY
  //TODO: Check whether s changed: Etpc = alphaPT * (s / (s + psi)) * RnCanopy
  var alphaPT = 1.26;
  var Etpc = ee.Image(alphaPT).multiply(s.divide(s.add(psi))).multiply(RnCanopy);

  //calPotET_SOIL
  var G = 0;
  var Etps = ee.Image(alphaPT).multiply(s.divide(s.add(psi))).multiply(Rn_SOIL.subtract(G));

  //PET
  PotEVTR = Etpc.add(Etps).rename('PET');

  return image.addBands(PotEVTR);
}

function merge_bands(element) {
  return ee.Image.cat(element.get('primary'), element.get('secondary'));
}

// Define the join and filter
var Join = ee.Join.inner();
var FilterOnStartTime = ee.Filter.equals({
                                        'leftField': 'system:time_start',
                                        'rightField': 'system:time_start'
                                        });

var TE_Joined = ee.ImageCollection(Join.apply(MODIS_TE_TERRA, MODIS_TE_AQUA, FilterOnStartTime)).map(merge_bands).map(maskTEcomposite);

print (TE_Joined);

// Join the collections, passing entries through the filter
var AL_Joined = ee.ImageCollection(Join.apply(MODIS_BRDFA, MODIS_LAIFPAR, FilterOnStartTime)).map(merge_bands);

print (AL_Joined);

var FinalDataset = ee.ImageCollection(Join.apply(AL_Joined, TE_Joined, FilterOnStartTime)).map(merge_bands);

print (FinalDataset);

var PotEVTR = FinalDataset.map(calc_albedo).select('Albedo','Lai', 'EMIS_comp', 'LST_comp', 'DVT_comp')
            .map(get_SWR_inputs)
            .map(calc_ShortWaveRadiation)
            .map(get_LWR_inputs)
            .map(calc_LongWaveRadiation)
            .map(get_meantemp)
            .map(get_LAI)
            .map(split_cannopy)
            .map(PET);
print (PotEVTR);
var disp = PotEVTR.select('PET').mean().clip(framed);
Map.addLayer(disp, {min:0, max: 1000}, 'PET');

// var disp = TE_Joined.select('LST_Day_1km_aqua').mean().clip(framed);
// Map.addLayer(disp, {min:270, max: 300}, 'lst');

// var disp = TE_Joined.select('LST_Day_1km_terra').mean().clip(framed);
// Map.addLayer(disp, {min:270, max: 300}, 'lst');

var disp = TE_Joined.select('LST_comp').mean().clip(framed);
Map.addLayer(disp, {min:290, max: 320}, 'lst');
