//define variable
var startYear = 2003;
var endYear = 2003;
var startDay = '-01-01';
var endDay = '-12-31';
var startDate = ee.Date(ee.String(ee.Number(startYear).toInt()).cat(startDay));
var endDate = ee.Date(ee.String(ee.Number(endYear).toInt()).cat(endDay));
//define a test area
var framed = ee.FeatureCollection('ft:1uZ5LdtjfS6fcrCpiKJw9vko91P6do30nFtt2FNp2', 'geometry');
//Map.addLayer(framed, {color: '800080'},'frame');
Map.centerObject(framed,4);

var I_sr = ['state_1km', 'SensorZenith', 'sur_refl_b01', 'sur_refl_b02']
var I_te = ['LST_Day_1km', 'QC_Day', 'Day_view_time', 'Day_view_angle', 'LST_Night_1km', 'QC_Night', 'Night_view_time', 'Night_view_angle', 'Emis_31', 'Emis_32']
var I_laifpar = ['Fpar', 'Lai', 'FparExtra_QC']
var I_ndvi = ['NDVI', 'SummaryQA']
var I_albedo = ['Albedo_BSA_shortwave', 'Albedo_WSA_shortwave', 'BRDF_Albedo_Band_Mandatory_Quality_shortwave']

var MODIS_SR = ee.ImageCollection("MODIS/006/MOD09GA")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_sr).map(maskSR),
    MODIS_TE = ee.ImageCollection("MODIS/006/MOD11A1")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_te).map(maskTE),
    MODIS_NDVI = ee.ImageCollection("MODIS/006/MOD13A1")
                  .filterDate(startDate, endDate).map(reprojMODIS).select(I_ndvi).map(maskNDVI),
    MODIS_LAIFPAR = ee.ImageCollection("MODIS/006/MCD15A3H")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_laifpar).map(maskLAIFPAR),
    MODIS_BRDFA = ee.ImageCollection("MODIS/006/MCD43A3")
        .filterDate(startDate, endDate).map(reprojMODIS).select(I_albedo).map(maskALBEDO);

// print (MODIS_SR)
// print (MODIS_TE)
// print (MODIS_NDVI)
// print (MODIS_LAIFPAR)
// print (MODIS_BRDFA)

var laimedian = MODIS_NDVI.select('NDVI').mean().clip(framed);
Map.addLayer(laimedian, {min:1000, max: 4000}, 'NDVI')

function reprojMODIS(image) {
  return image.reproject('EPSG:4326', null, 500).set('system:time_start', image.get('system:time_start'));
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


// Functions to keep desired pixels for MODIS products.

function maskSR(image) {
  // Select the QA band.
  var QA = image.select('state_1km');
  // Get the internal_cloud_algorithm_flag bit.
  var cloud = getQABits(QA, 0, 1, 'cloud');
  var cloudshadow = getQABits(QA, 2, 2, 'cloudshadow');
  // Return an image masking out cloudy areas.
  return image.updateMask(cloud.eq(0).and(cloudshadow.eq(0)));
}

function maskTE(image) {
  // Select the QA band.
  var QA = image.select('QC_Day');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 0, 7, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.eq(0));
}

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

function maskNDVI(image) {
  // Select the QA band.
  var QA = image.select('SummaryQA');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 0, 1, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.lte(1));
}

function maskALBEDO(image) {
  // Select the QA band.
  var QA = image.select('BRDF_Albedo_Band_Mandatory_Quality_shortwave');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA, 0, 0, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.eq(0));
}


var i=0;
var print_point = function(coords, map) {
  i++;

  var coord_array = Object.keys(coords).map(function (key) { return coords[key]; });
  var point = ee.Geometry.Point(coord_array);
  print('point ' + i, point);
  Map.addLayer(point);
  // put i somewhere near that point on the map
  // Plot the time series data at the ROI.
  print(ui.Chart.image.series(MODIS_LAIFPAR.select(['Lai']), point, ee.Reducer.mean(), 250)
    .setOptions({
      title: 'MODIS 16-DAY NDVI',
      lineWidth: 1,
      pointSize: 3,
  }));
};

Map.onClick(print_point);
