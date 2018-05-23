//define variable
var startYear = 2016;
var endYear = 2016;
var startDay = '-01-01';
var endDay = '-01-31';
var startDate = ee.Date(ee.String(ee.Number(startYear).toInt()).cat(startDay));
var endDate = ee.Date(ee.String(ee.Number(endYear).toInt()).cat(endDay));
//define a test area
var framed = ee.FeatureCollection('ft:1uZ5LdtjfS6fcrCpiKJw9vko91P6do30nFtt2FNp2', 'geometry');
//Map.addLayer(framed, {color: '800080'},'frame');
Map.centerObject(framed,4);
var I_gfs = ['temperature_2m_above_ground','u_component_of_wind_10m_above_ground','v_component_of_wind_10m_above_ground',
'relative_humidity_2m_above_ground','total_precipitation_surface','downward_shortwave_radiation_flux'];

var GFS = ee.ImageCollection("NOAA/GFS0P25").filterDate(startDate, endDate).select(I_cfsv2).map(dewtemp).map(wd);
print (GFS, 'GFS');

//dew point temperature calculation: unit-Celsius
//follow: https://iridl.ldeo.columbia.edu/dochelp/QA/Basic/dewpoint.html
function dewtemp(image) {
  var dew_point = image.expression(
      ' T - ((100 - RH)/5.)', {
        'T': image.select('temperature_2m_above_ground'),
        'RH': image.select('relative_humidity_2m_above_ground')
  }).rename('dewt');
  return image.addBands(dew_point);
}
//wind direction as the wind is coming from: unit-degrees
function wd(image) {
  var wd_degree = image.expression(
      'atan2(u ,v) * (180 / pi ) + 180', {
        'u': image.select('u_component_of_wind_10m_above_ground'),
        'v': image.select('v_component_of_wind_10m_above_ground'),
        'pi': Math.PI
  }).rename('wd');

  return image.addBands(wd_degree);
}

// //NCEP Climate Forecast System Version 2, 6-Hourly Products
// var I_cfsv2 = ['Temperature_height_above_ground','u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground'];
// var CFSV2 = ee.ImageCollection("NOAA/CFSV2/FOR6H").filterDate(startDate, endDate).select(I_cfsv2).map(t2m).map(wdcfsv);
// print (CFSV2);


// //CFSV2: Temperature 2m above ground, 6-hour interval, unit: k, convert to ceicius
// function t2m(image) {
//   var t2m_c = image.expression('temp - 273.15', {'temp': image.select('Temperature_height_above_ground')}).rename('t2mc');
//   return image.addBands(t2m_c);
// }

// //CFSV2: wind direction from u and v components
// function wdcfsv(image) {
//   var wd_degree = image.expression(
//       'atan2(u ,v) * (180 / pi ) + 180', {
//         'u': image.select('u-component_of_wind_height_above_ground'),
//         'v': image.select('v-component_of_wind_height_above_ground'),
//         'pi': Math.PI
//   }).rename('wd');

//   return image.addBands(wd_degree);
// }
