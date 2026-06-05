// 0️⃣ Study area
Map.addLayer(
  studyArea.style({
    color: '000000',   // black outline
    fillColor: '00000000', // transparent fill
    width: 2
  }),
  {},
  'Study Area'
);
Map.centerObject(studyArea, 13);
Map.addLayer(table, {color:'orange'}, 'p2013');
Map.centerObject(table, 13);
// 1️⃣ Minimal Landsat 8 cloud mask
function maskClouds(image) {

  var qa = image.select('QA_PIXEL');

  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var shadow = qa.bitwiseAnd(1 << 4).eq(0);

  return image.updateMask(cloud.and(shadow));
}

// 2️⃣ Apply reflectance scaling
function applyScaleFactors(image) {

  var optical = image.select([
    'SR_B2','SR_B3','SR_B4',
    'SR_B5','SR_B6','SR_B7'
  ])
  .multiply(0.0000275)
  .add(-0.2);

  return optical.copyProperties(image, image.propertyNames());
}

// 3️⃣ Load Landsat 8 (2013 only)
var landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(studyArea)
  .filterDate('2013-01-01','2013-12-31')
  .map(maskClouds)
  .map(applyScaleFactors)
  .select([
    'SR_B2','SR_B3','SR_B4',
    'SR_B5','SR_B6','SR_B7'
  ]);

// 4️⃣ Create median composite
var composite = landsat8.median().clip(studyArea);

// 5️⃣ Strong spatial interpolation to fill gaps
var gapFill = composite.unmask(
  composite.focal_mean({
    radius: 3,        // larger neighborhood
    units: 'pixels'
  })
);

// 6️⃣ Second pass interpolation (for larger gaps)
var gapFill2 = gapFill.unmask(
  gapFill.focal_mean({
    radius: 5,
    units: 'pixels'
  })
);

// 7️⃣ Final composite
var annualComposite = gapFill2.clip(studyArea);

// 8️⃣ Display
Map.addLayer(
  annualComposite,
  {bands:['SR_B4','SR_B3','SR_B2'], min:0, max:0.3},
  'Landsat 8 RGB 2013'
);

// 9️⃣ Diagnostic
print('Number of Landsat scenes used:', landsat8.size());
// 4️⃣ Compute indices for Landsat
function addIndices(image){

  // ==========================
  // Vegetation
  // ==========================
  var NDVI = image.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI');

  // ==========================
  // Built-up
  // ==========================
  var MNDBI = image.normalizedDifference(['SR_B6','SR_B5']).rename('MNDBI');

  // ==========================
  // Soil (BSI)
  // ==========================
  var BSI = image.expression(
    '((SWIR+RED)-(NIR+BLUE))/((SWIR+RED)+(NIR+BLUE))',{
      'SWIR': image.select('SR_B6'),
      'RED': image.select('SR_B4'),
      'NIR': image.select('SR_B5'),
      'BLUE': image.select('SR_B2')
  }).rename('BSI');

  // ==========================
  // Plastic Index (PI → PI1 → PI2)
  // ==========================
  var PI = image.expression(
    'NIR/(NIR+RED)',{
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4')
  }).rename('PI');

  // NDVI correction
  var PI1 = PI.where(NDVI.gt(0), PI.subtract(NDVI)).rename('PI1');

  // MNDBI correction
  var PI2 = PI1.where(MNDBI.gt(0), PI1.subtract(MNDBI)).rename('PI2');

  // ==========================
  // 🔧 BSI-based correction (DATA-DRIVEN)
  // ==========================

  // Clamp BSI (based on your distribution)
  var BSI_c = BSI.min(0.2);

  // Attenuation factor
  var attenuation = ee.Image(1).subtract(BSI_c.multiply(0.4));

  // Apply conditional correction
  var PI2_corrected = PI2.where(
    BSI_c.gt(0.05),
    PI2.multiply(attenuation)
  ).rename('PI2_corrected');

  // Water
  var MNDWI = image.normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI');

  // SWIR ratio
  var SWIRratio = image.select('SR_B6')
                       .divide(image.select('SR_B5'))
                       .rename('SWIRratio');

  // Blue–SWIR Plastic Separation
  var PSI = image.normalizedDifference(['SR_B2','SR_B6'])
                 .rename('PSI');


  // SWIR difference
  var SWIRdiff = image.normalizedDifference(['SR_B6','SR_B7'])
                      .rename('SWIRdiff');


  // ==========================
  // Return all bands
  // ==========================
  return image.addBands([
    NDVI,
    MNDBI,
    PI,
    PI1,
    PI2,
    PI2_corrected,   // ✅ NEW corrected index
    MNDWI,
    PSI,
    BSI,
    SWIRdiff
  ]);
}

// 5️⃣ Create image stack with indices
var image_stack = addIndices(annualComposite);

// 6️⃣ Map indices
Map.addLayer(image_stack.select('NDVI'),
  {min:-1, max:1, palette:['brown','yellow','green']},
  'NDVI');

Map.addLayer(image_stack.select('MNDBI'),
  {min:-1, max:1, palette:['white','grey','black']},
  'MNDBI');

Map.addLayer(image_stack.select('PI2'),
  {min:0, max:1, palette:['white','grey','blue']},
  'Plastic Index (PI2)');

Map.addLayer(image_stack.select('MNDWI'),
  {min:-1, max:1, palette:['brown','white','blue']},
  'MNDWI');

Map.addLayer(image_stack.select('BSI'),
  {min:-1, max:1, palette:['black','yellow','brown']},
  'BSI');
  // 7️⃣ Bands for classification
var bandsForClassification = [
  'SR_B2','SR_B3','SR_B4',
  'SR_B5','SR_B6','SR_B7','MNDWI','BSI',
  'NDVI','MNDBI','PI2_corrected','PSI','SWIRdiff'
];

// 8️⃣ Merge training polygons (assumes polygons exist)
var polygons = Vegetation.merge(Waterbody)
                         .merge(Debris)
                         .merge(Builtup)
                         .merge(Bareground);

// Split 75% training / 25% validation
var sample = polygons.randomColumn();

var trainingPolygons = sample.filter(ee.Filter.lte('random', 0.75));
var validationPolygons = sample.filter(ee.Filter.gt('random', 0.75));


// Print number of polygons
print('Training polygons:', trainingPolygons.size());
print('Validation polygons:', validationPolygons.size());

// Print class distribution
print('Training class distribution:',
      trainingPolygons.aggregate_histogram('Class'));

print('Validation class distribution:',
      validationPolygons.aggregate_histogram('Class'));


// Sample image values for training & validation
var training = image_stack.select(bandsForClassification).sampleRegions({
  collection: trainingPolygons,
  properties: ['Class'],
  scale: 30,
  tileScale: 2
});

var validation = image_stack.select(bandsForClassification).sampleRegions({
  collection: validationPolygons,
  properties: ['Class'],
  scale: 30,
  tileScale: 2
});


// 9️⃣ Train Random Forest classifier
var RFclassifier = ee.Classifier.smileRandomForest(300)
  .train({
    features: training,
    classProperty: 'Class',
    inputProperties: bandsForClassification
  });


// 10️⃣ Classify image
var Classified = image_stack.select(bandsForClassification)
  .classify(RFclassifier)
  .clip(studyArea)
   .focal_mode({
    radius: 2,
    units: 'pixels'
  });

// 11️⃣ Display classification
var palette = ['006400','1E3A8A','FF1A1A','FF8C00','A0522D']; // Vegetation, Waterbody, Built-up, Debris, Baregroundl
var classNames = ['Vegetation','Waterbody','Debris','Builtup','Bareground'];
Map.addLayer(Classified, {min:1, max:5, palette:palette}, 'RF Classified');


// 12️⃣ Add Legend
function addLegend(classNames, classPalette) {
  var legend = ui.Panel({style:{position:'bottom-right', padding:'8px 15px', backgroundColor:'white'}});
  legend.add(ui.Label('Classification Legend', {fontWeight:'bold', fontSize:'14px', margin:'0 0 4px 0'}));
  for (var i=0;i<classNames.length;i++){
    var colorBox = ui.Label({style:{backgroundColor:classPalette[i], padding:'8px', margin:'0 0 4px 0'}});
    var label = ui.Label({value:classNames[i], style:{margin:'0 0 4px 6px'}});
    legend.add(ui.Panel([colorBox,label], ui.Panel.Layout.Flow('horizontal')));
  }
  Map.add(legend);
}
addLegend(classNames, palette);

// 12️⃣ Create visualization image
var classifiedPNG = Classified.visualize({
  min: 1,
  max: 5,
  palette: palette
});

// 13️⃣ Create a bounded region for thumbnail export
var exportRegion = studyArea.geometry().bounds();

// 14️⃣ Generate PNG download link
var pngURL = classifiedPNG.getThumbURL({
  region: exportRegion,
  dimensions: 1500,   // Reduced to avoid size-limit errors
  format: 'png'
});

// 15️⃣ Print download link
print('PNG Download Link:', pngURL);

// ========================================================
// 14️⃣ Variable Importance
var explain = RFclassifier.explain();
var importance = ee.Feature(null, ee.Dictionary(explain).get('importance'));

var varImportanceChart = ui.Chart.feature.byProperty(importance)
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Random Forest Variable Importance 2013',
    legend: {position: 'none'},
    hAxis: {title: 'Bands / Indices'},
    vAxis: {title: 'Importance'}
  });

print(varImportanceChart);


// 13️⃣ Accuracy assessment
// ========================================================
// 13️⃣ Accuracy Assessment (Class 0 removed)

// --- 13A: Training Accuracy ---
var trainAccuracy = RFclassifier.confusionMatrix();
print('Training Confusion Matrix', trainAccuracy);
print('Training Overall Accuracy', trainAccuracy.accuracy());

// --- 13B: Correctly Classified Training Samples per Class ---
// Convert training confusion matrix to ee.Array
var trainConfArray = ee.Array(trainAccuracy.array());

// Extract diagonal skipping class 0 (indices 1–5)
var diagTrain = ee.List.sequence(1, 5).map(function(i){
  return trainConfArray.get([i, i]);
});

// Training class names
var classNamesTrain = [
  'Vegetation',
  'Waterbody',
  'Debris',
  'Builtup',
  'Bareground'
];

// Training chart
var trainingClassChart = ui.Chart.array.values({
  array: ee.Array(diagTrain),
  axis: 0,
  xLabels: classNamesTrain
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Correctly Classified Training Samples per Class',
  hAxis: {title: 'Land Cover Class'},
  vAxis: {title: 'Number of Correct Samples'},
  legend: {position: 'none'},
  colors: ['#1f78b4']
});
print(trainingClassChart);

// --- 13C: Validation Accuracy ---
var validated = validation.classify(RFclassifier);
var validationAccuracy = validated.errorMatrix('Class','classification');
print('Validation Confusion Matrix', validationAccuracy);
print('Validation Overall Accuracy', validationAccuracy.accuracy());
print('Kappa Coefficient', validationAccuracy.kappa());

// Producer’s Accuracy and User’s Accuracy (skip class 0)
var PA = validationAccuracy.producersAccuracy().slice(1);
var UA = validationAccuracy.consumersAccuracy().slice(1);
print('Producer Accuracy per Class', PA);
print('User Accuracy per Class', UA);

// --- 13D: Correctly Classified Validation Samples per Class ---
// Convert validation confusion matrix to array
var valConfArray = ee.Array(validationAccuracy.array());

// Extract diagonal skipping class 0 (indices 1–5)
var diagVal = ee.List.sequence(1, 5).map(function(i){
  return valConfArray.get([i, i]);
});

// Validation class names
var classNamesVal = [
  'Vegetation',
  'Waterbody',
  'Debris',
  'Builtup',
  'Bareground'
];

// Validation chart
var validationClassChart = ui.Chart.array.values({
  array: ee.Array(diagVal),
  axis: 0,
  xLabels: classNamesVal
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Correctly Classified Validation Samples per Class',
  hAxis: {title: 'Land Cover Class'},
  vAxis: {title: 'Number of Correct Samples'},
  legend: {position: 'none'},
  colors: ['red']
});
print(validationClassChart);


// 15️⃣ Export classified map with projection
//Export.image.toDrive({
//  image: Classified,
 // description: 'Landsat8_RF_2013',
 // region: studyArea,
 // scale: 30,
//  maxPixels: 1e13,
 //crs: 'EPSG:4326',   // WGS84 (lat/long)
//});
// ===============================


// ===============================
// 0️⃣ INPUTS
// ===============================
// Classified → RF classified image (debris = class 3)
// plastic, nonplastic → validation FeatureCollections
// studyArea → boundary

// ===============================
// 1️⃣ Predicted debris (binary)
// ===============================
var Debris = Classified.eq(3).rename('predicted');

var Debris_vis = Debris.focal_max({
  radius: 1,
  units: 'pixels'
});

// ===============================
// 2️⃣ MERGE VALIDATION
// ===============================
var validationPoints = plastic.merge(nonplastic);

print('Validation distribution (6/7):',
  validationPoints.aggregate_histogram('Class'));

// ===============================
// 3️⃣ CREATE BINARY CLASS
// ===============================
var validationBinary = validationPoints.map(function(f) {

  var cls = ee.Number.parse(f.get('Class'));

  return f.set('Class_bin',
    ee.Algorithms.If(cls.eq(6), 1, 0)
  );
});

print('Binary distribution (1/0):',
  validationBinary.aggregate_histogram('Class_bin'));

// ===============================
// 4️⃣ SAMPLE USING BINARY
// ===============================
var samples = Classified.sampleRegions({
  collection: validationBinary,
  properties: ['Class_bin'],
  scale: 10,
  geometries: true
});

// ===============================
// 5️⃣ CONFUSION TYPES
// ===============================
var samplesWithError = samples.map(function(f) {

  var actual = ee.Number(f.get('Class_bin'));
  var predictedClass = ee.Number(f.get('classification'));

  var pred = predictedClass.eq(3);
  pred = ee.Number(pred);

  var errorType = ee.Algorithms.If(
    pred.eq(1).and(actual.eq(1)), 1,   // TP
    ee.Algorithms.If(
      pred.eq(1).and(actual.eq(0)), 2, // FP
      ee.Algorithms.If(
        pred.eq(0).and(actual.eq(1)), 3, // FN
        4 // TN
      )
    )
  );

  return f.set({
    'predicted_binary': pred,
    'errorType': errorType
  });
});

// ===============================
// 6️⃣ DERIVED LAYERS
// ===============================
var ActualPlastic = samplesWithError.map(function(f) {
  var e = ee.Number(f.get('errorType'));
  return f.set('actual_plastic', e.eq(1).or(e.eq(3)));
});

var CriticalErrors = samplesWithError.map(function(f) {
  var e = ee.Number(f.get('errorType'));
  return f.set('critical_error', e.eq(2).or(e.eq(3)));
});

// ===============================
// 7️⃣ CONVERT TO IMAGE
// ===============================
var ErrorMap = samplesWithError.reduceToImage({
  properties: ['errorType'],
  reducer: ee.Reducer.first()
});

var ActualPlastic_img = ActualPlastic.reduceToImage({
  properties: ['actual_plastic'],
  reducer: ee.Reducer.first()
});

var CriticalErrors_img = CriticalErrors.reduceToImage({
  properties: ['critical_error'],
  reducer: ee.Reducer.first()
});

// ===============================
// 8️⃣ SMOOTHING (visual only)
// ===============================
var ErrorMap_vis = ErrorMap.focal_max({radius:4, units:'pixels'});
var ActualPlastic_vis = ActualPlastic_img.focal_max({radius:4, units:'pixels'});
var CriticalErrors_vis = CriticalErrors_img.focal_max({radius:4, units:'pixels'});

// ===============================
// 9️⃣ VALIDATION BINARY (ONE LAYER)
// ===============================
var validationBinary_img = validationBinary.reduceToImage({
  properties: ['Class_bin'],
  reducer: ee.Reducer.first()
});

var validationBinary_vis = validationBinary_img.focal_max({
  radius: 2,
  units: 'pixels'
});

// ===============================
// 🔟 DISPLAY
// ===============================

// Predicted debris
Map.addLayer(Debris_vis.clip(studyArea), {
  min: 0,
  max: 1,
  palette: [
    '000000', // black → no debris
    'FF0000'  // red → predicted debris
  ]
}, 'Predicted Debris');

// Confusion map
Map.addLayer(ErrorMap_vis.clip(studyArea), {
  min: 1,
  max: 4,
  palette: [
    '00FF00', // green → True Positive (TP)
    'FFFF00', // yellow → False Positive (FP)
    'FF0000', // red → False Negative (FN)
    '0000FF'  // blue → True Negative (TN)
  ]
}, 'Confusion Map');

// Actual plastic
Map.addLayer(
  ActualPlastic_vis.eq(1).selfMask().clip(studyArea),
  {palette: ['00FF00']},  // green → plastic (TP + FN)
  'Actual Plastic'
);

// Critical errors
Map.addLayer(
  CriticalErrors_vis.eq(1).selfMask().clip(studyArea),
  {palette: ['FF0000']},  // red → critical errors (FP + FN)
  'Critical Errors'
);
// ✅ SINGLE validation binary layer
Map.addLayer(validationBinary_vis.clip(studyArea), {
  min: 0,
  max: 1,
  palette: [
    '0000FF', // blue → non-plastic (0)
    '00FF00'  // green → plastic (1)
  ]
}, 'Validation Binary');

// ===============================
// CENTER MAP
// ===============================
Map.centerObject(studyArea, 13);
// ===============================
// CALCULATE CONFUSION COUNTS
// ===============================
var TP = samplesWithError.filter(ee.Filter.eq('errorType', 1)).size();
var FP = samplesWithError.filter(ee.Filter.eq('errorType', 2)).size();
var FN = samplesWithError.filter(ee.Filter.eq('errorType', 3)).size();
var TN = samplesWithError.filter(ee.Filter.eq('errorType', 4)).size();

// Print counts
print('TP:', TP);
print('FP:', FP);
print('FN:', FN);
print('TN:', TN);

// ===============================
// METRICS
// ===============================

// Precision = TP / (TP + FP)
var Precision = TP.divide(TP.add(FP));

// Recall = TP / (TP + FN)
var Recall = TP.divide(TP.add(FN));

// F1 Score
var F1 = Precision.multiply(Recall)
  .multiply(2)
  .divide(Precision.add(Recall));

// Print results
print('Precision:', Precision);
print('Recall:', Recall);
print('F1 Score:', F1);
