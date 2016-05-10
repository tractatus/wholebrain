#' Generate a web-map from a sample section
#'
#' Generates multiple image tiles in JPEG format.
#' @param img input a character vector consisting of the full path name to 16-bit raw tif image files.
#' @param alpha filename of the stitched output, default is the directory where the tiles are situated with the added prefix stitched_
#' @param beta type of motorized microscope stage. Will define the order of acquisition. Arguments are either row.by.row or default snake.by.row.
#' @param cells the order each tile is acquired. Default is left.&.up which starts from bottom left corner and moves up. Other options include right.&.up, left.&.down and right.&.down.
#' @param outputfolder name of output folder to save the image in. 
#' @param verbose boolean value. If true diagnostic output is written to the R console. Deafult is true.
#' @examples
#' #folder where image tiles are stored
#'makewebmap('/Users/danielfurth/Documents/mouse_atlas/cropped_D1_59_s4.lif_TileScan_006_Merging001_ch00.tif', alpha=60, beta=10, verbose=TRUE)
#' #stitch images
#' stitch(images, type = 'snake.by.row', order = 'left.&.up', tilesize=2048, overlap=0.1, show.image=TRUE)

makewebmap<-function(img, filter, alpha=1, beta=0, scale = 0.64, bregmaX = 9280, bregmaY = 3120, fluorophore = 'Rabies-EGFP', enable.drawing=TRUE, verbose=FALSE, registration = FALSE, AP=0,ML=0,DV=0){
    file <- as.character(img)[1]
    if(!file.exists(file))
    stop(file, "not found")
    
    
    ## expand path
    file <- path.expand(file)
    outputfile<-basename(file)
    outputfile<-sub("^([^.]*).*", "\\1", outputfile)
    create.output.directory(paste('Web', outputfile, sep='_'))
    setwd(paste('Web', outputfile, sep='_'))
    create.output.directory(paste('Tiles', outputfile, sep='_'))
    setwd('../')
    #if(is.null(overlap)){overlap<-(-999)}
    verbose<-as.numeric(verbose)


    tiled_imagefolder<-paste('Tiles', outputfile, sep='_')
    cat('Registering to atlas')
    if(!missing(filter)){
      alpha<-filter$Max
      beta<-filter$Min
    }
    a <- .Call("createWeb", file, alpha, beta, verbose, outputfile)


    headerP01<-'<!DOCTYPE html>
<html>
<head>
    <title>'


headerP02<-"</title>

    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no\">

    <link href=\"http://maxcdn.bootstrapcdn.com/font-awesome/4.1.0/css/font-awesome.min.css\" rel=\"stylesheet\">
    <link rel=\"stylesheet\" href=\"./dist/leaflet.css\" />
    <!-- <link rel=\"stylesheet\" href=\"http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.css\" /> -->
    <!--[if lte IE 8]><link rel=\"stylesheet\" href=\"http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.ie.css\" /><![endif]-->
  <link href=\"http://www.openbrainmap.org/cdn/css/lato.css\" rel=\"stylesheet\" type=\"text/css\">
    <link rel=\"stylesheet\" href=\"http://www.openbrainmap.org/cdn/css/leaflet-sidebar.css\" />
   <link rel=\"stylesheet\" type=\"text/css\" media=\"all\" href=\"http://www.openbrainmap.org/cdn/css/whhg.css\" />

    <link href=\"http://www.openbrainmap.org/cdn/css/ionicons.min.css\" rel=\"stylesheet\" type=\"text/css\">
      <link rel=\"stylesheet\" href=\"http://www.openbrainmap.org/cdn/css/Control.MiniMap.min.css\" /> 
      <link href=\"http://netdna.bootstrapcdn.com/bootstrap/3.1.1/css/bootstrap.min.css\" rel=\"stylesheet\">

    <style type=\"text/css\">
    
    .leaflet-control-layers {
  box-shadow: 0 1px 5px rgba(0,0,0,0.4);
  background: #383838;
  border-radius: 10px;
  color: white;
  font-family: 'Lato', sans-serif;
  font-weight: 300;
  }

  .leaflet-control-layers .hr{
  border-color: black;
  }

  
  .marker-delete-button {
    background-color: #64A375;
    -moz-transition: all 0.15s linear;
    -webkit-transition: all 0.15s linear;
    -ms-transition: all 0.15s linear;
    transition: all 0.15s linear;
    border-radius: 6px;
    background-color: transparent;
    color: white;
    border: 2px solid white;
  }
    .the-icons {
        list-style: none outside none;
        margin-left: 0;
    }
    .the-icons li {
        float: left;
        line-height: 25px;
        width: 25%;
    }
    .the-icons i:hover {
        background-color: rgba(255, 0, 0, 0.25);
    }
  .btn-file {
    position: relative;
    overflow: hidden;
  }
  .btn-file input[type=file] {
    position: absolute;
    top: 0;
    right: 0;
    min-width: 100%;
    min-height: 100%;
    font-size: 100px;
    text-align: right;
    filter: alpha(opacity=0);
    opacity: 0;
    background: red;
    cursor: inherit;
    display: block;
  }
  input[readonly] {
    background-color: white !important;
    cursor: text !important;
  }
  
  .leaflet-container {
          background-color: #64A375;
      }
    
        html, body { width: 100%; height: 100%; margin: 0; background-color: #212121;}
    hr {
        padding: 0;
        border: none;
        border-top: medium double #fff;
        color: #fff;
        text-align: center;
    }
    hr:after {
        content: \"&#167;\";
        display: inline-block;
        position: relative; 
        top: -0.7em;  
        font-size: 1.5em;
        padding: 0 0.25em;
        background: #64A375;
    }
        #map, #container { width: 100%; height: 100%; }
        #map { float: left; }
        #container { float: right;}
        #container .map { width: 100%; height: 50%;}
    .leaflet-popup-content {
    margin: 13px 19px;
    line-height: 1.4;
    font-family: 'Lato', sans-serif;
    font-weight: 300;
    color: white;
    background-color: #64A375;
    }
    .leaflet-container a.leaflet-popup-close-button {
    position: absolute;
    top: 0;
    right: 0;
    padding: 4px 4px 0 0;
    text-align: center;
    width: 18px;
    height: 14px;
    font: 16px/14px Tahoma, Verdana, sans-serif;
    color: #ffffff;
    text-decoration: none;
    font-weight: bold;
    background: transparent;
    }
    .leaflet-popup-content-wrapper, .leaflet-popup-tip {
    background: #64A375;
    box-shadow: 0 3px 14px rgba(0,0,0,0.4);
    }
    .leaflet-container {
        background: #000000;
    }
    .leaflet-container .leaflet-control-attribution, .leaflet-container .leaflet-control-scale {
      font-size: 14px;
      font-family: 'Lato', sans-serif;
      font-weight: 300;
      color: white;
      border: 2px solid white;
      background-color: #64A375;
      border-radius: 6px;
      /* unvisited link */
    }
    a:link {
        color: #ddd;
    }

    /* visited link */
    a:visited {
        color: #ddd;
    }

    /* mouse over link */
    a:hover {
        color: #fff;
    }

    /* selected link */
    a:active {
        color: #fff;
    }
    
    a.leaflet-control-zoom-in {
       color: #000;
    }
    
    a.leaflet-control-zoom-out {
       color: #000;
    }

    input[type=\"radio\"], input[type=\"checkbox\"] {
      line-height: normal;
      margin: 4px 0px 0px;
      height: 1.6em;
      width: 1.6em;
    }

    </style>
    <script src=\"http://code.jquery.com/jquery-1.11.0.min.js\"></script>
  <script src=\"https://rawgit.com/jeresig/jquery.hotkeys/master/jquery.hotkeys.js\"></script>
  <script src=\"http://netdna.bootstrapcdn.com/bootstrap/3.1.1/js/bootstrap.min.js\"></script>

  <script src=\"http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.js\"></script>
    <script src=\"http://cdn.openbrainmap.org/0.0.1/src/Leaflet.draw.js\"></script>
  <link rel=\"stylesheet\" href=\"http://cdn.openbrainmap.org/0.0.1/dist/leaflet.draw.css\" />
  <link rel=\"stylesheet\" href=\"http://cdn.openbrainmap.org/0.0.1/src/easy-button.css\" />

  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/Toolbar.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/Tooltip.js\"></script>

  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/ext/GeometryUtil.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/ext/LatLngUtil.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/ext/LineUtil.Intersect.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/ext/Polygon.Intersect.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/ext/Polyline.Intersect.js\"></script>


  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/DrawToolbar.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/handler/Draw.Feature.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/handler/Draw.SimpleShape.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/handler/Draw.Polyline.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/handler/Draw.Circle.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/handler/Draw.Marker.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/handler/Draw.Polygon.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/draw/handler/Draw.Rectangle.js\"></script>


  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/EditToolbar.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/handler/EditToolbar.Edit.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/handler/EditToolbar.Delete.js\"></script>

  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/Control.Draw.js\"></script>

  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/handler/Edit.Poly.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/handler/Edit.SimpleShape.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/handler/Edit.Circle.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/handler/Edit.Rectangle.js\"></script>
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/edit/handler/Edit.Marker.js\"></script>
  
  <script src=\"http://cdn.openbrainmap.org/0.0.1/src/easy-button.js\"></script>

  <script src=\"https://ajax.googleapis.com/ajax/libs/jquery/2.1.4/jquery.min.js\"></script>
    
    <script src=\"http://www.jqueryscript.net/demo/jQuery-Plugin-To-Print-Any-Part-Of-Your-Page-Print/jQuery.print.js\"></script>
  
</head>
<body>"

JSONfilename<-'http://www.openbrainmap.org/cdn/js/mupp.js'
JSONdata<-paste('<script src=\"', JSONfilename,'\" type=\"text/javascript\"></script>',sep='') 

footer<-sprintf('<div id=\"map\"></div>

    <script src=\"http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.js\"></script>
  <script type=\"text/javascript\" src=\"http://www.openbrainmap.org/cdn/js/L.TileLayer.Zoomify.js\"></script> 
  <script type=\"text/javascript\" src=\"http://www.openbrainmap.org/cdn/js/L.Control.MousePosition.js\"></script> 
  <script src=\"http://www.openbrainmap.org/cdn/js/Control.MiniMap.min.js\" type=\"text/javascript\"></script>
  <script src=\"http://www.openbrainmap.org/cdn/js/hashing.js\" type=\"text/javascript\"></script>
  <link rel=\"stylesheet\" href=\"http://www.openbrainmap.org/cdn/css/L.Control.MousePosition.css\" />
    <script>

    var scale = %1.2f;
    var bregmaX = %d;
    var bregmaY = %d;

    function onEachFeature(feature, layer) {
    var popupContent = \"<h3>Neuron ID: #\" +
        feature.id + \"<br>Data type: \" + feature.geometry.type + \"</h3><hr>\";

    if (feature.properties && feature.properties.popupContent) {
      popupContent += feature.properties.popupContent;
    }

    layer.bindPopup(popupContent);
  }
  
    var center = [180, 10];

  var original = L.tileLayer.zoomify(\'./%s/\', {
              width: %d,
              height: %d,
              tolerance: 0.8,
          });

      var minimap = L.tileLayer.zoomify(\'./%s/\', {
                  width: %d,
                  height: %d,
                  tolerance: 0.8,
              });
          
          /***  little hack to be able to have multiple popups open starts here ***/
          L.Map = L.Map.extend({
              openPopup: function(popup) {
                  //        this.closePopup();  // just comment this
                  this._popup = popup;

                  return this.addLayer(popup).fire(\'popupopen\', {
                      popup: this._popup
                  });
              }
          }); /***  end of hack ***/
            
    var map = L.map(\'map\', {
        layers: [original],
        center: center,
    zoomControl: false,
        zoom: 20,
    crs: L.CRS.Simple,
    inertia: true,
    inertiaDeceleration: 2000
    }).setView([0, 0], 1);
    map.attributionControl.setPrefix(\'Image name: | %s\');
  
  //Plugin magic goes here! Note that you cannot use the same layer object again, as that will confuse the two map controls
  var rect1 = {color: \"#64A375\", weight: 1};
  var rect2 = {color: \"#A54B6C\", weight: 1, opacity:0, fillOpacity:0};    
  var miniMap = new L.Control.MiniMap(minimap, { toggleDisplay: true, position: \"topright\", width: 200, height: 150,  aimingRectOptions : rect1, shadowRectOptions: rect2, collapsedWidth: 36, collapsedHeight: 36}).addTo(map);

  
  var southWest = map.unproject([0, %d], map.getMaxZoom());
  var northEast = map.unproject([%d,0], map.getMaxZoom());
  //map.setMaxBounds(new L.LatLngBounds(southWest, northEast));
  
  map.setView( map.unproject([%d*0.65, %d*0.65], map.getMaxZoom()) , 2);

       
    
    L.control.mousePosition().addTo(map); //.divideBy(256).floor();

    var hash = new L.Hash(map);
    
    var drawnItems = L.featureGroup().addTo(map);

    map.addControl(new L.Control.Draw({
      edit: { featureGroup: drawnItems }
    }));

    var shapeJSON;
    var markercounter = 0;

    map.on(\'draw:created\', function (e) {
    var type = e.layerType,
        layer = e.layer;

    if (type === \'marker\') {
      markercounter = markercounter +1;
      var markerPoint = map.project( layer.getLatLng(), map.getMaxZoom());
        markerstr = \"<h5>Marker ID: #\"+ markercounter+ \"</h5><hr><br><b>Stereotactic coordinates (mm):</b><br> c( \"+ ((markerPoint.x-bregmaX)/scale)/1000 +\", \" + ((markerPoint.y-bregmaY)/scale)/1000 +\" )<br><b>Pixel coordinates:</b><br> c( \" + markerPoint.x + \", \" + markerPoint.y + \" )\";
        layer.bindPopup(markerstr); 
    }

    drawnItems.addLayer(layer);

     
});
    

        var AllenOut = new L.geoJson(allenoutlines, {
            coordsToLatLng: function (latlng) {
                return (map.unproject([latlng[1], latlng[0]], map.getMaxZoom()));
            },

          style: function (feature) {
            return feature.properties && feature.properties.style;
          },

        });/* .addTo(map); */

        
        
        var geoJsonTest = new L.geoJson(neurons, {
            coordsToLatLng: function (latlng) {
                return (map.unproject([latlng[1], latlng[0]], map.getMaxZoom()));
            },
            pointToLayer: function (feature, latlng) {      
                return L.circleMarker(latlng,  {
              radius: 4,
              fillColor: feature.hexcolor,
              color: feature.linecolor,
              weight: 1,
              opacity: 1,
              fillOpacity: 0.8
            });
            },
          
          onEachFeature: onEachFeature
        });/*.addTo(map);*/
        
        var neurongroup = L.layerGroup([geoJsonTest]);
        var allengroup = L.layerGroup([AllenOut]);
        
        var baseMaps = {
            \"%s\": original,
        };
        ', scale, bregmaX, bregmaY, tiled_imagefolder, a$width, a$height, tiled_imagefolder, a$width, a$height, outputfile, a$height, a$width, a$width, a$height, fluorophore)


    footer2<-'
        var overlayMaps = {
            \"Segmented cell bodies\": neurongroup,
          \"Region boundaries\": allengroup
        };
        
        L.control.layers(baseMaps, overlayMaps, {position: \'bottomright\'}).addTo(map);

        $(document).bind(\'keydown\', \'p\', function printoutputlog() {shapeJSON = drawnItems.toGeoJSON();console.log(JSON.stringify(shapeJSON));});


        L.easyButton( \'fa-download\', function printstuff(){ var shapes = \'manual.output<-list(\\n\'; var k = 1;

    drawnItems.eachLayer(function(layer) {
        if(k > 1){shapes = shapes + \',\\n\'};
        // Note: Rectangle extends Polygon. Polygon extends Polyline.
        // Therefore, all of them are instances of Polyline
        if ((layer instanceof L.Polyline) && ! (layer instanceof L.Polygon)) {
            shapes = shapes + \" list(type = \\"PolyLine\\", coords = c(\";
            var numPoints = layer.getLatLngs().length;
            var coordinateString = \'\';
            var stereoCoordinateString = \'\';

            for (i=0; i<(numPoints-1); i++){

              var pixelpoint = map.project(layer.getLatLngs()[i], map.getMaxZoom());

              coordinateString = coordinateString + pixelpoint.x.toString() + \", \" + pixelpoint.y.toString() + \", \";
              stereoCoordinateString = stereoCoordinateString + parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString() + \", \";
            }

            var pixelpoint = map.project(layer.getLatLngs()[(numPoints-1)], map.getMaxZoom());

            coordinateString = coordinateString + pixelpoint.x.toString() + \", \" + pixelpoint.y.toString();
            stereoCoordinateString = stereoCoordinateString + parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString();

            shapes = shapes + coordinateString + \"), stereo.coords = c(\"+ stereoCoordinateString + \"))\"
        }

        if ((layer instanceof L.Polygon) && ! (layer instanceof L.Rectangle)) {
            shapes = shapes + \" list(type = \\"Polygon\\", coords = c(\";
            var numPoints = layer.getLatLngs().length;
            var coordinateString = \'\';
            var stereoCoordinateString = \'\';

            for (i=0; i<(numPoints-1); i++){

              var pixelpoint = map.project(layer.getLatLngs()[i], map.getMaxZoom());

              coordinateString = coordinateString + pixelpoint.x.toString() + \", \" + pixelpoint.y.toString() + \", \";
              stereoCoordinateString = stereoCoordinateString + parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString() + \", \";
            }

            var pixelpoint = map.project(layer.getLatLngs()[(numPoints-1)], map.getMaxZoom());

            coordinateString = coordinateString + pixelpoint.x.toString() + \", \" + pixelpoint.y.toString();
            stereoCoordinateString = stereoCoordinateString + parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString();

            shapes = shapes + coordinateString + \"), stereo.coords = c(\"+ stereoCoordinateString + \"))\"
        }

        if (layer instanceof L.Rectangle) {
            shapes = shapes + \" list(type = \\"Rectangle\\", coords = c(\";
            var numPoints = layer.getLatLngs().length;
            var coordinateString = \'\';
            var stereoCoordinateString = \'\';

            for (i=0; i<(numPoints-1); i++){

              var pixelpoint = map.project(layer.getLatLngs()[i], map.getMaxZoom());

              coordinateString = coordinateString + pixelpoint.x.toString() + \", \" + pixelpoint.y.toString() + \", \";
              stereoCoordinateString = stereoCoordinateString + parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString() + \", \";
            }

            var pixelpoint = map.project(layer.getLatLngs()[(numPoints-1)], map.getMaxZoom());

            coordinateString = coordinateString + pixelpoint.x.toString() + \", \" + pixelpoint.y.toString();
            stereoCoordinateString = stereoCoordinateString + parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString();

            shapes = shapes + coordinateString + \"), stereo.coords = c(\"+ stereoCoordinateString + \"))\"
        }


        if (layer instanceof L.Circle) {
            shapes = shapes + \" list(type = \\"Circle\\", coords = c(\";
            var pixelpoint = map.project(layer.getLatLng(), map.getMaxZoom());
            var coordinateString = pixelpoint.x.toString() + \", \" + pixelpoint.y.toString();
            var stereoCoordinateString = parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString();
            shapes = shapes + coordinateString + \"), stereo.coords = c(\"+ stereoCoordinateString + \"))\"
        }

        if (layer instanceof L.Marker) {
            shapes = shapes + \" list(type = \\"Marker\\", coords = c(\";
            var pixelpoint = map.project(layer.getLatLng(), map.getMaxZoom());
            var coordinateString = pixelpoint.x.toString() + \", \" + pixelpoint.y.toString();
            var stereoCoordinateString = parseFloat( ((pixelpoint.x-bregmaX)/scale)/1000 ).toFixed(3).toString() + \", \" + parseFloat( ((pixelpoint.y-bregmaY)/scale)/1000 ).toFixed(3).toString();
            shapes = shapes + coordinateString + \"), stereo.coords = c(\"+ stereoCoordinateString + \"))\"

            //shapes.push(coordinateString);
        }
        k++;
    });
    shapes = shapes + \"\\n)\";    

    console.log(shapes);
var blob = new Blob([shapes], {type: \'application/octet-binary\'});
    // pass a useful mime type here
var url = URL.createObjectURL(blob);
window.open(url, \'_blank\');


  }).addTo(map);




   L.easyButton( \'fa-print\', function printPage(){

//$(\'#map\').print();
  if (this.elementsToHide){
    var htmlElementsToHide = document.querySelectorAll(this.elementsToHide);  

    for (var i = 0; i < htmlElementsToHide.length; i++) {
      htmlElementsToHide[i].className = htmlElementsToHide[i].className + \' _epHidden\';
    }
  }
  window.print();

  if (this.elementsToHide){
    var htmlElementsToHide = document.querySelectorAll(this.elementsToHide);  

    for (var i = 0; i < htmlElementsToHide.length; i++) {
      htmlElementsToHide[i].className = htmlElementsToHide[i].className.replace(\' _epHidden\',\'\');
    }
  } 

}).addTo(map);
    
    </script>
</body>
</html>'
    
    setwd(paste('Web', outputfile, sep='_'))
    sink(paste(outputfile,".html",sep=''), append=FALSE)
    cat(headerP01)
    cat(headerP02)
    cat(JSONdata)
    cat(footer)
    cat(footer2)
    sink()
    setwd('../')
    return(a)
}