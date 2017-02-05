
printJSONneuron<-function(neuX, neuY, intensity, area, structureAN, ytop){
#probably should remove ytop
cat('var neurons = {\n
  
    "type": "FeatureCollection",\n
    "features": [\n')
for(i in 1:length(neuX)){
    if(structureAN[i]==0){
    line.color<-"#ff0000"
    
    cat(paste('{\n
            "geometry": {\n
                "type": "Point",\n
                "coordinates": [\n',
          neuY[i],',\n',
          neuX[i],
                '\n]\n',
            '},\n "type": "Feature",\n "properties": {\n "popupContent": "<b><u>Pixel properties:</u></b> <br> <b>Mean intensity:</b> ', round(log(intensity[i],2),3), ' bits<br> <b>Soma size:</b> ', round(area[i]), ' &mu;m<sup>2</sup> <br> </br> <b><u>Allen P56 reference atlas:</u></b> <br> <b>Region: </b>', '?', ': ', 'Undefined/unspecific', '<br>"\n},\n "hexcolor": ', paste('"',ontology$allen.color[which(ontology$id== structureAN[i])], '"', sep=''), ',\n"linecolor": ', paste('"',line.color, '"', sep=''),', \n"id": ', i, '\n },\n'))
    
  }else{
  line.color<-"#000"
  
  cat(paste('{\n
            "geometry": {\n
                "type": "Point",\n
                "coordinates": [\n',
          neuY[i],',\n',
          neuX[i],
                '\n]\n',
            '},\n "type": "Feature",\n "properties": {\n "popupContent": "<b><u>Pixel properties:</u></b> <br> <b>Mean intensity:</b> ', round(log(intensity[i],2),3), ' bits<br> <b>Soma size:</b> ', round(area[i]), ' &mu;m<sup>2</sup> <br> </br> <b><u>Allen P56 reference atlas:</u></b> <br> <b>Region: </b>', ontology$acronym[which(ontology$id== structureAN[i])], ': ',  ontology$name[which(ontology$id== structureAN[i])], '<br>"\n},\n "hexcolor": ', paste('"',ontology$allen.color[which(ontology$id== structureAN[i])], '"', sep=''), ',\n"linecolor": ', paste('"',line.color, '"', sep=''),', \n"id": ', i, '\n },\n'))
  
  }
  
  
}
cat(']\n};')
}



printJSONoutlines<-function(registration){
scaleup<-registration$transformationgrid$width/dim(registration$transformationgrid$mx)[2]
k=registration$atlas$numRegions
line.type<-rep(2,k)
line.type[which(registration$atlas$col%in%c('#cccccc', '#aaaaaa'))]<-1
for(i in 1:k){

if(i==1){
  cat('var allenoutlines = [\n')
}

if(!is.na(registration$atlas$outlines[[i]]$xrT)) { 
###############   
cat('   {\n')
cat('    "type": "Feature",\n')
cat('    "properties": {"typeofregion": "graymatter",')
cat('             "popupContent": ', paste('"<b><u>Allen P56 reference atlas:</u></b> <br> <b>Region: </b>', 'ACRONYM', ': ',  'NAME', '<br>"'), ',\n')  
cat(paste('         "id":', i ,',\n'))
cat('             "style": {\n')
cat('                 weight: 1.2,\n')
cat('                 color: "#F59B24",\n')
cat('                 opacity: 1,\n')
cat('                 fillColor: "#ffffff",\n')
if(line.type[i]==2){cat('         dashArray: "1, 5",\n')} 
cat('                 fillOpacity: 0.0\n')
cat('              }\n\n')      
cat('    },\n')
cat('    "geometry": {\n')
cat('        "type": "Polygon",\n')
cat('        "coordinates": [[\n')

for(j in 1:length(registration$atlas$outlines[[i]]$xrT)){
  if(j != length(registration$atlas$outlines[[i]]$xrT)){
    cat(paste('            [', registration$atlas$outlines[[i]]$yrT[j]*scaleup,', ', registration$atlas$outlines[[i]]$xrT[j]*scaleup, '],\n', sep=''))
  }else{
    cat(paste('            [', registration$atlas$outlines[[i]]$yrT[j]*scaleup,', ', registration$atlas$outlines[[i]]$xrT[j]*scaleup, ']\n', sep=''))
  }
}
cat('        ]]\n')
cat('    }\n')
cat('},\n')
}

###########

if(!is.na(registration$atlas$outlines[[i]]$xlT)) { 
###############   
cat('   {\n')
cat('    "type": "Feature",\n')
cat('    "properties": {"typeofregion": "graymatter",')
cat('             "popupContent": ', paste('"<b><u>Allen P56 reference atlas:</u></b> <br> <b>Region: </b>', 'ACRONYM', ': ',  'NAME', '<br>"'), ',\n')  
cat(paste('         "id":', i ,',\n'))
cat('             "style": {\n')
cat('                 weight: 1.2,\n')
cat('                 color: "#F59B24",\n')
cat('                 opacity: 1,\n')
cat('                 fillColor: "#ffffff",\n')
if(line.type[i]==2){cat('         dashArray: "1, 5",\n')} 
cat('                 fillOpacity: 0.0\n')
cat('              }\n\n')      
cat('    },\n')
cat('    "geometry": {\n')
cat('        "type": "Polygon",\n')
cat('        "coordinates": [[\n')

for(j in 1:length(registration$atlas$outlines[[i]]$xlT)){
  if(j != length(registration$atlas$outlines[[i]]$xlT)){
    cat(paste('            [', registration$atlas$outlines[[i]]$ylT[j]*scaleup,', ', registration$atlas$outlines[[i]]$xlT[j]*scaleup, '],\n', sep=''))
  }else{
    cat(paste('            [', registration$atlas$outlines[[i]]$ylT[j]*scaleup,', ', registration$atlas$outlines[[i]]$xlT[j]*scaleup, ']\n', sep=''))
  }
}
cat('        ]]\n')
cat('    }\n')
cat('},\n')
}
###########
    
    

if(i==k){
  cat('];')
} 

}

}


write.section<-function(imgURL){
  tobewritten<-sprintf("<!-- BEGIN -->  
<div class=\"element-item ctxpl str d1\" data-category=\"d1\">
   <iframe src=\'http://www.wholebrainsoftware.org/interface/brainmap.php?image=%27../sections/1/1/71/7.1-1.02_ch00_modified/%27&width=20281&height=14834\' frameborder=\"1\" style=\"border-width: 8px;   -moz-border-radius: 12px;
      -webkit-border-radius: 12px;
      border-radius: 10px;box-shadow: 0px 0px 0px 4px #6e84bd;\" width=287px height=240px></iframe>
   <h3 class=\"name\">Mercury</h3>
   <div style=\"float:left;width: 100%\">
      <div style=\"float:left;width:30%;\">
         <ul class=\"coordinates\" style=\"visibility: visible; padding: 0px;margin: 0 0;\">
            <li class=\"subject\">
               %d
            </li>
            <li class=\"AP\">5.8</li>
            <li class=\"ML\">
               0
            </li>
            <li class=\"DV\">
               0
            </li>
            <li class=\"cellcount\">
               3450
            </li>
         </ul>
      </div>
      <div style=\"float:left;width:40%;\">
         <ul class=\"coordinates\" style=\"visibility: visible; padding: 0px;margin: 0 0;\">
            <li class=\"transgenic\">
               %s
            </li>
            <li class=\"labeling\">
               %s
            </li>
            <li class=\"target\">
               %s
            </li>
         </ul>
      </div>
      <ul class=\"tools\" >
         <li class=\"expand\">
            <a href=\"%s\" class=\"opener\">Expand</a>
         </li>
      </ul>
   </div>
</div>
<!-- END -->", imgURL, subjectID, line, probe, target, imgURL)

}

slidetray.HTML<-function(){
  header01<-"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" 
\"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">
<html xmlns=\"http://www.w3.org/1999/xhtml\">
<head>  
  <meta charset=\"utf-8\">
  <meta http-equiv=\"X-UA-Compatible\" content=\"IE=edge,chrome=1\">
  <title>WholeBrain - webinterface</title>
  <meta name=\"description\" content=\"\">
  <meta name=\"author\" content=\"\">
  
  <meta http-equiv=\"Cache-Control\" content=\"no-cache\" />
  <meta http-equiv=\"Pragma\" content=\"no-cache\" />
  <meta http-equiv=\"Expires\" content=\"0\" />
  
  <script src=\"http://cdn.openbrainmap.org/0.0.1/slidetray/pace/pace.min.js\"></script>
  <link href=\"http://cdn.openbrainmap.org/0.0.1/slidetray/pace/pace.css\" rel=\"stylesheet\" />
  
  
  <link rel=\"stylesheet\" href=\"http://cdn.leafletjs.com/leaflet-0.5.1/leaflet.css\" />
      <link href='http://fonts.googleapis.com/css?family=Open+Sans:400,300,800,800italic,400italic|Adamina' rel='stylesheet' type='text/css'>

    <link href=\"http://netdna.bootstrapcdn.com/font-awesome/3.2.1/css/font-awesome.css\" rel=\"stylesheet\"> 
  <style type=\"text/css\">
        
    </style>
  
  <link rel=\"stylesheet\" href=\"http://cdn.openbrainmap.org/0.0.1/slidetray/css/normalise.css\">
  <link rel=\"stylesheet\" href=\"http://cdn.openbrainmap.org/0.0.1/slidetray/css/custom.css\">
  
    <link rel=\"stylesheet\" href=\"http://cdn.openbrainmap.org/0.0.1/slidetray/css/style.css\" media=\"screen\" type=\"text/css\" />
    <script type=\"text/javascript\" src=\"http://cdn.openbrainmap.org/0.0.1/slidetray/js/jquery-latest.js\"></script>

  
  <link rel=\"stylesheet\" href=\"http://cdn.leafletjs.com/leaflet-0.5.1/leaflet.css\" />
</head>"

header02<-"<body>
<div class=\'filterpanel\'>
  Filter menu<br>   
      <div id=\"filters\" class=\"button-group\">  
        <button class=\"button is-checked\" data-filter=\"*\">show all<br><div class=\'img-circular\' style=\"background-image: url(\'./img/ALL_parts01.png\');\"></div><p>ALL</button>
        <button class=\"button\" data-filter=\".ctxpl\"><font style=\'font-size: 12px;\'>Cortical plate</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/CTXpl_parts01.png\');\"></div><p>CTXpl</button>
        <button class=\"button\" data-filter=\".ctxsp\"><font style=\'font-size: 12px;\'>Cortical subplate</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/CTXsp_parts01.png\');\"></div><p>CTXsp</button>
        <button class=\"button\" data-filter=\".str\"><font style=\'font-size: 12px;\'>Striatum</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/STR_parts01.png\');\"></div><p>STR</button>
        <button class=\"button\" data-filter=\".pal\"><font style=\'font-size: 12px;\'>Pallidum</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/PAL_parts01.png\');\"></div><p>PAL</button>
        <button class=\"button\" data-filter=\".th\"><font style=\'font-size: 12px;\'>Thalamus</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/TH_parts01.png\');\"></div><p>TH</button>
        <button class=\"button\" data-filter=\".hy\"><font style=\'font-size: 12px;\'>Hypothalamus</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/HY_parts01.png\');\"></div><p>HY</button>
        <button class=\"button\" data-filter=\".mb\"><font style=\'font-size: 12px;\'>Midbrain</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/MB_parts01.png\');\"></div><p>MB</button>  
        <button class=\"button\" data-filter=\".hb\"><font style=\'font-size: 12px;\'>Hindbrain</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/HB_parts01.png\');\"></div><p>HB</button>
        <button class=\"button\" data-filter=\".cb\"><font style=\'font-size: 12px;\'>Cerebellum</font><br><div class=\'img-circular\' style=\"background-image: url(\'http://cdn.openbrainmap.org/0.0.1/slidetray/img/CB_parts01.png\');\"></div><p>CB</button>
    </div>
    <div style=\"margin: -10px; display: inline;\">
      <div id=\"sorts\" class=\"button-group\"  style=\"display:inline-block; margin: 10px; vertical-align: top;\">  
        <div class=\'sorters\'>
          Order:</br>
        <button class=\"button\" data-sort-by=\"name\" ascending=\"true\">image name</button>
        <button class=\"button\" data-sort-by=\"subject\" ascending=\"false\">animal</button>
        <button class=\"button\" data-sort-by=\"AP\" ascending=\"false\">AP</button>
        <button class=\"button\" data-sort-by==\"ML\" ascending=\"false\">ML</button>
        <button class=\"button\" data-sort-by=\"DV\" ascending=\"false\">DV</button>
        <button class=\"button\" data-sort-by=\"probe\" ascending=\"false\">probe</button>
        <button class=\"button\" data-sort-by=\"transgenic\" ascending=\"false\">genotype</button>
          </div>
      <!--</div> -->
                
      <!-- <div id=\"sorts\" class=\"button-group\" style=\"display:inline-block; margin: 10px; vertical-align: top;\"> -->
        <div class=\'sorters2\'>
        Cell count:</br>
        <button class=\"button\" data-sort-by=\"cellcount\" ascending=\"false\">total</button>
        <button class=\"button\" data-sort-by=\"cellregions\" ascending=\"false\">regions</button>
          <button class=\"button\" data-sort-by=\"cellcount\" ascending=\"true\">min</button> 
          <button class=\"button\" data-sort-by=\"cellcount\" ascending=\"false\">max</button> 
          </div>
      </div>"

customGenotypeMenu<-"<div id=\"filtersGeno\" class=\"button-group\" style=\"display:inline-block; margin: 10px; vertical-align: top;\">
        
        <div class=\'sorters\'>
        Project specific tags:</br>
      <ul class=\"tags\" >
        <li><a href=\"#\" data-filter=\".d1\">Input to D1 dorsal striatum</a></li>
      </ul>
      <ul class=\"tags\" >
        <li><a href=\"#\" class=\"second\" data-filter=\".chat\">Input to ChAT dorsal striatum</a></li>
      </ul>

      </div>
      </div>
      <br>
      <br>
      <br>
      
    </div>
      

      
    </div>

      
      
      <div class=\"isotope\">"




}


slidetray<-function(input, filter, folder.prefix=NULL, start.at=1, dont.run=NULL, coordinates=NULL){


    if(length(input)>1){
      #process bunch of images
      file <- as.character(input)[1]
      if(!file.exists(file))
      stop(file, "not found")

      folder<-input
      folder<-strsplit(folder, "\\/")[[1]]
      folder<-folder[-length(folder)]
      folder<-paste(folder, collapse='/')
      images<-input

    }else{
      if(is.null(folder.prefix)){
        #process single folder with images
        folder<-input
        images<-get.images(input)
      }else{
        #process folder with folder that have a prefix
        all.section.folder<-list.dirs(input, recursive=FALSE, full.names=FALSE)
        index<-c(which(substr(all.section.folder, 1,nchar(folder.prefix))%in%c(folder.prefix)))
        if(length(index)>0){
          all.section.folder<-all.section.folder[index]
        }
        #start at
        if(start.at>1){
          if(length(start.at)>1){
            all.section.folder<-all.section.folder[start.at[1]:start.at[2]]
          }else{
            all.section.folder<-all.section.folder[start.at:length(all.section.folder)]
          }
        }
        #dont run
        if(!is.null(dont.run)){
          all.section.folder<-all.section.folder[-dont.run]
        }
        images<-get.images(paste(input, all.section.folder, sep='/') )
     }
     #end if else
    }

    image.name<-basename(images)
    image.id<-formatC(1:length(images),digits=nchar(length(images)),flag="0")
    if(!is.null(coordinates)){
      if(length(coordinates)!=length(image.id)){stop("Length of coordinates is not the same as number of images")}
    }
    print("Making web maps... this can take some time.")
    if(is.missing(filter)){
      invisible( lapply(images, function(x){makewebmap(x, filter)}) )
    }else{
      invisible( lapply(images, function(x){makewebmap(x)}) )
    }


  

  #FFC.folder<-paste('FFC',basename(section.folder), sep='_')
      #FFC.folder<-paste(folder, FFC.folder, sep='/')
      #images<-get.images(FFC.folder)
      #images<-images[index]
      

}



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

makewebmap<-function(img, filter, registration=NULL, dataset=NULL, folder.name = NULL, scale = 0.64, bregmaX = 0, bregmaY = 0, fluorophore = 'Rabies-EGFP', combine = NULL, enable.drawing=TRUE, write.tiles=TRUE, verbose=FALSE){
    file <- as.character(img)[1]
    if(!file.exists(file))
    stop(file, "not found")

    if(length(img)>1){

    }else{
      ## expand path
      file <- path.expand(file)
      outputfile<-basename(file)
      #outputfile<-sub("^([^.]*).*", "\\1", outputfile) #this cannot handle filenames with punctiation
      outputfile<-strsplit(outputfile, "\\.")[[1]]
      outputfile<-paste(outputfile[-length(outputfile)], collapse='.')
      folder.name<-paste('Web', outputfile, sep='_')
      create.output.directory(folder.name, verbose=verbose)

    }

    setwd(folder.name)
    lapply(outputfile, function(x){create.output.directory(paste('Tiles', x, sep='_'), verbose=verbose)} )
    setwd('../')
    
    #if(is.null(overlap)){overlap<-(-999)}
    verbose<-as.numeric(verbose)


    tiled_imagefolder<-paste('Tiles', outputfile, sep='_')
    
    if(!missing(filter)){
      alpha<-filter$Max
      beta<-filter$Min
      if(is.null(filter$Min)){
        beta<-0
      }
    }else{
      if(verbose){
        cat('No filter, creating 8-bit range from quantiles. \n')
      }
      alpha<-0
      beta<-0
    }
    if(write.tiles){
      a <- .Call("createWeb", file, alpha, beta, verbose, outputfile)
    }
    if(!is.null(registration)){
      bregmaX<-round(stereotactic.coordinates(0, 0, registration, inverse = TRUE)$x)
      bregmaY<-round(stereotactic.coordinates(0, 0, registration, inverse = TRUE)$y) 
    }


    headerP01<-'<!DOCTYPE html>
<html>
<head>
    <title>'


headerP02<-"</title>

    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no\">

    <link href=\"http://maxcdn.bootstrapcdn.com/font-awesome/4.1.0/css/font-awesome.min.css\" rel=\"stylesheet\">
    <link rel=\"stylesheet\" href=\"http://cdn.openbrainmap.org/0.0.1/dist/leaflet.css\" />
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

JSONfilename<-paste('./', outputfile,'.js', sep='')
JSONdata<-paste('<script src=\"', JSONfilename,'\" type=\"text/javascript\"></script>',sep='') 

footer<-sprintf('<div id=\"map\"></div>

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
        feature.id + \"</h3><hr>\";

    if (feature.properties && feature.properties.popupContent) {
      popupContent += feature.properties.popupContent;
    }

    layer.bindPopup(popupContent);
  }
  
    var center = [180, 10];

  var original = L.tileLayer.zoomify(\'./%s/\', {
              width:  imageProperties.width,
              height: imageProperties.height,
              tolerance: 0.8,
          });

      var minimap = L.tileLayer.zoomify(\'./%s/\', {
                  width:  imageProperties.width,
                  height: imageProperties.height,
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

  
  var southWest = map.unproject([0,  imageProperties.height], map.getMaxZoom());
  var northEast = map.unproject([imageProperties.width,0], map.getMaxZoom());
  //map.setMaxBounds(new L.LatLngBounds(southWest, northEast));
  
  map.setView( map.unproject([imageProperties.width*0.65,  imageProperties.height*0.65], map.getMaxZoom()) , 2);

       
    
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

     
});', scale, bregmaX, bregmaY, tiled_imagefolder, tiled_imagefolder, outputfile)

  // registration=NULL, dataset=NULL
    if(is.null(registration)&is.null(dataset)){
    footer1<-sprintf('
        
        var baseMaps = {
            \"%s\": original,
        };

         var overlayMaps = {
          \"Manually defined\": drawnItems
        };
        
        ', fluorophore)
    }else{
      if(is.null(registration)){
    footer1<-sprintf('
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
        });
        
        var neurongroup = L.layerGroup([geoJsonTest]);
        
        var baseMaps = {
            \"%s\": original,
        };

         var overlayMaps = {
          \"Segmented cell bodies\": neurongroup,
          \"Manually defined\": drawnItems
        };
        ', fluorophore)

        }else{
          if(is.null(dataset)){
            footer1<-sprintf('var AllenOut = new L.geoJson(allenoutlines, {
            coordsToLatLng: function (latlng) {
                return (map.unproject([latlng[1], latlng[0]], map.getMaxZoom()));
            },

          style: function (feature) {
            return feature.properties && feature.properties.style;
          },

        });

      
        
        var neurongroup = L.layerGroup([geoJsonTest]);
        
        var baseMaps = {
            \"%s\": original,
        };

         var overlayMaps = {
          \"Region boundaries\": allengroup,
          \"Manually defined\": drawnItems
        };
        
        ', fluorophore)
          }else{
            footer1<-sprintf('var AllenOut = new L.geoJson(allenoutlines, {
            coordsToLatLng: function (latlng) {
                return (map.unproject([latlng[1], latlng[0]], map.getMaxZoom()));
            },

          style: function (feature) {
            return feature.properties && feature.properties.style;
          },

        });

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
        });
        
        var neurongroup = L.layerGroup([geoJsonTest]);
        var allengroup = L.layerGroup([AllenOut]);
        
        var baseMaps = {
            \"%s\": original,
        };

         var overlayMaps = {
          \"Segmented cell bodies\": neurongroup,
          \"Region boundaries\": allengroup,
          \"Manually defined\": drawnItems
        };
        
        ', fluorophore)
          }

        }
    } 

    footer2<-'
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
    sink(paste(outputfile,".js",sep=''), append=FALSE)
cat('var imageProperties = {
    "width": ',a$width,',
    "height": ', a$height,'};\n')
if( (!is.null(dataset))&(!is.null(registration)) ){
  printJSONneuron(dataset$x, dataset$y, dataset$intensity, dataset$area*scale, dataset$id, registration$transformationgrid$height)
  printJSONoutlines(registration)
}else{
  if((is.null(dataset))&(!is.null(registration))){
      printJSONoutlines(registration)
  }
  cat('var neurons = {};')
  cat('var allenoutlines = {};')
}
    sink()
    setwd('../')

    return(paste(getwd(), paste('Web', outputfile, sep='_'), sep='' ))
}