#' Segment image using multiple thresholds
#'
#' Segment brain section.
#' @param input input a character vector consisting of the full path name to 16-bit raw tif image files.
#' @param numthresh numbe rof binary thresholds to use (default is 8).
#' @examples
#' #folder where image tiles are stored
#'images<-get.images('/Volumes/microscope/animal001/slide001/section001')
#' #stitch images
#' segment(system.file('sample_tiles/rabiesEGFP.tif', package='wholebrain')) 

segment<-function(input, numthresh=8, downsample=0.25, filter=NULL, post=NULL, pre=NULL, get.contour=FALSE, display=TRUE){
  inputfile<-character()
  for(i in 1:length(input)){
    inputfile <- as.character(input[i])
    ## check for existence
    if(!file.exists(inputfile))
      stop(inputfile, ", file not found")
    inputfile <- path.expand(inputfile)
  #  files<-append(files, file)
  }

  file <- system.file('double_slider.png', package='wholebrain')
  fileslider <- system.file('slider.png', package='wholebrain')
  filebackground <- system.file('GUI_background.png', package='wholebrain')
  resizeP = as.integer(downsample*100)
  #check for filter
  if(is.null(filter)){
    areaMin<-(-999)
    areaMax<-(-999)
    threshMin<-(-999)
    threshMax<-(-999)
    eccent<-(-999)
    renderMax<-(-999)
    renderMin<-(-999)
    bThresh<-(-999)
    resizeB<-(-999)
    gaussBlur<-(-999)
  }else{
    areaMin<-filter$alim[1]
    areaMax<-filter$alim[2]
    threshMin<-filter$threshold.range[1]
    threshMax<-filter$threshold.range[2]
    eccent<-filter$eccentricity[1]
    renderMax<-filter$Max
    renderMin<-filter$Min
    bThresh<-filter$brain.threshold
    resizeB<-filter$resize
    resizeP<-as.integer(filter$downsample*100)
    gaussBlur<-filter$blur
    downsample<-filter$downsample
  }
  #check for 
  if(is.null(pre)){
    prefiltering<-FALSE
    transformation<-0
    iterations<-0   
  }else{
    prefiltering<-TRUE
    transformation<-pre$trans
    iterations<-pre$iter   
  }

  #check for post-filtering
  if(is.null(post)){
    postfiltering<-TRUE
    transformation<-0
    iterations<-0   
  }else{
    postfiltering<-TRUE
    transformation<-pre$trans
    iterations<-pre$iter   
  }  

  a<-.Call("GUI",inputfile,numthresh, resizeP,file,fileslider,filebackground, display, areaMin, areaMax, threshMin, threshMax, eccent, renderMin, renderMax, bThresh, resizeB, gaussBlur, as.integer(get.contour))
  a$x<-(1/ downsample)*a$x
  a$y<-(1/ downsample)*a$y
  a$contour.x<-(1/ downsample)*a$contour.x
  a$contour.y<-(1/ downsample)*a$contour.y
  a$soma.area <-(1/ downsample)*a$soma.area
  outputlist<-list(filter=list(alim= a$alim, threshold.range = a$threshold.range, eccentricity = a$eccentricity,  Max = a$Max, Min = a$Mina, brain.threshold=a$brain.threshold, resize=a$resize, blur=a$blur, downsample=a$downsample), soma = list(x =a$x, y=a$y, intensity = a$intensity, area = a$soma.area, contour.x= a$contour.x, contour.y=a$contour.y, contour.ID=a$contour.ID))
  if(is.null(outputlist$filter$Min)){
    outputlist$filter$Min<-0
  }
  return(outputlist)
}


#' Get contour
#'
#' Get the contours of a binary segmented image
#' @param input input a character vector consisting of the full path name to 16-bit raw tif image files.
#' @param coordinate matching coordinates in the main plane to be registered to (in millimeters).
#' @param plane the main plane to register too: "coronal", "sagital", "".
#' @param brain.threshold a integer value, which determien sthe segmentation of the brain slice.
#' @param verbose boolean value. If true diagnostic output is written to the R console. Deafult is true.
#' @examples
#' #path to image
#' image<-'/Volumes/microscope/animal001/slide001/section001.tif'
#' #register the image
#' registration(image, AP=1.05, brain.threshold=220)
nuclear.segment<-function(input, thresh=0.1, kernel=3, iter=20, threshold = 'otsu', invert=FALSE, get.largest.object = TRUE, num.nested.objects = 2, blur=0, watershed=FALSE, resize=0.25, display=TRUE, verbose=TRUE, savefilename=FALSE){
    file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)
    if(is.character(threshold)){threshold<-0}
    resizeP = resize
    saveoutput<-0
     if(savefilename==TRUE){
      savefilename<-paste(getwd(),'/miniature_',basename(file), sep='')
      saveoutput<-1
     }

    a<-.Call("nuclearSegment", file, as.numeric(thresh), as.integer(kernel), as.integer(iter), as.integer(threshold), as.integer(invert), as.integer(get.largest.object), as.integer(num.nested.objects), as.integer(display), resizeP, blur, as.integer(verbose), savefilename, saveoutput)
    if(saveoutput==1){
    print(savefilename)
    }
    return(a)
}

imshow<-function(input, auto.range = FALSE, quantile = c(0.001,0.999), resize=0.25){
  inputfile<-character()
  for(i in 1:length(input)){
    inputfile <- as.character(input[i])
    ## check for existence
    if(!file.exists(inputfile))
      stop(inputfile, ", file not found")
    inputfile <- path.expand(inputfile)
  #  files<-append(files, file)
  }

  file <- system.file('double_slider.png', package='wholebrain')
  fileslider <- system.file('slider.png', package='wholebrain')
  filebackground <- system.file('GUI_background.png', package='wholebrain')
  resizeP = as.integer(resize*100)
  a<-.Call("imageshow",inputfile, as.integer(auto.range), quantile[1], quantile[2], resizeP,file,fileslider,filebackground)
  if(auto.range){b<-a$quantile.values
    names(b)<-paste(round(a$quantile*100,2), '%', sep='')
    return(b)
  }
}