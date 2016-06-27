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

segment<-function(input, numthresh=8, resize=0.25, filter=NULL, display=TRUE){
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
  a<-.Call("GUI",inputfile,numthresh, resizeP,file,fileslider,filebackground, display)
  a$x<-(1/ resize)*a$x
  a$y<-(1/ resize)*a$y
  a$soma.area <-(1/ resize)*a$soma.area
  outputlist<-list(filter=list(alim= a$alim, threshold.range = a$threshold.range, eccentricity = a$eccentricity,  Max = a$Max, Min = a$Mina), soma = list(x =a$x, y=a$y, intensity = a$intensity, area = a$soma.area))
  return(outputlist)
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