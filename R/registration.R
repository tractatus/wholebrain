data(EPSatlas, envir=environment())
data(atlasIndex, envir=environment())

#' Get atlas binary image
#'
#' 
#' @param coordinate matching coordinates in the main plane to be registered to (in millimeters).
#' @param directory either an absolute path or the default TEMPORARY which will save the tif output to .
#' @param plate if true coordinate indicates plate number out of 132 plates.
#' @param width the pixel width to save as TIFF, default is 456px, but 11300 is high res.
#' @param verbose boolean value. If true diagnostic output is written to the R console. Deafult is true.
#' @examples
#' #path to image
#' image<-'/Volumes/microscope/animal001/slide001/section001.tif'
#' #register the image
#' registration(image, AP=1.05, brain.threshold=220)

get.atlas.image<-function(coordinate, directory='TEMPORARY', width=456, plate=FALSE, close.image=TRUE, save.image=TRUE){

if(plate){
  k<-coordinate
}else{
  k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))
plate<-atlasIndex$plate.id[which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))]
}
  

  
ventricles <-c(which(EPSatlas$plate.info[[k]]$style=='#aaaaaa')) #117 is optic chiasm och   which(EPSatlas$plate.info[[k]]$structure_id== 117)
if(atlasIndex$plate.id[k] %in%c(100960309, 100960312, 100960316, 100960320)){
  fibertracts<-c(which(EPSatlas$plate.info[[k]]$style=='#cccccc'))
  fibertracts <-fibertracts[which.min(sapply(fibertracts, function(x){min(EPSatlas$plates[[k]][[x]]@paths$path@y )}))]
  ventricles<-append(ventricles, fibertracts) 
}


xmin<-min(EPSatlas$plates[[k]][[1]]@paths$path@x)-97440/2

quartz(width= 11.300 , height= 7.900)
par(mar=c(0,0,0,0), bg='black', xaxs='i', yaxs='i')
plot(EPSatlas$plates[[k]][[1]]@paths$path@x, EPSatlas$plates[[k]][[1]]@paths$path@y, col=0, xlim=c(0,97440), ylim=c(0, 68234.56), axes=F)
polygon(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin, EPSatlas$plates[[k]][[1]]@paths$path@y, col='white', border='white' )
polygon(-(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin - 97440/2)+97440/2 , EPSatlas$plates[[k]][[1]]@paths$path@y, col='white', border='white')
if(coordinate <1.3){
for(i in ventricles){
polygon(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, col='black', border='black', lwd=4)
polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin- 97440/2)+97440/2, EPSatlas$plates[[k]][[i]]@paths$path@y, col='black', border='black', lwd=4)
}
}
if(directory == 'TEMPORARY'){
  full.filename<-paste(tempdir(),'/',plate,'.tif', sep='')
}else{
  full.filename<-paste(directory,'/',plate,'.tif', sep='')
} 
cat(paste('Plate',k,'out of',132,'\n') )
if(save.image){
dev.copy(tiff, full.filename, width = width, height = 7900*(width/11300), units = "px", res = 150)
dev.off()
cat(paste('Saved as:', full.filename,'\n') )
}
if(close.image){dev.off()}
return(full.filename)
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
get.contour<-function(input, threshold = 'otsu', get.largest.object = TRUE, watershed=FALSE, verbose=TRUE){
    file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)
    .Call("segment", file, as.integer(threshold), as.integer(verbose))
}

#' Register
#'
#' An image sensor of a camera is composed of a two dimensional array of light sensitive detectors or pixels. The sesnor array is #'mechanically quite stable with the pixels retaining a rigidly fixed geometric relationship. Each pixel within the array, however, #'has its own unique light sensitivity characteristics. As these characteristics affect camera performance, they must be removed #'through calibration. The process by which a camera is calibrated is known as "Flat Fielding" or "Shading Correction".
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

registration<- function(input, coordinate=5.04, plane="coronal", brain.threshold = 200, scale.factor=0.28, plateimage = FALSE, verbose=TRUE){

    file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)

    threshold<-brain.threshold

    if(plane!='coronal'){
      stop(plane, " plane is not supported yet.")
    }

    if(plateimage!=FALSE){
      plateimage <- as.character(plateimage)
      ## check for existence
      if(!file.exists(plateimage))
        stop(plateimage, ", file not found")
      plateimage <- path.expand(plateimage)

    }

  
  .Call("ThinPlateRegistration", file, as.integer(threshold), as.integer(verbose))
} 

testregistration<-function(input, brain.threshold = 200, verbose=TRUE){
	file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)

	threshold<-brain.threshold
  .Call("ThinPlateRegistration", file, as.integer(threshold), as.integer(verbose))
} 