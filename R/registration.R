data(EPSatlas, envir=environment())
data(SAGITTALatlas, envir=environment())

data(atlasIndex, envir=environment())
#data(atlasOntology, envir=environment())
data(ontology, envir=environment())


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

get.atlas.image<-function(coordinate, directory='TEMPORARY', width=456, plane='coronal', plate=FALSE, right.hemisphere=NULL, close.image=TRUE, save.image=TRUE){
plate.width<-1
#get cutting plane
if(plane=="sagittal"){
 EPSatlas<-SAGITTALatlas
 atlasIndex<-atlasIndex[atlasIndex$plane=="sagittal", ]
 plate.width<-1.159292
}else{
   atlasIndex<-atlasIndex[atlasIndex$plane=="coronal", ]
}

#get closest coordinate in millimeters
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


xmin<-(plane!='sagittal')*(min(EPSatlas$plates[[k]][[1]]@paths$path@x)-97440/2)

quartz(width= plate.width*11.300 , height= 7.900)
par(mar=c(0,0,0,0), bg='black', xaxs='i', yaxs='i')
plot(EPSatlas$plates[[k]][[1]]@paths$path@x, EPSatlas$plates[[k]][[1]]@paths$path@y, col=0, xlim=c(0, plate.width*97440), ylim=c(0, 68234.56), axes=F)
if(is.null(right.hemisphere)&(plane!='sagittal') ){
  polygon(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin, EPSatlas$plates[[k]][[1]]@paths$path@y, col='white', border='white' )
  polygon(-(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin - (plate.width*97440)/2)+(plate.width*97440)/2 , EPSatlas$plates[[k]][[1]]@paths$path@y, col='white', border='white')
}else{
  #check if right or left hemisphere
  if( ( (plane!='sagittal')&(right.hemisphere==TRUE) ) | ( (plane=='sagittal')&(right.hemisphere==FALSE) ) ){
      polygon(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin, EPSatlas$plates[[k]][[1]]@paths$path@y, col='white', border='white' )
    }else{
      polygon(-(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin - (plate.width*97440)/2)+(plate.width*97440)/2 , EPSatlas$plates[[k]][[1]]@paths$path@y, col='white', border='white')
    }
}

if(plane!='sagittal'){

if( coordinate <1.3 ){
for(i in ventricles){
  polygon(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, col='black', border='black', lwd=4)
  polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin- (plate.width*97440)/2)+(plate.width*97440)/2, EPSatlas$plates[[k]][[i]]@paths$path@y, col='black', border='black', lwd=4)
}
}else{
  for(i in ventricles){
    polygon(EPSatlas$plates[[k]][[i]]@paths$path@x, EPSatlas$plates[[k]][[i]]@paths$path@y, col='black', border='black', lwd=4)
  }
}
}

if(directory == 'TEMPORARY'){
  full.filename<-paste(tempdir(),'/',plate,'.tif', sep='')
}else{
  full.filename<-paste(directory,'/',plate,'.tif', sep='')
} 

cat(paste('Plate',k,'out of',132,'\n') )
if(save.image){
dev.copy(tiff, full.filename, width = plate.width*width, height = 7900*(width/11300), units = "px", res = 150)
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
get.contour<-function(input, threshold = 'otsu', invert=FALSE, get.largest.object = TRUE, num.nested.objects = 2, blur=0, watershed=FALSE, resize=0.25, display=TRUE, verbose=TRUE, savefilename=FALSE){
    file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)
    if(is.character(threshold)){threshold<-0}
    resizeP = resize
    saveoutput<-0
     if(savefilename==TRUE){
      savefilename<-paste(tempdir(),'/miniature_',basename(file), sep='')
      saveoutput<-1
     }

    a<-.Call("getContour", file, as.integer(threshold), as.integer(invert), as.integer(get.largest.object), as.integer(num.nested.objects), as.integer(display), resizeP, blur, as.integer(verbose), savefilename, saveoutput)
    if(saveoutput==1){
    print(savefilename)
    }
    return(a)
}


atlasToimage<-function(atlas, affineX, affineY, downsmpX, downsmpY, translateX, translateY){
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


  scale.factor<-(1/97440)*4

  imgWidth<-scale.factor*456
  imgHeight<-scale.factor*320



  affineX<- matfile$affineX*(97440/11300)*(11300/(456*9.75))
  affineY<- matfile$affineY*(97440/11300)*(11300/(456*9.75))
  
  xmax<-(97440/11300)*downsmpX*11300/(456*9.75)
  ymax<-(97440/11300)*downsmpY*11300/(456*9.75)
  
  center<-min(EPSatlas$plates[[k]][[1]]@paths$path@x)-97440/2
  #for loop here
  numPaths<-atlas@summary@numPaths
  outlines<-list()
  for(i in 1:numPaths){
    xr<-(atlas[[i]]@paths$path@x-center+affineX)*scale.factor
    yr<-(atlas[[i]]@paths$path@y)*scale.factor+affineY
    xl<-(-atlas[[i]]@paths$path@x+97440+ center )*scale.factor+affineX
    yl<-(atlas[[i]]@paths$path@y)*scale.factor+affineY
  
    index<-cbind(downsmpY-as.integer(round(yr)), as.integer(round(xr)))
    xrT<-(xr+translateX[index] )
    yrT<-(yr-translateY[index] )

    index<-cbind(downsmpY-as.integer(round(yl)), as.integer(round(xl)))
    xlT<-(xl+translateX[index] )
    ylT<-(yl-translateY[index] )
    
    outlines[[i]]<-list(xrT  = xrT, yrT = yrT, xlT = xlT, ylT = ylT)
  }
  return(outlines)  
}


get.forward.warp<-function(registration){
  trans<-list()
  trans$mx<-registration$transformationgrid$mx
  trans$my<-registration$transformationgrid$my

  mx<-matrix(rep(NA,prod(dim(trans$mx))), nrow=dim(trans$mx)[1])
  my<-matrix(rep(NA,prod(dim(trans$mx))), nrow=dim(trans$mx)[1])

  for(i in 1:nrow(mx)){
    for(j in 1:ncol(mx)){
      if((trans$mx[i,j]>=0)&(trans$mx[i,j]<dim(trans$mx)[2])&(trans$my[i,j]>=0)&(trans$my[i,j]<dim(trans$mx)[1])){
        mx[trans$my[i,j]+1,trans$mx[i,j]+1]<-j
        my[trans$my[i,j]+1,trans$mx[i,j]+1]<-i
      }
    }
  }

  index<-which(is.na(my), arr.ind=TRUE)
  for(i in 1:nrow(index)){
    if(index[i,1]>=dim(mx)[1]){
      top<-0
    }else{
      top<-1
    }
    if(index[i,1]<=1){
      bottom<-0
    }else{
      bottom<-1
    }
    if(index[i,2]>=dim(mx)[2]){
      right<-0
    }else{
      right<-1
    }
    if(index[i,2]<=1){
      left<-0
    }else{
      left<-1
    }
    mx[index[i,1],index[i,2]]<-mean(c(mx[index[i,1]-bottom,index[i,2]], mx[index[i,1]+top,index[i,2]], mx[index[i,1],index[i,2]+right], mx[index[i,1],index[i,2]-left] ),na.rm=TRUE)
    my[index[i,1],index[i,2]]<-mean(c(my[index[i,1]-bottom,index[i,2]], my[index[i,1]+top,index[i,2]], my[index[i,1],index[i,2]+right], my[index[i,1],index[i,2]-left] ),na.rm=TRUE)
  }

  registration$transformationgrid<-list(mx=registration$transformationgrid$mx, my=registration$transformationgrid$my, mxF=mx, myF=my, width=registration$transformationgrid$width  ,height=registration$transformationgrid$height)

  return(registration)
}


get.forward.warpRCPP<-function(registration){
  trans<-list()
  trans$mx<-registration$transformationgrid$mx
  trans$my<-registration$transformationgrid$my

  mx<-matrix(rep(-999,prod(dim(trans$mx))), nrow=dim(trans$mx)[1])
  my<-matrix(rep(-999,prod(dim(trans$mx))), nrow=dim(trans$mx)[1])
  a<-.Call("forwardWarp", mx, my, trans$mx, trans$my)
  mx<-a$mx
  my<-a$my

  #index<-which(is.na(my), arr.ind=TRUE)
  #for(i in 1:nrow(index)){
  #  if(index[i,1]>=dim(mx)[1]){
  #    top<-0
  #  }else{
  #    top<-1
  #  }
  #  if(index[i,1]<=1){
  #    bottom<-0
  #  }else{
  #    bottom<-1
  #  }
  #  if(index[i,2]>=dim(mx)[2]){
  #    right<-0
  #  }else{
  #    right<-1
  #  }
  #  if(index[i,2]<=1){
  #    left<-0
  #  }else{
   #   left<-1
  #  }
   # mx[index[i,1],index[i,2]]<-mean(c(mx[index[i,1]-bottom,index[i,2]], mx[index[i,1]+top,index[i,2]], mx[index[i,1],index[i,2]+right], mx[index[i,1],index[i,2]-left] ),na.rm=TRUE)
   # my[index[i,1],index[i,2]]<-mean(c(my[index[i,1]-bottom,index[i,2]], my[index[i,1]+top,index[i,2]], my[index[i,1],index[i,2]+right], my[index[i,1],index[i,2]-left] ),na.rm=TRUE)
  #}


  registration$transformationgrid<-list(mx=registration$transformationgrid$mx, my=registration$transformationgrid$my, mxF=mx, myF=my, width=registration$transformationgrid$width  ,height=registration$transformationgrid$height)

  return(registration)
}


custom.registration<-function(input, atlas){

  transformationgrid<-.Call("ThinPlateRegistration", file, targetP.x, targetP.y, referenceP.x, referenceP.y, resizeP, MaxDisp, MinDisp, outputfile)
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

registration<- function(input, coordinate=NULL, plane="coronal", right.hemisphere=NULL, brain.threshold = 200, blurring=c(4,15), pixel.resolution=0.64, resize=(1/8)/4, correspondance=NULL, resolutionLevel=c(4,2), num.nested.objects=0, display=TRUE, plateimage = FALSE, forward.warp=FALSE, filter=NULL, output.folder='../', batch.mode=FALSE, verbose=TRUE){
    if(.Platform$OS.type=="windows") {
      
      batch.mode=TRUE
  }

    if(is.null(coordinate)){
      if(is.null(correspondance)){
        coordinate<-0
      }else{
        coordinate<-correspondance$coordinate
      }
    }

    file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)


    outputfile<-basename(file)
   # outputfile<-sub("^([^.]*).*", "\\1", outputfile)
    outputfile<-strsplit(outputfile, "\\.")[[1]]
    outputfile<-paste(outputfile[-length(outputfile)], collapse='.')

    defaultwd<-getwd()

    if(output.folder=='./'){
      parentpath<-dirname(input)[1]
    }
    if(output.folder=='../'){
      parentpath<-dirname(dirname(input))[1]
    }

    outfolder<-paste('output', outputfile, sep='_')
    setwd(parentpath)
    create.output.directory(outfolder, verbose=verbose)
    setwd(defaultwd)

    outputfile<-paste(parentpath, outfolder,paste('Registration', outputfile, sep='_'), sep='/')


    plate.width<-1
    SAGITTAL<-TRUE#right.hemisphere change this if sgaittal not working
  #get cutting plane
  if(plane=="sagittal"){
    EPSatlas<-SAGITTALatlas
     atlasIndex<-atlasIndex[atlasIndex$plane=="sagittal", ]
    plate.width<-1.159292
    SAGITTAL<-!SAGITTAL
  }else{
   atlasIndex<-atlasIndex[atlasIndex$plane=="coronal", ]
}



    if(plateimage!=FALSE){
      plateimage <- as.character(plateimage)
      ## check for existence
      if(!file.exists(plateimage))
        stop(plateimage, ", file not found")
      plateimage <- path.expand(plateimage)

    }

      if(is.null(filter)){
        MaxDisp<-0
        MinDisp<-0
      }else{
        MaxDisp<-filter$Max
        MinDisp<-filter$Min
        brain.threshold<-filter$brain.threshold
        resize<-filter$resize
        blurring[1]<-filter$blur
      }
      if(is.null(MinDisp)){
        MinDisp<-0
      }


    if(is.null(correspondance)){
    #scale.factor<-pixel.resolution*0.390625

    #get correspondance points for input image
    contourInput<-get.contour(file, thresh=brain.threshold, resize=resize, blur=blurring[1], num.nested.objects=num.nested.objects, display=FALSE)
    contoursI<-as.numeric(names(sort(tapply(contourInput$x,contourInput$contour.ID, min))))

    #get correspondance points for atlas
    filename<-get.atlas.image(coordinate, right.hemisphere=right.hemisphere, plane=plane)
    contourAtlas<-get.contour(filename, resize=1, blur=blurring[2], num.nested.objects=num.nested.objects, display=FALSE)

    contours<-as.numeric(names(sort(tapply(contourAtlas$x,contourAtlas$contour.ID, min))))
    cor.pointsInput<-list()
    cor.pointsAtlas<-list()
    for(i in 1:length(contours) ){
      if(contours[i]==0){
        resLevel<-resolutionLevel[1]
      }else{
        resLevel<-resolutionLevel[2]
      }
      cor.pointsAtlas[[i]]<-automatic.correspondences(cbind(contourAtlas$x[which(contourAtlas$contour.ID==contours[i])],contourAtlas$y[which(contourAtlas$contour.ID==contours[i])]), resLevel, plot=FALSE)
      cor.pointsInput[[i]]<-automatic.correspondences(cbind(contourInput$x[which(contourInput$contour.ID==contoursI[i])],contourInput$y[which(contourInput$contour.ID==contoursI[i])]), resLevel, plot=FALSE)
    }

    cor.points<-list(atlas=cor.pointsAtlas, input=cor.pointsInput)

    centroidAtlas<-cor.points$atlas[[1]]$q[1,]
    offsetAtlas<-centroidAtlas-456/2

    targetP.x<-numeric()
    referenceP.x<-numeric()
    targetP.y<-numeric()
    referenceP.y<-numeric()
    shape<-numeric()

    #scale.reference<-mean(c((max(contourInput$x)-min(contourInput$x))/(max(contourAtlas$x)-min(contourAtlas$x)),(max(contourInput$y)-min(contourInput$y))/(max(contourAtlas$y)-min(contourAtlas$y)) ) )

    #centroidNorm<-scale.reference*centroidAtlas-cor.points$input[[1]]$q[1,]
    centroidNorm<-centroidAtlas-cor.points$input[[1]]$q[1,]
    for(i in 1:length(contours) ){
      #referenceP.x<-append(referenceP.x, (scale.reference*cor.points$atlas[[i]]$p[,1]-centroidNorm[1]))
      #referenceP.y<-append(referenceP.y, (scale.reference*cor.points$atlas[[i]]$p[,2]-centroidNorm[2]))
      referenceP.x<-append(referenceP.x, 4*(cor.points$atlas[[i]]$p[,1]-centroidNorm[1]))
      referenceP.y<-append(referenceP.y, 4*(cor.points$atlas[[i]]$p[,2]-centroidNorm[2]))

      targetP.x<-append(targetP.x, 4*cor.points$input[[i]]$p[,1])
      targetP.y<-append(targetP.y, 4*cor.points$input[[i]]$p[,2])
      shape<-append(shape, rep(i, length(cor.points$input[[i]]$p[,2])))
    }

    }else{
      correspondance$correspondance<-na.omit(correspondance$correspondance)
      targetP.x<-correspondance$correspondance[,1]
      targetP.y<-correspondance$correspondance[,2]
      referenceP.x<-correspondance$correspondance[,3]
      referenceP.y<-correspondance$correspondance[,4]
      shape<-correspondance$correspondance[,5]
      centroidNorm<-correspondance$centroidNorm
      centroidAtlas<-correspondance$centroidAtlas
    }

    resizeP<-resize*4
    #resizeP<-resize*4
    #resizeP<-(25/resizeP)
    #resizeP<-(pixel.resolution/resizeP)
   #(1/8)
   #/4)*9.75

  transformationgrid<-.Call("ThinPlateRegistration", file, targetP.x, targetP.y, referenceP.x, referenceP.y, resizeP, MaxDisp, MinDisp, outputfile)


  #transform outlines
  k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))
  xmin<-(plane!='sagittal')*(min(EPSatlas$plates[[k]][[1]]@paths$path@x)-97440/2)


  numPaths<-EPSatlas$plates[[k]]@summary@numPaths
  scale.factor<-456/97440
  style<-EPSatlas$plate.info[[k]]$style

  outlines<-list()
  for(i in 1:numPaths){
    if(is.null(right.hemisphere)){
      xr<-4*((EPSatlas$plates[[k]][[i]]@paths$path@x-xmin)*scale.factor-centroidNorm[1])
      yr<-4*(((-EPSatlas$plates[[k]][[i]]@paths$path@y)*scale.factor+320)-centroidNorm[2])
      xl<-4*((-(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin - (plate.width*97440)/2)+(plate.width*97440)/2)*scale.factor-centroidNorm[1])
      yl<-4*(((-EPSatlas$plates[[k]][[i]]@paths$path@y)*scale.factor+320)-centroidNorm[2])


      index<-cbind(as.integer(round(yr)), as.integer(round(xr)))
      xrT<-(xr+(transformationgrid$mx[index]-xr) )
      yrT<-(yr+(transformationgrid$my[index]-yr) )

      index<-cbind(as.integer(round(yl)), as.integer(round(xl)))
      xlT<-(xl+(transformationgrid$mx[index]-xl) )
      ylT<-(yl+(transformationgrid$my[index]-yl) )
    }else{
      if(right.hemisphere==SAGITTAL){
        xr<-4*((EPSatlas$plates[[k]][[i]]@paths$path@x-xmin)*scale.factor-centroidNorm[1])
        yr<-4*(((-EPSatlas$plates[[k]][[i]]@paths$path@y)*scale.factor+320)-centroidNorm[2])
        xl<-NA
        yl<-NA


        index<-cbind(as.integer(round(yr)), as.integer(round(xr)))
        xrT<-(xr+(transformationgrid$mx[index]-xr) )
        yrT<-(yr+(transformationgrid$my[index]-yr) )

        
        xlT<-NA
        ylT<-NA

      }else{
        xr<-NA
        yr<-NA
        xl<-4*((-(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin - (plate.width*97440)/2)+(plate.width*97440)/2)*scale.factor-centroidNorm[1])
        yl<-4*(((-EPSatlas$plates[[k]][[i]]@paths$path@y)*scale.factor+320)-centroidNorm[2])

        xrT<-NA
        yrT<-NA

        index<-cbind(as.integer(round(yl)), as.integer(round(xl)))
        xlT<-(xl+(transformationgrid$mx[index]-xl) )
        ylT<-(yl+(transformationgrid$my[index]-yl) )
      }
    }
  
    outlines[[i]]<-list(xr  = xr, yr = yr, xl = xl, yl = yl, xrT = xrT, yrT= yrT, xlT = xlT, ylT = ylT)
  }



    if(display){
    par(yaxs='i', xaxs='i')
        img<-paste(outputfile,'_undistorted.png', sep='')
        img <- readPNG(img)

        img = as.raster(img[,])

        #img <- apply(img, 2, rev)
        if(batch.mode){img <- apply(img, 2, rev)}

        par(xaxs='r', yaxs='r')
        plot(c(0, dim(img)[2]*2),c(0, dim(img)[1]), axes=F, asp=1, col=0, xlab='', ylab='', ylim=c(dim(img)[1],0))
        polygon(c(-5, dim(img)[2]+5, dim(img)[2]+5, -5), c(-5, -5, dim(img)[1]+5, dim(img)[1]+5) , col='orange')
        polygon(c(-5+ dim(img)[2], 2*dim(img)[2]+5, 2*dim(img)[2]+5, -5+ dim(img)[2]), c(-5, -5, dim(img)[1]+5, dim(img)[1]+5) , col='purple')
        rasterImage(img,0,0, dim(img)[2], dim(img)[1])
        rasterImage(img,dim(img)[2],0, 2*dim(img)[2], dim(img)[1])
        abline(v=dim(img)[2], lwd=2, col='white')
             
        if(is.null(right.hemisphere)){     
          lapply(1:numPaths, function(x){polygon(outlines[[x]]$xr,outlines[[x]]$yr, border='orange' )})
          lapply(1:numPaths, function(x){polygon(outlines[[x]]$xl,outlines[[x]]$yl, border='orange' )})

          lapply(1:numPaths, function(x){polygon(dim(img)[2]+outlines[[x]]$xrT,outlines[[x]]$yrT, border='purple' )})
          lapply(1:numPaths, function(x){polygon(dim(img)[2]+outlines[[x]]$xlT,outlines[[x]]$ylT, border='purple' )})
        }else{
          if(right.hemisphere==SAGITTAL){
              lapply(1:numPaths, function(x){polygon(outlines[[x]]$xr,outlines[[x]]$yr, border='orange' )})
              lapply(1:numPaths, function(x){polygon(dim(img)[2]+outlines[[x]]$xrT,outlines[[x]]$yrT, border='purple' )})
            }else{
              lapply(1:numPaths, function(x){polygon(outlines[[x]]$xl,outlines[[x]]$yl, border='orange' )})
              lapply(1:numPaths, function(x){polygon(dim(img)[2]+outlines[[x]]$xlT,outlines[[x]]$ylT, border='purple' )})
            }
        }

        #lapply(unique(shape), function(x){polygon(referenceP.x[which(shape==x)], referenceP.y[which(shape==x)], col=rgb(1,1,0.1,0.3), border=rgb(1,1,0.1))})
        #lapply(unique(shape), function(x){polygon(targetP.x[which(shape==x)], targetP.y[which(shape==x)], col=rgb(0.1,1,1,0.3), border=rgb(0.1,1,1))})
        # lines(c(targetP.x[x],referenceP.x[x]) , c(targetP.y[x],referenceP.y[x]), col=x ); 
        lapply(1:length(targetP.x), function(x){points(c(targetP.x[x]+dim(img)[2],referenceP.x[x]) , c(targetP.y[x],referenceP.y[x]), pch=c(19), col='black', cex=1.8); text(c(targetP.x[x]+dim(img)[2],referenceP.x[x]) , c(targetP.y[x],referenceP.y[x]), label=x, col='white', cex=0.7) } )
        #legend('topright', c('Target section', 'Reference atlas', 'Overlap'), pch=c(21,23,22), col='black', bg='white', horiz=TRUE,pt.bg=c(rgb(0.1,1,1,0.3), rgb(1,1,0.1,0.3), rgb(0.72,1,0.8)))
    }

    returnlist<-list(atlas=list(outlines=outlines, numRegions=numPaths, col=style), transformationgrid=transformationgrid, correspondance=data.frame(targetP.x, targetP.y, referenceP.x, referenceP.y, shape), centroidAtlas=centroidAtlas, centroidNorm = centroidNorm, coordinate=coordinate, resize= resize, outputfile=outputfile, plane=plane )
  if(forward.warp){
  returnlist<-get.forward.warp(returnlist)
  }
  return(returnlist)
  #return(data.frame(r.x=referenceP.x, r.y=referenceP.y, t.x=targetP.x, t.y=targetP.y))
}

add.corrpoints<-function(registration, n.points=NULL){
  print("Begin with Target image (right) then add same point in reference atlas (left)")
  height.n.width<-dim(registration$transformationgrid$mx)
  if(is.null(n.points)){
    p<-locator(, type='p', pch=c(21), bg=c('lightblue'))
  }else{
    p<-locator(n.points*2, type='p', pch=c(21), bg=c('lightblue'))
  }
  new.p<-cbind(rbind(cbind(p$x[seq(1,length(p$x),by=2)], p$y[seq(1,length(p$y),by=2)])),
        rbind(cbind(p$x[seq(2,length(p$x),by=2)], p$y[seq(2,length(p$y),by=2)]))
  )
  new.corr<-data.frame(cbind(new.p,rep(4, nrow(new.p))))
  names(new.corr)<-names(registration$correspondance)
  new.corr[,1]<-(new.corr[,1]-height.n.width[2])
  registration$correspondance<-rbind(registration$correspondance, new.corr )
  return(registration)
}

remove.corrpoints<-function(registration, n.points){
  registration$correspondance<-registration$correspondance[-n.points,]
  return(registration)
}

change.corrpoints<-function(registration, n.points, target.only=TRUE){
  height.n.width<-dim(registration$transformationgrid$mx)
  if(target.only){
    print(paste("Change where you want to position points:", n.points, 'in the target image only'))
    p<-locator(length(n.points), type='p', pch=c(21), bg=c('lightblue'))
    new.p<-cbind(p$x,p$y,registration$correspondance[n.points,3:5])
    new.corr<-data.frame(new.p)
    new.corr[,1]<-(new.corr[,1]-height.n.width[2])
    names(new.corr)<-names(registration$correspondance)
    registration$correspondance[n.points,]<-new.corr
  }else{
    print("Begin with Target image (right) then add same point in reference atlas (left)")
    p<-locator(length(n.points)*2, type='p', pch=c(21), bg=c('lightblue'))
    new.p<-cbind(rbind(cbind(p$x[seq(1,length(p$x),by=2)], p$y[seq(1,length(p$y),by=2)])),
        rbind(cbind(p$x[seq(2,length(p$x),by=2)], p$y[seq(2,length(p$y),by=2)]))
    )
    new.corr<-data.frame(cbind(new.p,registration$correspondance[n.points,5] ) )
    names(new.corr)<-names(registration$correspondance)
    new.corr[,1]<-(new.corr[,1]-height.n.width[2])
    registration$correspondance[n.points,]<-new.corr
  }
  
  return(registration)
}

shift.corrpoints<-function(registration, steps=1, clockwise=TRUE, num.points=32){
  clockwise<-1:2+(!clockwise)*2
  registration$correspondance[1:num.points,clockwise]<-registration$correspondance[c((2+steps):num.points,1:(1+steps) ), clockwise]
  return(registration)
}


stereotactic.coordinates <-function(x,y,registration,inverse=FALSE){
  scale.factor<-mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height, registration$transformationgrid$width) )
   #get cutting plane
  if('plane'%in%names(registration)){
    if(registration$plane=="sagittal"){
      EPSatlas<-SAGITTALatlas
      atlasIndex<-atlasIndex[atlasIndex$plane=="sagittal", ]
      plate.width<-1.159292
      #get plate in atlas
      coordinate<-registration$coordinate
      k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))


      width<-(528-214)*2 #this will create abug because of hemisphere.

      atlas.fit.scale<-diff(range(registration$atlas$outlines[[1]]$yl/scale.factor))
      atlas.y.range<-range(EPSatlas$plates[[k]][[1]]@paths$path@y)
      atlas.scale<-diff(atlas.y.range)
    
      distance.to.atlas.center<-(atlas.y.range-68234.56/2)

      xpos<-min(-(EPSatlas$plates[[k]][[1]]@paths$path@x- (plate.width*97440)/2)+(plate.width*97440)/2)/(plate.width*97440)* plate.width*456 - scale.factor*min(registration$atlas$outlines[[1]]$xl/scale.factor)/4

      ypos<- -(range(registration$atlas$outlines[[1]]$yl/scale.factor) - (atlas.fit.scale/atlas.scale)*distance.to.atlas.center )[1]
      ypos<-scale.factor*ypos/4

    }else{
      #coronal
      width<-456
      xpos<-registration$centroidNorm[1]
      ypos<-registration$centroidNorm[2]
    }
  }

   if(inverse){
        x<-x*1000/25+width/2
        y<-(-y*1000/25)
    
        x<-4*(x-xpos)/scale.factor
        y<-4*(y-ypos)/scale.factor
        return(list(x=x, y=y))
      
      }else{
        x<-scale.factor*(x)/4
        y<-scale.factor*(y)/4
        x<-x+xpos
        y<-y+ypos
        x<-(x-width/2)*25/1000
        y<-(-y/1000*25)
        #xref<-(4*(456/2-registration$centroidNorm[1]))/scale.factor
        #yref<-(4*(320/2-registration$centroidNorm[2]))/scale.factor
  
        #xref<-(registration$centroidAtlas[1]/xref)
        #yref<-(registration$centroidAtlas[2]/yref)
        return(list(x=x, y=y))
        #return(list(x=(x*xref-456/2)*25/1000, y=-((y*yref)*25)/1000))
     }

}

plot.registration<-function(registration, main=NULL, border=rgb(154,73,109,maxColorValue=255), draw.trans.grid=FALSE, batch.mode=FALSE){

  scale.factor<-mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height,registration$transformationgrid$width) )
  
  xMax<-max(c(registration$transformationgrid$mx,registration$transformationgrid$mxF),na.rm=TRUE)*(1/scale.factor)
    xMin<-min(c(registration$transformationgrid$mx,registration$transformationgrid$mxF),na.rm=TRUE)*(1/scale.factor)
  yMax<-max(c(registration$transformationgrid$my,registration$transformationgrid$myF),na.rm=TRUE)*(1/scale.factor)
    yMin<-min(c(registration$transformationgrid$my,registration$transformationgrid$myF),na.rm=TRUE)*(1/scale.factor)

if(is.null(main)){
  main.title<-basename(registration$outputfile)
}
  plot(c(xMin, xMax), c(yMin, yMax), ylim=c(yMax,yMin), xlim=c(xMin, xMax), asp=1, axes=F, xlab='', ylab='', col=0, main=main.title, font.main=1)
  polygon(c(0,rep(registration$transformationgrid$width,2),0),c(0, 0,rep(registration$transformationgrid$height,2)))

  img<-paste(registration$outputfile,'_undistorted.png', sep='')
        img <- readPNG(img)

        img = as.raster(img[,])

        if(batch.mode){img <- apply(img, 2, rev)}


        rasterImage(img,0,0, registration$transformationgrid$width, registration$transformationgrid$height)


  lapply(1:registration$atlas$numRegions, function(x){polygon(registration$atlas$outlines[[x]]$xrT/scale.factor, registration$atlas$outlines[[x]]$yrT/scale.factor, border=border ) })
        lapply(1:numPaths, function(x){polygon(registration$atlas$outlines[[x]]$xlT/scale.factor, registration$atlas$outlines[[x]]$ylT/scale.factor, border=border )})
        
if(draw.trans.grid){
lapply(seq(1,hei,by=100), function(x){lines(registration$transformationgrid$mx[x,]/scale.factor,registration$transformationgrid$my[x,]/scale.factor, col='lightblue')})
lines(registration$transformationgrid$mx[hei,]/scale.factor,registration$transformationgrid$my[hei,]/scale.factor, col='lightblue')
lapply(seq(1, wid,by=100), function(x){lines(registration$transformationgrid$mx[,x]/scale.factor,registration$transformationgrid$my[,x]/scale.factor, col='lightblue')})
lines(registration$transformationgrid$mx[,wid]/scale.factor,registration$transformationgrid$my[, wid]/scale.factor, col='lightblue')
}

}


inspect.registration<-function(registration,segmentation,soma=TRUE, forward.warps=FALSE, batch.mode=FALSE, cex=0.5, draw.trans.grid=TRUE, width= 12.280488, height=  6.134146){
  if(.Platform$OS.type=="windows") {
      batch.mode=TRUE
  }
  quartz(width= width, height=  height)
par(yaxs='i',xaxs='i', mfrow=c(1,2), mar=c(4,4,1,1))

if(is.null(registration$transformationgrid$myF)){
  registration.nm <-deparse(substitute(registration))

  cat(paste('Forward warps has not been not computed, have to compute that first!\n If you want to avoid this from happening either turn forward.warp==TRUE in registration() or sink the output of this object into the registration output, i.e.:\n', paste(registration.nm,'<-plot.registration(', registration.nm ,',forward.warp=TRUE)', sep='') ))
  if(forward.warps){
    cat('\nCOMPUTING FORWARD WARPS... this might take some time')
    registration<-get.forward.warpRCPP(registration)
  }
}
cat('\nGetting anatomical assignment for segmented objects')
dataset<-get.cell.ids(registration, segmentation, forward.warp=forward.warps)

  scale.factor<-mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height,registration$transformationgrid$width) )
  
  xMax<-max(c(registration$transformationgrid$mx,registration$transformationgrid$mxF),na.rm=TRUE)*(1/scale.factor)
    xMin<-min(c(registration$transformationgrid$mx,registration$transformationgrid$mxF),na.rm=TRUE)*(1/scale.factor)
  yMax<-max(c(registration$transformationgrid$my,registration$transformationgrid$myF),na.rm=TRUE)*(1/scale.factor)
    yMin<-min(c(registration$transformationgrid$my,registration$transformationgrid$myF),na.rm=TRUE)*(1/scale.factor)

    
plot(c(xMin, xMax), c(yMin, yMax), ylim=c(yMax,yMin), xlim=c(xMin, xMax), asp=1, axes=F, xlab='', ylab='', col=0, main=paste('Bregma:',registration$coordinate,'mm'),font.main = 1
 )
polygon(c(0,rep(registration$transformationgrid$width,2),0),c(0, 0,rep(registration$transformationgrid$height,2)))
 numPaths<-registration$atlas$numRegions
 outlines<-registration$atlas$outlines
 
mtext('Dorso-ventral (mm)',side=2,line=1.5)
mtext('Medio-lateral (mm)',side=1,line=-1.5)

 lapply(1:numPaths, function(x){polygon(outlines[[x]]$xr/scale.factor,outlines[[x]]$yr/scale.factor, border='black', col=as.character(registration$atlas$col[x]) )})
        lapply(1:numPaths, function(x){polygon(outlines[[x]]$xl/scale.factor,outlines[[x]]$yl/scale.factor, border='black', col=as.character(registration$atlas$col[x]) )})

hei<-dim(registration$transformationgrid$mx)[1]
wid<-dim(registration$transformationgrid$mx)[2]

if(draw.trans.grid){
lapply(seq(1, hei,by=100), function(x){lines(registration$transformationgrid$mxF[x,]/scale.factor,registration$transformationgrid$myF[x,]/scale.factor, col='lightblue')})
lines(registration$transformationgrid$mxF[hei,]/scale.factor,registration$transformationgrid$myF[hei,]/scale.factor, col='lightblue')
lapply(seq(1, wid,by=100), function(x){lines(registration$transformationgrid$mxF[,x]/scale.factor,registration$transformationgrid$myF[,x]/scale.factor, col='lightblue')})
lines(registration$transformationgrid$mxF[, wid]/scale.factor,registration$transformationgrid$myF[, wid]/scale.factor, col='lightblue')
}
 
 index<-round(scale.factor*cbind(dataset$y, dataset$x))
 somaX<-registration$transformationgrid$mxF[index]/scale.factor
 somaY<-registration$transformationgrid$myF[index]/scale.factor
 #for(x in id){
  #index<-round(scale.factor*cbind(rois[[x]]$coords[,2],rois[[x]]$coords[,1]))
  #  rois2[[x]]$coords[,1]<-registration$transformationgrid$mxF[index]/scale.factor
 #   rois2[[x]]$coords[,2]<-registration$transformationgrid$myF[index]/scale.factor

#} 
#lapply(id, function(x){polygon(rois2[[x]]$coords, col=rgb(100,163,117,120,maxColorValue=255), border=gray(0.95));text(apply(rois2[[x]]$coords,2,mean)[1],apply(rois2[[x]]$coords,2,mean)[2],x, cex=0.7, col='white')})
circle.color<-rep('', nrow(dataset))
circle.color[dataset$id>0]<-'black'
circle.color[dataset$id==0]<-'red'
if(soma){
  points(somaX,somaY,pch=21,bg= dataset$color, col= circle.color, cex=cex)
}
axis(1, at=stereotactic.coordinates(seq(-4,4,by=0.1),NA,registration, inverse=TRUE)$x, line=-4, labels=FALSE, tck=-0.01, col.ticks='lightblue')
axis(1, at=stereotactic.coordinates(seq(-4,4,by=0.5),NA,registration, inverse=TRUE)$x, line=-4, labels=FALSE, tck=-0.02, col.ticks='coral')
axis(1, at=stereotactic.coordinates(c(-4:4),c(0:-6),registration, inverse=TRUE)$x, line=-4, labels=c(-4:4))

axis(2, at=stereotactic.coordinates(NA,seq(0,-6,by=-0.1),registration, inverse=TRUE)$y, line=-0.5, labels=FALSE, tck=-0.01, las=1, col.ticks='lightblue')
axis(2, at=stereotactic.coordinates(NA,seq(0,-6,by=-0.5),registration, inverse=TRUE)$y, line=-0.5, labels=FALSE, tck=-0.02, las=1, col.ticks='coral')
axis(2, at=stereotactic.coordinates(NA,c(0:-6),registration, inverse=TRUE)$y, line=-0.5, labels=c(0:-6), las=1)

plot(c(xMin, xMax), c(yMin, yMax), ylim=c(yMax,yMin), xlim=c(xMin, xMax), asp=1, axes=F, xlab='', ylab='', col=0, main=basename(registration$outputfile), font.main=1)
polygon(c(0,rep(registration$transformationgrid$width,2),0),c(0, 0,rep(registration$transformationgrid$height,2)))

img<-paste(registration$outputfile,'_undistorted.png', sep='')
        img <- readPNG(img)

        img = as.raster(img[,])

        if(batch.mode){img <- apply(img, 2, rev)}


        rasterImage(img,0,0, registration$transformationgrid$width, registration$transformationgrid$height)


lapply(1:numPaths, function(x){polygon(outlines[[x]]$xrT/scale.factor,outlines[[x]]$yrT/scale.factor, border=rgb(154,73,109,maxColorValue=255))})
        lapply(1:numPaths, function(x){polygon(outlines[[x]]$xlT/scale.factor,outlines[[x]]$ylT/scale.factor, border=rgb(154,73,109,maxColorValue=255) )})
        
if(draw.trans.grid){
lapply(seq(1,hei,by=100), function(x){lines(registration$transformationgrid$mx[x,]/scale.factor,registration$transformationgrid$my[x,]/scale.factor, col='lightblue')})
lines(registration$transformationgrid$mx[hei,]/scale.factor,registration$transformationgrid$my[hei,]/scale.factor, col='lightblue')
lapply(seq(1, wid,by=100), function(x){lines(registration$transformationgrid$mx[,x]/scale.factor,registration$transformationgrid$my[,x]/scale.factor, col='lightblue')})
lines(registration$transformationgrid$mx[,wid]/scale.factor,registration$transformationgrid$my[, wid]/scale.factor, col='lightblue')
}
if(soma){
 points(dataset$x, dataset$y, pch=21, bg= dataset$color, col= circle.color, cex=cex)
}
#lapply(id, function(x){polygon(rois[[x]]$coords, col=rgb(100,163,117,120,maxColorValue=255));text(apply(rois[[x]]$coords,2,mean)[1],apply(rois[[x]]$coords,2,mean)[2],x, cex=0.7, col='white')})
return(dataset)
}

get.cell.ids<-function(registration, segmentation, forward.warp=FALSE){
  if('plane'%in%names(registration)){
    if(registration$plane=="sagittal"){
      EPSatlas<-SAGITTALatlas
      atlasIndex<-atlasIndex[atlasIndex$plane=="sagittal", ]
    }else{
      atlasIndex<-atlasIndex[atlasIndex$plane=="coronal", ]
    } 
  }else{
      atlasIndex<-atlasIndex[atlasIndex$plane=="coronal", ]
    } 

  coordinate<-registration$coordinate

  k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))

 plate.info<-EPSatlas$plate.info[[k]] 

  #scale up transformation grid
  scale.factor<-mean(c(dim(registration$transformationgrid$mx)[1]/registration$transformationgrid$height,
   dim(registration$transformationgrid$mx)[2]/registration$transformationgrid$width) )

  outlines<-registration$atlas$outlines
  
  namelist<-as.numeric(as.vector(plate.info$structure_id))
  colorlist<-plate.info$style
  neuronregion<-rep(0,length(segmentation$soma$x))
  neuroncolor<-rep('#000000',length(segmentation$soma$x))
  right.hemisphere<-rep(NA, length(segmentation$soma$x))

  for(i in 1:length(outlines)){
    temp<-point.in.polygon(segmentation$soma$x, segmentation$soma$y, c(outlines[[i]]$xlT/scale.factor),c(outlines[[i]]$ylT/scale.factor) )
    neuronregion[which(temp==1)]<-namelist[i]
    neuroncolor[which(temp==1)]<-as.character(colorlist[i])
    right.hemisphere[which(temp==1)]<-FALSE
  
    temp<-point.in.polygon(segmentation$soma$x, segmentation$soma$y, c(outlines[[i]]$xrT/scale.factor),c(outlines[[i]]$yrT/scale.factor) )
    neuronregion[which(temp==1)]<-namelist[i]
    neuroncolor[which(temp==1)]<-as.character(colorlist[i])
    right.hemisphere[which(temp==1)]<-TRUE
  }
  
  segmentation$soma$id<- neuronregion
  segmentation$soma$color<- neuroncolor
  segmentation$soma$right.hemisphere<-right.hemisphere
  #the first four list items are x,y, area and intensity and have equal length. rest are contours with greater length.
  dataset<-data.frame(segmentation$soma[c(1:4,8:10)])  
  if(forward.warp==TRUE){
    if(!(length(registration$transformationgrid$mxF)>0) ){
      registration<-get.forward.warpRCPP(registration)
    }
    index<-round(scale.factor*cbind(dataset$y, dataset$x))
    #ensure that 0 indexes are 1
    index[index==0]<-1
    #check if point is outside image
    if( length( which(index[,1]>dim(registration$transformationgrid$mxF)[1]) ) ){
        index[which(index[,1]>dim(registration$transformationgrid$mxF)[1]) ,1]<-dim(registration$transformationgrid$mxF)[1]
    }
    if( length( which(index[,2]>dim(registration$transformationgrid$mxF)[2]) ) ){
      index[which(index[,2]>dim(registration$transformationgrid$mxF)[2]),2]<-dim(registration$transformationgrid$mxF)[2]
    }

    somaX<-registration$transformationgrid$mxF[index]/scale.factor
    somaY<-registration$transformationgrid$myF[index]/scale.factor
    tomecoord<-stereotactic.coordinates(somaX,somaY,registration, inverse=FALSE)
    dataset$ML<-tomecoord$x
    dataset$DV<-tomecoord$y
  }
  dataset$acronym<-rep(NA, length(dataset$id))
  class(dataset$acronym)<-'character'
  dataset$acronym[dataset$id>0]<-as.character(acronym.from.id(dataset$id[dataset$id>0]))
  dataset$name<-rep(NA, length(dataset$id))
  class(dataset$name)<-'character'
  dataset$name[dataset$id>0]<-as.character(name.from.id(dataset$id[dataset$id>0]))
  imagename<-substr(basename(registration$outputfile), nchar('Registration_')+1 ,nchar(basename(registration $outputfile)))
  dataset$image<-rep(imagename, nrow(dataset))
 dataset<-data.frame(animal=rep(NA, nrow(dataset)), AP=rep(registration$coordinate, nrow(dataset)), dataset)
  dataset$color<-as.character(dataset$color)
  return(dataset)
}

get.parts.of<-function(coordinate=0.5, acronym='Isocortex'){
  k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))
  plate.info<-EPSatlas$plate.info[[k]]
  unique.regions<-acronym.from.id(plate.info$structure_id)
  roi<-get.sub.structure(acronym)
  roi<-roi[which(roi%in%unique.regions)]  
  return(roi) 
}  

get.region<-function(acronym, registration){
  coordinate<-registration$coordinate
  k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))
  plate.info<-EPSatlas$plate.info[[k]]  
  get.outline<-function(acronym){
    id<-id.from.acronym(acronym)
    index<-which(plate.info$structure_id==id)
    
    if(length(index)==0){
      while(length(index)==0){
        id<-ontology$id[which(ontology$parent%in%id)]
        index<-which(plate.info$structure_id%in%id)
      }
      
    
      region<-registration$atlas$outlines[index]
      region<-lapply(1:length(region), function(x){
        data.frame(xT=c(region[[x]]$xlT, region[[x]]$xrT), yT=c(region[[x]]$ylT, region[[x]]$yrT), x=c(region[[x]]$xl, region[[x]]$xr), y=c(region[[x]]$yl, region[[x]]$yr), right.hemisphere = c(rep(FALSE, length(region[[x]]$xl) ), rep(TRUE, length(region[[x]]$xr) ) ), name=acronym.from.id(id[x]) )
      } )
      region <- do.call("rbind", region)
      
      
    }else{
      region<-registration$atlas$outlines[index]
      region<-data.frame(xT=c(region[[1]]$xlT, region[[1]]$xrT), yT=c(region[[1]]$ylT, region[[1]]$yrT), x=c(region[[1]]$xl, region[[1]]$xr), y=c(region[[1]]$yl, region[[1]]$yr), right.hemisphere = c(rep(FALSE, length(region[[1]]$xl) ), rep(TRUE, length(region[[1]]$xr) ) ) )
    region<-data.frame(region, name=rep(acronym, nrow(region)))
    }
    return(region)
  }
  if(length(acronym)==1){
    region <-get.outline(acronym)
    
  }else{
    region<-lapply(acronym, get.outline)
    region <- do.call("rbind", region)
  }
  if(length(region)>1)
  
  scale.factor<-mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height,registration$transformationgrid$width) )
  
  region[,1:4]<-region[,1:4]*(1/scale.factor)
  
  return(region)
}

update.regi<-function(regi, manual.output, filter=NULL, display=FALSE){
  if(!is.list(manual.output)){
    ee <- new.env()
    sys.source(manual.output, ee)
    manual.output<-ee$manual.output
  }

  #read in manual.output markers
  coords<-data.frame(x=numeric(),y=numeric())
  for(i in seq_along(manual.output)){
    if(manual.output[[i]]$type=='Marker'){
      coords<-rbind(coords, manual.output[[i]]$coords)
      names(coords)<-c('x', 'y')
    }
  }

  #indexed on odd versus even
  odd<-seq(1, nrow(coords), by=2)
  even<-seq(2, nrow(coords), by=2)
  
  #scale fatcor becuase we often downsample and resize the original image for spped and to meet the atlas.
  scale.factor<-dim(regi$transformationgrid$mx)[2]/regi$transformationgrid$width

  #target is the image that the atlas should transform to
  targetP.x<-scale.factor*coords$x[odd]
  targetP.y<-scale.factor*coords$y[odd]
  #reference is what you want it to transform to
  referenceP.x<-scale.factor*coords$x[even]
  referenceP.y<-scale.factor*coords$y[even]  
  
  if(is.null(regi$transformationgrid$myF)){

  cat(paste('Forward warps has not been not computed, have to compute that first!', sep='') )
    cat('\nCOMPUTING FORWARD WARPS... this might take some time')
    regi <-get.forward.warpRCPP(regi)
    
  }
  
  referenceP.y <-regi$transformationgrid$myF[round(cbind(referenceP.y, referenceP.x))]  
  referenceP.x <-regi$transformationgrid$mxF[round(cbind(referenceP.y, referenceP.x))]  
  
  new.corrpoints<-data.frame(cbind(targetP.x, targetP.y, referenceP.x, referenceP.y))
  new.corrpoints$shape<-1
  
  regi$correspondance <-rbind(regi$correspondance , new.corrpoints )
  
  
#get tiff filename from registration file
  file <- paste0(gsub("Registration_", "", gsub(c("output_"), "", regi$outputfile)), '.tif')
  cat("Computing registration with updated corrpoints...")
  regi<-registration(file, coordinate= regi$coordinate, filter= filter, correspondance=regi, display=display)
  
  return(regi)


} 


get.cellcounts<-function(formula = acronym ~ right.hemisphere + animal, roi=NULL, dataset, exclude=TRUE){
  dataset$acronym<-as.character(dataset$acronym)
  

  for(i in roi){
    id<-id.from.acronym(i)
    index<-which(dataset$id==id)
    
    if(length(index)==0){
      while(length(index)==0){
        id<-ontology$id[which(ontology$parent%in%id)]
        index<-which(dataset$id%in%id)
      }
      }
      
      dataset$acronym[index]<-i
    } 
    
    if(exclude){
       dataset<-dataset[which(dataset$acronym%in%roi), ]
    } 
      
      mf <- model.frame(formula=formula, data= dataset)
      
    
     return( table(mf) )
 }



get.XML.info<-function(XML.info.file){
data <- xmlParse(XML.info.file) 

plateinfo <-c(id=numeric(),parent_id=numeric() , order=numeric()  , structure_id=numeric(), style=character())

root <- xmlRoot(data)

objects<-xmlChildren(xmlChildren(root)[[1]])


paths<-xmlChildren( objects[[1]] )
for(i in 1:length(paths)){
plateinfo<-rbind(plateinfo, xmlAttrs( paths[[i]] ))
}

plateinfo <-cbind(plateinfo[,c(1,3,4)], substr(plateinfo[,6], nchar(plateinfo[1,6])-6, nchar(plateinfo[1,6]) ) )
plateinfo <-data.frame(plateinfo)
names(plateinfo)<-c('id', 'order', 'structure_id', 'style')
return(plateinfo)
}


create.registration.plate<-function(eps.file="test.eps", svg.file="test.svg"){
eps<-PostScriptTrace(eps.file) # 100960324_2.eps
mainAtlas <- readPicture(paste0(eps.file , '.xml'))


plateInfo<-get.XML.info(svg.file)

return(list(plates=mainAtlas, plate.info= plateInfo))

}




custom.get.cell.ids<-function(registration, segmentation){

  #scale up transformation grid
  scale.factor<-mean(c(dim(registration$transformationgrid$mx)[1]/registration$transformationgrid$height,
   dim(registration$transformationgrid$mx)[2]/registration$transformationgrid$width) )

  outlines<-registration$atlas$outlines
  
  namelist<-as.numeric(as.vector(plate.info$structure_id))
  colorlist<-plate.info$style
  neuronregion<-rep(0,length(segmentation$soma$x))
  neuroncolor<-rep('#000000',length(segmentation$soma$x))
  right.hemisphere<-rep(NA, length(segmentation$soma$x))

  for(i in 1:length(outlines)){
    temp<-point.in.polygon(segmentation$soma$x, segmentation$soma$y, c(outlines[[i]]$xlT/scale.factor),c(outlines[[i]]$ylT/scale.factor) )
    neuronregion[which(temp==1)]<-namelist[i]
    neuroncolor[which(temp==1)]<-as.character(colorlist[i])
    right.hemisphere[which(temp==1)]<-FALSE
  
    temp<-point.in.polygon(segmentation$soma$x, segmentation$soma$y, c(outlines[[i]]$xrT/scale.factor),c(outlines[[i]]$yrT/scale.factor) )
    neuronregion[which(temp==1)]<-namelist[i]
    neuroncolor[which(temp==1)]<-as.character(colorlist[i])
    right.hemisphere[which(temp==1)]<-TRUE
  }
  
  segmentation$soma$id<- neuronregion
  segmentation$soma$color<- neuroncolor
  segmentation$soma$right.hemisphere<-right.hemisphere
  dataset<-data.frame(segmentation$soma)
  if(forward.warp==TRUE){
    if(!(length(registration$transformationgrid$mxF)>0) ){
      registration<-get.forward.warpRCPP(registration)
    }
    index<-round(scale.factor*cbind(dataset$y, dataset$x))
    #ensure that 0 indexes are 1
    index[index==0]<-1
    #check if point is outside image
    if( length( which(index[,1]>dim(registration$transformationgrid$mxF)[1]) ) ){
        index[which(index[,1]>dim(registration$transformationgrid$mxF)[1]) ,1]<-dim(registration$transformationgrid$mxF)[1]
    }
    if( length( which(index[,2]>dim(registration$transformationgrid$mxF)[2]) ) ){
      index[which(index[,2]>dim(registration$transformationgrid$mxF)[2]),2]<-dim(registration$transformationgrid$mxF)[2]
    }

    somaX<-registration$transformationgrid$mxF[index]/scale.factor
    somaY<-registration$transformationgrid$myF[index]/scale.factor
    tomecoord<-stereotactic.coordinates(somaX,somaY,registration, inverse=FALSE)
    dataset$ML<-tomecoord$x
    dataset$DV<-tomecoord$y
  }
  dataset$acronym<-rep(NA, length(dataset$id))
  class(dataset$acronym)<-'character'
  dataset$acronym[dataset$id>0]<-as.character(acronym.from.id(dataset$id[dataset$id>0]))
  dataset$name<-rep(NA, length(dataset$id))
  class(dataset$name)<-'character'
  dataset$name[dataset$id>0]<-as.character(name.from.id(dataset$id[dataset$id>0]))
  imagename<-substr(basename(registration$outputfile), nchar('Registration_')+1 ,nchar(basename(registration $outputfile)))
  dataset$image<-rep(imagename, nrow(dataset))
 dataset<-data.frame(animal=rep(NA, nrow(dataset)), AP=rep(registration$coordinate, nrow(dataset)), dataset)
  dataset$color<-as.character(dataset$color)
  return(dataset)
}




custom.registration<-function(input='/Volumes/Seagate Backup Plus Drive/transgenic/copy_transgenic_Lhx6_EGFP.tif' , EPSfile='~/Documents/Annotation2014_141_0254-01.eps', output.folder='../', directory='TEMPORARY', close.image=TRUE, save.image=TRUE, filter=seg$filter, width=456, resize=(1/8)/4, correspondance=NULL, resolutionLevel=c(4,2)){


file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)

eps<-PostScriptTrace(EPSfile) 
atlas <- readPicture(paste0(getwd(),'/', tools::file_path_sans_ext(basename(EPSfile)), '.eps.xml')
) #100960324_2
plate.width<-1 #1.159292

quartz(width= plate.width*11.300 , height= 7.900)
par(mar=c(0,0,0,0), bg='black', xaxs='i', yaxs='i')
plot(atlas[[1]]@paths$path@x, atlas[[1]]@paths$path@y, col=0, xlim=c(0, plate.width*97440), ylim=c(0, 68234.56), axes=F)
polygon(atlas[[1]]@paths$path@x, atlas[[1]]@paths$path@y, col='white', border='white' )

plate<-basename(tools::file_path_sans_ext(EPSfile))

if(directory == 'TEMPORARY'){
  full.filename<-paste(tempdir(),'/',plate,'.tif', sep='')
}else{
  full.filename<-paste(directory,'/',plate,'.tif', sep='')
} 

if(save.image){
dev.copy(tiff, full.filename, width = plate.width*width, height = 7900*(width/11300), units = "px", res = 150)
dev.off()
cat(paste('Saved as:', full.filename,'\n') )
}

if(close.image){dev.off()}


if(is.null(correspondance)){
#get the contours
contourInput<-get.contour(file, thresh=filter$brain.threshold, resize=filter$resize, blur=filter$blur, num.nested.objects=0, display=FALSE)  
contourAtlas<-get.contour(full.filename, resize=1, blur=15, num.nested.objects=0, display=FALSE)

contoursI<-as.numeric(names(sort(tapply(contourInput$x,contourInput$contour.ID, min))))
contours<-as.numeric(names(sort(tapply(contourAtlas$x,contourAtlas$contour.ID, min))))

#get the resolution level
    cor.pointsInput<-list()
    cor.pointsAtlas<-list()
    for(i in 1:length(contours) ){
      if(contours[i]==0){
        resLevel<-resolutionLevel[1]
      }else{
        resLevel<-resolutionLevel[2]
      }
      cor.pointsAtlas[[i]]<-automatic.correspondences(cbind(contourAtlas$x[which(contourAtlas$contour.ID==contours[i])],contourAtlas$y[which(contourAtlas$contour.ID==contours[i])]), resLevel, plot=FALSE)
      cor.pointsInput[[i]]<-automatic.correspondences(cbind(contourInput$x[which(contourInput$contour.ID==contoursI[i])],contourInput$y[which(contourInput$contour.ID==contoursI[i])]), resLevel, plot=FALSE)
    }

    cor.points<-list(atlas=cor.pointsAtlas, input=cor.pointsInput)
    
  targetP.x<-numeric()
    referenceP.x<-numeric()
    targetP.y<-numeric()
    referenceP.y<-numeric()
    shape<-numeric()
    
   centroidAtlas<-cor.points$atlas[[1]]$q[1,]
    centroidNorm<-centroidAtlas-cor.points$input[[1]]$q[1,]
    for(i in 1:length(contours) ){
      referenceP.x<-append(referenceP.x, 4*(cor.points$atlas[[i]]$p[,1]-centroidNorm[1]))
      referenceP.y<-append(referenceP.y, 4*(cor.points$atlas[[i]]$p[,2]-centroidNorm[2]))

      targetP.x<-append(targetP.x, 4*cor.points$input[[i]]$p[,1])
      targetP.y<-append(targetP.y, 4*cor.points$input[[i]]$p[,2])
      shape<-append(shape, rep(i, length(cor.points$input[[i]]$p[,2])))
    }

    }else{
      correspondance$correspondance<-na.omit(correspondance$correspondance)
      targetP.x<-correspondance$correspondance[,1]
      targetP.y<-correspondance$correspondance[,2]
      referenceP.x<-correspondance$correspondance[,3]
      referenceP.y<-correspondance$correspondance[,4]
      shape<-correspondance$correspondance[,5]
      centroidNorm<-correspondance$centroidNorm
      centroidAtlas<-correspondance$centroidAtlas
    }
    
    
    #SET OUTPUT FILE
     outputfile<-basename(file)
    outputfile<-strsplit(outputfile, "\\.")[[1]]
    outputfile<-paste(outputfile[-length(outputfile)], collapse='.')

    defaultwd<-getwd()

    if(output.folder=='./'){
      parentpath<-dirname(input)[1]
    }
    if(output.folder=='../'){
      parentpath<-dirname(dirname(input))[1]
    }

    outfolder<-paste('output', outputfile, sep='_')
    setwd(parentpath)
    create.output.directory(outfolder, verbose=FALSE)
    setwd(defaultwd)

    outputfile<-paste(parentpath, outfolder,paste('Registration', outputfile, sep='_'), sep='/')
       if(is.null(filter)){
        MaxDisp<-0
        MinDisp<-0
      }else{
        MaxDisp<-filter$Max
        MinDisp<-filter$Min
        brain.threshold<-filter$brain.threshold
        resize<-filter$resize
      }
      if(is.null(MinDisp)){
        MinDisp<-0
      } 

    resizeP<-resize*4

    
   transformationgrid<-.Call("ThinPlateRegistration", file, targetP.x, targetP.y, referenceP.x, referenceP.y, resizeP, MaxDisp, MinDisp, outputfile)

numPaths<-atlas@summary@numPaths
  scale.factor<-456/97440

  outlines<-list()
  for(i in 1:numPaths){
    
        x<-4*((atlas[[i]]@paths$path@x)*scale.factor-centroidNorm[1])
        y<-4*(((-atlas[[i]]@paths$path@y)*scale.factor+320)-centroidNorm[2])


        index<-cbind(as.integer(round(y)), as.integer(round(x)))
        xT<-(x+(transformationgrid$mx[index]-x) )
        yT<-(y+(transformationgrid$my[index]-y) )


      outlines[[i]]<-list(x  = x, y = y, xT = xT, yT= yT)
  }
  


 img<-paste(outputfile,'_undistorted.png', sep='')
        img <- readPNG(img)

        img = as.raster(img[,])

        par(xaxs='r', yaxs='r')
        plot(c(0, dim(img)[2]*2),c(0, dim(img)[1]), axes=F, asp=1, col=0, xlab='', ylab='', ylim=c(dim(img)[1],0))
        polygon(c(-5, dim(img)[2]+5, dim(img)[2]+5, -5), c(-5, -5, dim(img)[1]+5, dim(img)[1]+5) , col='orange')
        polygon(c(-5+ dim(img)[2], 2*dim(img)[2]+5, 2*dim(img)[2]+5, -5+ dim(img)[2]), c(-5, -5, dim(img)[1]+5, dim(img)[1]+5) , col='purple')
        
        rasterImage(img,0,0, dim(img)[2], dim(img)[1])
        rasterImage(img,dim(img)[2],0, 2*dim(img)[2], dim(img)[1])
        abline(v=dim(img)[2], lwd=2, col='white')
             
    lapply(1:numPaths, function(x){polygon(outlines[[x]]$x, outlines[[x]]$y, border='orange' )})
        lapply(1:numPaths, function(x){polygon(dim(img)[2]+outlines[[x]]$xT,outlines[[x]]$yT, border='purple' )})
        lapply(1:length(targetP.x), function(x){points(c(targetP.x[x]+dim(img)[2],referenceP.x[x]) , c(targetP.y[x],referenceP.y[x]), pch=c(19), col='black', cex=1.8); text(c(targetP.x[x]+dim(img)[2],referenceP.x[x]) , c(targetP.y[x],referenceP.y[x]), label=x, col='white', cex=0.7) } )
    
    
    
    returnlist<-list(atlas=list(outlines=outlines, numRegions=numPaths, col=unlist(lapply(1:numPaths, function(x){atlas[[x]]@paths$path@rgb})) ), transformationgrid=transformationgrid, correspondance=data.frame(targetP.x, targetP.y, referenceP.x, referenceP.y, shape), centroidAtlas=centroidAtlas, centroidNorm = centroidNorm, coordinate=NA, resize= resize, outputfile=outputfile, plane=NA )
    
   return(returnlist)   
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