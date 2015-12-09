conv <- function(a, b)
       .C("convolve",
          as.double(a),
          as.integer(length(a)),
          as.double(b),
          as.integer(length(b)),
          ab = double(length(a) + length(b) - 1))$ab

showImage <- function(file) {
  file <- as.character(file)[1]
  ## check for existence
  if(!file.exists(file))
      stop(file, "not found")
  ## expand path
  file <- path.expand(file)
  a <- .Call("loadImage", file)
  cat("SHOWING...\n")
  .Call("showImage", a)
  return(a)
}

liftingscheme <- function(file) {
    file <- as.character(file)[1]
    ## check for existence
    if(!file.exists(file))
    stop(file, "not found")
    ## expand path
    file <- path.expand(file)
    a <- .Call("liftScheme", file)
    return(a)
}

#create output directories
create.output.directory<-function(subDir, mainDir=getwd()){

    if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
        #cat(paste(subDir, "allready exists\n"))
    } else if (file.exists(paste(mainDir, subDir, sep = "/", collapse = "/"))) {
      cat(paste(subDir, "exists in", mainDir, "but is a file\n"))
      # you will probably want to handle this separately
    } else {
      cat(paste(subDir, "folder created\n"))
      dir.create(file.path(mainDir, subDir))
    }

    if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
      # By this point, the directory either existed or has been successfully created
      #setwd(file.path(mainDir, subDir))
    } else {
      cat(paste(subDir, "does not exist\n"))
      # Handle this error as appropriate
    }
}

createfilter<- function(input, numthresh=8){
  file <- as.character(input)[1]
  ## check for existence
  if(!file.exists(file))
    stop(file, "not found")
  file <- path.expand(file)
  a<-.Call("interactiveThreshold", file, numthresh)
  return(a)
}    

stackapply<- function(input){
  files<-character()
  for(i in 1:length(input)){
    file <- as.character(input[i])
    ## check for existence
    if(!file.exists(file))
      stop(file, "not found")
    file <- path.expand(file)
    files<-append(files, file)
  }
  .Call("runonstack", files, file, file)
}


#flat.field.correction<- function(input, output.folder='../', output.prefix='FFC', kernel=301, show.image=FALSE, gain.image.name = 'gain_{output.folder}.tif'){
#  files<-character()
#  for(i in 1:length(input)){
#    file <- as.character(input[i])
#    ## check for existence
#    if(!file.exists(file))
#      stop(file, ", file not found")
#    file <- path.expand(file)
#    files<-append(files, file)
#  }
#  if(show.image){show.image<-1}else{show.image<-0}
#  show.image<-as.integer(show.image)
#  
#  if(output.folder=='../'){
#    defaultwd<-getwd()
#    parentpath<-dirname(dirname(input))[1]
#    outputfolder<-paste(output.prefix, basename(dirname(input))[1], sep='_')
#    setwd(parentpath)
#    create.output.directory(outputfolder)
#    setwd(defaultwd)
#
#    output.folder<-paste(parentpath,outputfolder, sep='/')
#  }
#  outname<-basename(input)
#  if(gain.image.name == 'gain_{output.folder}.tif'){
#    gain.image.name<-paste('gain_', outputfolder, '.tif', sep='')
#    gain.image.name<-paste(parentpath, gain.image.name, sep='/')
#  }
#
#  .Call("posteriorFFC", files, output.folder, outname, kernel, show.image, gain.image.name)
#} 


imgSWT <- function(input, filter.parameters=NULL, scales=6, cell.bodies = 3, processes = 5, family='db2', sigma=10, processLength=12, alim=c(200, 500), pch=21, bg="white", cex=2, lwd=2, illustrator=F, output='waveletoutput') {
    file <- as.character(input)[1]
    family <- as.character(family)[1]
    outputfile <- as.character(output)[1]
    scales<-as.integer(scales)
    cellBodies<-as.integer(cell.bodies)
    processes <-  as.integer(processes)
    sigmaR<-as.integer(sigma)
    processLength<-as.integer(processLength)
    areaSize<-alim
    # initialize the filter parameter list.
    if(is.null(filter.parameters)){
      noFilter<-1

      filter_minArea<-as.integer(0)
      filter_maxArea<-as.integer(0)
      filter_minThresh<-as.integer(0)
      filter_maxThresh<-as.integer(0)
      filter_eccentricity<-as.integer(0)
      filter_eccentricity<-as.integer(0)
      filter_alpha<-as.numeric(1.0)
      filter_beta<-as.integer(0)

    }else{
      noFilter<-0

      filter_minArea<-as.integer(filter.parameters$alim[1])
      filter_maxArea<-as.integer(filter.parameters$alim[2])
      filter_minThresh<-as.integer(filter.parameters$threshold.range[1])
      filter_maxThresh<-as.integer(filter.parameters$threshold.range[2])
      filter_eccentricity<-as.integer(filter.parameters$eccentricity)
      filter_alpha<-as.numeric(filter.parameters$alpha)
      filter_beta<-as.integer(filter.parameters$beta)
      filter_Max<-as.numeric(filter.parameters$Max)
      filter_Min<-as.integer(filter.parameters$Min)
    }


    ## check for existence
    if(!file.exists(file))
    stop(file, "not found")
    
    create.output.directory('approximations_coefficents')
    lapply(1:scales-1, function(x){create.output.directory( paste('d',x, sep='') ) } )
    create.output.directory('trace')
    create.output.directory('coherency')
    create.output.directory('orientation')
    create.output.directory('scharr')
    create.output.directory('cellbodies')
    create.output.directory('processes_labeled')
    create.output.directory('processes_anisotropic')
    create.output.directory('output')
    create.output.directory('blended')
    create.output.directory('analyzed')
    create.output.directory('rangecorrected')
    

    ## expand path
    file <- path.expand(file)
    #############################
    #           CALL            #
    #############################
    a <- .Call("filterImage", file, noFilter, filter_minArea, filter_maxArea, filter_minThresh, filter_maxThresh, filter_eccentricity, filter_alpha, filter_beta, filter_Max, filter_Min, scales, cellBodies, processes, family, sigmaR, areaSize, processLength, outputfile)
    cat("IMAGE PROCESSED\n")
    #############################
    #           PLOT            #
    #############################
    confocal.image<-readPNG(paste('output/d', cell.bodies, '_output_', output, '.png', sep='') )
    confocal.image<-as.raster(confocal.image[,,1:3])
    process.color<-rgb(a$process.r,  a$process.g, a$process.b, 255, maxColorValue = 255)

    if(illustrator){
      if(dim(confocal.image)[2]==dim(confocal.image)[1]){
        quartz(width=8, height=8)
        par(xaxs='i', yaxs='i', pty='s')
      }else{
        if(dim(confocal.image)[2]<dim(confocal.image)[1]){
          quartz(width=8*dim(confocal.image)[2]/dim(confocal.image)[1], height=8)
        }else{
          quartz(width=8, height=8*dim(confocal.image)[1]/dim(confocal.image)[2])
        }
        par(xaxs='i', yaxs='i', pty='m')
      }
      
      par(mar=c(0,0,0,0))
    
      plot(c(0,dim(confocal.image)[2]), c(0,dim(confocal.image)[1]), axes=F, ylab='', xlab='', col=0)
      rasterImage(confocal.image,0,0,dim(confocal.image)[2],dim(confocal.image)[1])
      points(a$process.x, dim(confocal.image)[1] - a$process.y, col=process.color, pch=16, cex=0.35) #cex=0.15
    
      dev.copy(png,  width = dim(confocal.image)[2], height = dim(confocal.image)[1], paste('analyzed/', output, '_illustrator.png', sep=''))
      dev.off()
    }

    if(dim(confocal.image)[2]==dim(confocal.image)[1]){
      quartz(width=8, height=8)
      par(xaxs='i', yaxs='i', pty='s')
    }else{
      if(dim(confocal.image)[2]<dim(confocal.image)[1]){
        quartz(width=8*dim(confocal.image)[2]/dim(confocal.image)[1], height=8)
      }else{
        quartz(width=8, height=8*dim(confocal.image)[1]/dim(confocal.image)[2])
      }
      par(xaxs='i', yaxs='i', pty='m')
    }

    par(mar=c(0,0,0,0))
    
    plot(c(0,dim(confocal.image)[2]), c(0,dim(confocal.image)[1]), axes=F, ylab='', xlab='', col=0)
    rasterImage(confocal.image,0,0,dim(confocal.image)[2],dim(confocal.image)[1])
    points(a$process.x, dim(confocal.image)[1] - a$process.y, col=process.color, pch=16, cex=0.15)  

    points(a$centroid.x[a$area>areaSize[1]], dim(confocal.image)[1] - a$centroid.y[a$area>areaSize[1]], bg=bg, pch=pch, cex=cex, lwd=lwd)
    dev.copy(png,  width = dim(confocal.image)[2], height = dim(confocal.image)[1], paste('analyzed/', output, '_analyzed.png', sep=''))
    dev.off()

    if(illustrator){
      confocal.image<-readPNG(paste('analyzed/', output, '_illustrator.png', sep='') )
      confocal.image<-as.raster(confocal.image[,,1:3])
      plot(c(0,dim(confocal.image)[2]), c(0,dim(confocal.image)[1]), axes=F, ylab='', xlab='', col=0)
      rasterImage(confocal.image,0,0,dim(confocal.image)[2],dim(confocal.image)[1])
      points(a$centroid.x[a$area>areaSize[1]], dim(confocal.image)[1] - a$centroid.y[a$area>areaSize[1]], bg=bg, pch=pch, cex=cex*1.1, lwd=lwd)
      dev.copy2pdf(file=paste('analyzed/', output, '_illustrator.pdf', sep=''))
    }  
    dev.off()

    a$process.type[which(a$process.type==1)]<-"slab"
    a$process.type[which(a$process.type==2)]<-"end.point"
    a$process.type[which(a$process.type==3)]<-"joint3"
    a$process.type[which(a$process.type==4)]<-"joint4"
    a$process.type[which(a$process.type==5)]<-"jointNA"

    a$eccentricity[which(a$eccentricity==0)]<-NA #contours with less than 5 pixels cannot have eccenticity therefore lets put them to NA instead of 0.

    #reorganize output
    somas<-list(centroid.x=a$centroid.x, centroid.y=a$centroid.y, area=a$area, intensity=a$intensity, eccentricity=a$eccentricity, id=a$id, contour.x=a$x, contour.y=a$y)
    processes<-list(id=a$process.id, type=a$process.type, orientation=a$process.theta, x=a$process.x,  y=a$process.y, col= process.color, coherence= a$process.coherence)
    b<-list(somas=somas, processes=processes)

    return(b)
}

unstitch<-function(img, tilesize, overlap=NULL, position='bottomright'){
    file <- as.character(img)[1]
    if(!file.exists(file))
    stop(file, "not found")
    
    
    ## expand path
    file <- path.expand(file)
    outputfile<-basename(file)
    outputfile<-sub("^([^.]*).*", "\\1", outputfile)
    create.output.directory(paste('Tiled', outputfile, sep='_'))
    if(is.null(overlap)){overlap<-(-999)}
    if(position=='topleft'){pos=1}
    if(position=='topright'){pos=2}
    if(position=='bottomright'){pos=3}
    if(position=='bottomleft'){pos=4}
    a <- .Call("createTiles", file, tilesize, overlap, pos, outputfile)
    return(a)
}


laplacetest<-function(){
  .Call("testLaplace")
}


browsersection<-function(imgfilename, adjustments, neurons){
    file <- as.character(img)[1]
    if(!file.exists(file))
    stop(file, "not found")
    ## expand path
    file <- path.expand(file)
    a <- .Call("Zoomify", file)
    return(a)
}


getcontours<-function(img, threshold){
    file <- as.character(img)[1]
    if(!file.exists(file))
    stop(file, "not found")
    
    
    ## expand path
    file <- path.expand(file)
    
     a <- .Call("getContours", file, thresh)
     return(a)
}