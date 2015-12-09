#' Flat-field correction
#'
#' An image sensor of a camera is composed of a two dimensional array of light sensitive detectors or pixels. The sesnor array is #'mechanically quite stable with the pixels retaining a rigidly fixed geometric relationship. Each pixel within the array, however, #'has its own unique light sensitivity characteristics. As these characteristics affect camera performance, they must be removed #'through calibration. The process by which a camera is calibrated is known as "Flat Fielding" or "Shading Correction".
#' @param input input a character vector consisting of the full path name to 16-bit raw tif image files. For microscopy the most common systematic error is vignetting where dark shades are stronger among the pixels on the edges.
#' @param outputfolder name of output folder to save the image in. Default is ../ which means a directory will be created (if not allready there) in the parent directory to where the tiles are placed. Alternative is either ./ which will create output folder in the current directory in which the tiles are placed. The user might also provide the full system path to an already existing directory.
#' @param output.prefix the prefix for the generated image tiles, default is FFC_.
#' @param kernel smoothing kernel for the Gaussian smothing of the gain image, higher means more bluring, the number needs to be an odd number default is 301.
#' @param show.image a boolean value, if true the stitched image will be displayed in a display window. Default is false.
#' @param gain.image.name the name of the gain image to be generated default is adding the prefix gain_ to the folder name.
#' @param verbose boolean value. If true diagnostic output is written to the R console. Deafult is true.
#' @examples
#' #folder where image tiles are stored
#'images<-get.images('/Volumes/microscope/animal001/slide001/section001')
#' #stitch images
#' flat.field.correction(images)

flat.field.correction<- function(input, output.folder='../', output.prefix='FFC', kernel=301, show.image=FALSE, gain.image.name = 'gain_{output.folder}.tif', verbose=TRUE){
  files<-character()
  for(i in 1:length(input)){
    file <- as.character(input[i])
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)
    files<-append(files, file)
  }
  if(show.image){show.image<-1}else{show.image<-0}
  show.image<-as.integer(show.image)
  
  if(output.folder=='../'){
    defaultwd<-getwd()
    parentpath<-dirname(dirname(input))[1]
    outputfolder<-paste(output.prefix, basename(dirname(input))[1], sep='_')
    setwd(parentpath)
    create.output.directory(outputfolder)
    setwd(defaultwd)

    output.folder<-paste(parentpath,outputfolder, sep='/')
  }
  if(output.folder=='./'){
    defaultwd<-getwd()
    parentpath<-dirname(input)[1]
    outputfolder<-paste(output.prefix, basename(dirname(input))[1], sep='_')
    setwd(parentpath)
    create.output.directory(outputfolder)
    setwd(defaultwd)

    output.folder<-paste(parentpath,outputfolder, sep='/')
  }
  outname<-basename(input)
  if(gain.image.name == 'gain_{output.folder}.tif'){
    gain.image.name<-paste('gain_', outputfolder, '.tif', sep='')
    gain.image.name<-paste(parentpath, gain.image.name, sep='/')
  }

  .Call("posteriorFFC", files, output.folder, outname, kernel, show.image, gain.image.name, as.integer(verbose))
} 