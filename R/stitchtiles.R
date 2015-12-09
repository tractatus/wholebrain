#' Stitch image using grid layout
#'
#' Stitches multiple image tiles.
#' @param input input a character vector consisting of the full path name to 16-bit raw tif image files.
#' @param stitched.image.name filename of the stitched output, default is the directory where the tiles are situated with the added prefix stitched_
#' @param type type of motorized microscope stage. Will define the order of acquisition. Arguments are either row.by.row or default snake.by.row.
#' @param order the order each tile is acquired. Default is left.&.up which starts from bottom left corner and moves up. Other options include right.&.up, left.&.down and right.&.down.
#' @param outputfolder name of output folder to save the image in. Default is ../ which means a directory will be created (if not allready there) in the parent directory to where the tiles are placed. Alternative is either ./ which will create output folder in the current directory in which the tiles are placed. The user might also provide the full system path to an already existing directory.
#' @param tilesize an integer value which gives the length of the side of the tiles in pixels, default is 2048 pixels.
#' @param overlap either  a real-valued numbe ron the unit interval [0,1] designating the proportion of overlap or an integer number which gives the border thickness of the overlap.
#' @param show.image a boolean value, if true the stitched image will be displayed in a display window. Default is false.
#' @param micromanager micromanager defines proportion, or rathe rpercentage, overlap as the thickness of the overlapping borderr ther than proportion of overlapping area. In practise this will only double the border thickness, i.e. 0.1 becomes 0.2 overlap. Default is true.
#' @param verbose boolean value. If true diagnostic output is written to the R console. Deafult is true.
#' @param contrast a real-valued number which gives the contrast 1 means no contrast change >1 increase <1 decerase. Default is increase of contrast with a value of 3.0.
#' @param brightness brightness value, default is 30.
#' @examples
#' #folder where image tiles are stored
#'images<-get.images('/Volumes/microscope/animal001/slide001/section001')
#' #stitch images
#' stitch(images, type = 'snake.by.row', order = 'left.&.up', tilesize=2048, overlap=0.1, show.image=TRUE)

stitch<-function(input, stitched.image.name = 'stitched_{default.folder}.tif', type = 'snake.by.row', order = 'left.&.up', output.folder='../', tilesize=2048, overlap=0.1, show.image=FALSE, micromanager = TRUE, verbose=TRUE, brightness=30, contrast=3.0){
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
    outputfolder<-paste('stitched', basename(dirname(input))[1], sep='_')
    setwd(parentpath)
    create.output.directory(outputfolder)
    setwd(defaultwd)

    output.folder<-paste(parentpath,outputfolder, sep='/')
  }
  if(output.folder=='./'){
    defaultwd<-getwd()
    parentpath<-dirname(input)[1]
    outputfolder<-paste('stitched', basename(dirname(input))[1], sep='_')
    setwd(parentpath)
    create.output.directory(outputfolder)
    setwd(defaultwd)

    output.folder<-paste(parentpath,outputfolder, sep='/')
  }
  outname<-basename(input)
  if(stitched.image.name == 'stitched_{default.folder}.tif'){
    stitched.image.name<-paste(outputfolder, '.tif', sep='')
    #stitched.image.name<-paste(parentpath, stitched.image.name, sep='/')
  }
  
  image.order<-1:length(files)
  
  if(micromanager){overlap<-overlap*2}
  
  grid.coordinates<-get.grid.coordinates(image.order, tilesize, overlap, plotgrid=F)
  numcols<-find.dim(length(files))[1]
  numrows<-find.dim(length(files))[2]
  
  image.grid<-matrix(image.order, nrow= numrows,ncol= numcols, byrow=T)
  
  if(order=='left.&.up'){
  	image.grid <-length(image.grid)-(image.grid-1)
  }
  if(order=='right.&.up'){
  	image.grid <-image.grid[nrow(image.grid):1,]
  }
  if(order=='left.&.down'){
  	image.grid <-image.grid[,ncol(image.grid):1]
  }
  if(type == 'snake.by.row'){
  	image.grid[rev(seq(nrow(image.grid)-1,1,,by=-2)),]<-image.grid[rev(seq(nrow(image.grid)-1,1,,by=-2)),ncol(image.grid):1]
  }


  final.image.order<-as.vector(t(image.grid))
  
  files<-files[final.image.order]
  
  if(overlap<1){
	overlappixels<-(overlap* tilesize)/2
}else{
	overlappixels<-overlap
}

overlappixels<-round(overlappixels)
.Call("LaplacianBlendPipe", files,
   									output.folder, 
   									numrows, 
   									numcols, 
   									overlappixels, 
   									grid.coordinates$width, 
   									grid.coordinates$height, 
   									grid.coordinates$grid.layout$top, 
   									grid.coordinates$grid.layout$bottom, 
   									grid.coordinates$grid.layout$left, 
   									grid.coordinates$grid.layout$right, 
   									grid.coordinates$grid.layout$x0,   
   									grid.coordinates$grid.layout$x1,  
   									grid.coordinates$grid.layout$y0, 
   									grid.coordinates$grid.layout$y1, 
   									grid.coordinates$overlaps$horizontal$leftImage, 
   									grid.coordinates$overlaps$horizontal$x0, 
   									grid.coordinates$overlaps$horizontal$y0, 
   									grid.coordinates$overlaps$horizontal$y1, 
   									grid.coordinates$overlaps$vertical$topImage, 
   									grid.coordinates$overlaps$vertical$bottomImage, 
   									grid.coordinates$overlaps$vertical$x0, 
   									grid.coordinates$overlaps$vertical$x1, 
   									grid.coordinates$overlaps$vertical$y0,
   									grid.coordinates$overlaps$small$topleftImage, 
   									grid.coordinates$overlaps$small$toprightImage, 
   									grid.coordinates$overlaps$small$bottomleftImage, 
   									grid.coordinates$overlaps$small$bottomrightImage,    
   									grid.coordinates$overlaps$small$x0,      
   									grid.coordinates$overlaps$small$y0,
   									as.integer(show.image),
   									as.integer(verbose),
                                    as.numeric(contrast),
                                    as.integer(brightness),
   			stitched.image.name)
}