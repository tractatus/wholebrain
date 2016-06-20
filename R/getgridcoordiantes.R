#' Find dimensions
#'
#' Find dimensions of stitched image from number of tiles.
#' @param x an integer value which give sthe number of tiles.
#' @keywords grid
#' @export
#' @examples
#' find.dim(63)

#find.dim<-function(x){
#    width<-ceiling(sqrt(x))
#    height<-floor(sqrt(x))
#    if( (x-width*height)>0){
#        if( (x-width*height)%%width == 0){
#            height<-height+(x-width*height)/width
#        }else{
#            width<-width+(x-width*height)/height
#        }
#    }
#    
#    if( (x-width*height)<0){
#        if(width>height){
#            area<-width*height
#            width<-width+abs((x-area))/2
#            height<-height-abs((x-area))/2
#        }
#    }
#    
#    return(c(width, height))
#}

find.dim<-function(x){
    possible.integer.values<-(x/seq(1,x))
    suggestions<-possible.integer.values[possible.integer.values%%1 == 0]
    #remove too high and too low sides
    suggestions<-suggestions[!(suggestions==x|suggestions==1)]
    #return a matrix of possible suggestions
    suggestions<-matrix(suggestions[which(suggestions%o%suggestions == x, arr.ind=TRUE)], ncol=2)
    suggestions<-suggestions[!duplicated(rowSums(suggestions)),]
    if(!is.null(dim(suggestions))){
        suggestions<-suggestions[which.min(rowSums(abs(suggestions-sqrt(x)))),]
    }
    suggestions<-sort(suggestions, decreasing=TRUE)
    width<-suggestions[1]
    height<-suggestions[2]
    return(c(width, height))
}


#' Get Grid Coordinates
#'
#' Given a list of file names as input or just a vector this function retrieves the coordinates for each tile as well as the overlaps. It returns a list object with two sublists; width is the image width in pixels, height is the image height in pixels, grid.layout is a data frame where each row is an image and grid.layout$top, grid.layout$bottom, grid.layout$left, grid.layout$right, indicates the position if any at the boundary (0 or 1). grid.layout$x0, grid.layout$x1, grid.layout$y0, grid.layout$y1 gives the coordinates of the tile in the stitched image. The object overlaps is a list with three data frames; overlaps$horizontal give sthe overlap in the horizontal direction of the image plane (right corner of each tile), overlaps$vertical give sthe bottom corner overlap, and overlaps$small gives the small intersection (normally perfectly square) between horizontal and vertical overlap.
#' @param input vector of files, can be eithe rnumeric or character only importance is that length(input) is defined.
#' @param tilesize a single integer which givesthe size of the side of the image tile in pixels (assumes square).
#' @param overlap overlap expressed either as proportion (if overlap<1) of in pixels (if overlap>1).
#' @keywords grid
#' @export
#' @examples
#' grid.coordinates<-get.grid.coordinates(1:63, 2048, plotgrid=F)

get.grid.coordinates<-function(input, tilesize, overlap=0.1, rotate=0, plotgrid=FALSE){

numfiles<-length(input)
if(abs(rotate)<90){
    colTiles<-find.dim(numfiles)[1]
    rowTiles<-find.dim(numfiles)[2]
}else{
    if(abs(rotate)==90){
        colTiles<-find.dim(numfiles)[2]
        rowTiles<-find.dim(numfiles)[1]
    }
}

if(overlap<1){
overlappixels<-(overlap* tilesize)/2
}else{
    overlappixels<-overlap
}

overlappixels<-round(overlappixels)

tileNonoverlap<-tilesize-overlappixels


width<-(colTiles-1)* tileNonoverlap+ tileNonoverlap +overlappixels
height<-(rowTiles-1)* tileNonoverlap+ tileNonoverlap +overlappixels


k=1
grid.layout<-data.frame(top=numeric(), bottom=numeric(), left=numeric(), right=numeric(), x0=numeric(), x1=numeric(), y0=numeric(), y1=numeric())
horizontal.overlaps<-data.frame(leftImage=numeric(), rightImage=numeric(), x0=numeric(), x1=numeric(), y0=numeric(), y1=numeric())
vertical.overlaps<-horizontal.overlaps
small.overlaps<-data.frame(topleftImage=numeric(), toprightImage=numeric(), bottomleftImage=numeric(), bottomrightImage=numeric(), x0=numeric(), x1=numeric(), y0=numeric(), y1=numeric())


if(plotgrid==TRUE){
    plot(c(0, width, width, 0,0), c(0,0, height, height, 0), type='l', ylim=c(height, 0), xlab='', ylab='', axes=F)
    axis(2, las=1)
    axis(1)
    box()
    wholebraincolors<-c('black', '#808080', '#a76776', '#3e6688', '#0e8a8a', '#ff8659')
    wholebraincolors<-rep(wholebraincolors, length.out=colTiles*rowTiles)
}

for (rowTile in 0:(rowTiles-1))
    {
        for (colTile in 0:(colTiles-1))
        {
            #check edge image
            if(rowTile==0){
                top<-1
            }else{
                top<-0
            }
            if(colTile==0){
                left<-1
            }else{
                left<-0
            }
            if(rowTile==(rowTiles-1)){
                bottom<-1
            }else{
                bottom <-0
            }
            if(colTile ==(colTiles-1)){
                right<-1
            }else{
                right <-0
            }
            #check corner image
            
            
            x0<-colTile* tileNonoverlap
            x1<-colTile* tileNonoverlap+ tileNonoverlap +overlappixels
            y0<-rowTile* tileNonoverlap
            y1<-rowTile* tileNonoverlap+ tileNonoverlap +overlappixels
            
            
            if(right==1){
                #do not add any vertical overlap
            }else{
                if(top==1){
                    h.x0<-x1-overlappixels
                    h.x1<-x1
                    h.y0<-y0
                    h.y1<-y1-overlappixels
                 
                }
                if(bottom==1){
                    h.x0<-x1-overlappixels
                    h.x1<-x1
                    h.y0<-y0+overlappixels
                    h.y1<-y1
                }
                if(top==0&bottom==0){
                    h.x0<-x1-overlappixels
                    h.x1<-x1
                    h.y0<-y0+overlappixels
                    h.y1<-y1-overlappixels
                }
                
                horizontal.overlaps <-rbind(horizontal.overlaps, c(k, k+1, h.x0, h.x1, h.y0, h.y1))
            }
            
            if(bottom==1){
                #do not add any vertical overlap
            }else{
                if(left==1){
                    v.x0<-x0
                    v.x1<-x1-overlappixels
                    v.y0<-y1-overlappixels
                    v.y1<-y1
                 
                }
                if(right==1){
                    v.x0<-x0+overlappixels
                    v.x1<-x1
                    v.y0<-y1-overlappixels
                    v.y1<-y1
                }
                if(left ==0& right ==0){
                    v.x0<-x0+overlappixels
                    v.x1<-x1-overlappixels
                    v.y0<-y1-overlappixels
                    v.y1<-y1
                }
                
                vertical.overlaps <-rbind(vertical.overlaps, c(k, k+colTiles, v.x0, v.x1, v.y0, v.y1))
            }
            
            if(right==0&bottom==0){
                s.x0<-x1-overlappixels
                s.x1<-x1
                s.y0<-y1-overlappixels
                s.y1<-y1
                
                small.overlaps <-rbind(small.overlaps, c(k, k+1, k+colTiles, k+colTiles+1, s.x0, s.x1, s.y0, s.y1))
            }
            
            grid.layout <-rbind(grid.layout, c(top,bottom,left,right, x0, x1, y0, y1))

               if(plotgrid==TRUE){
                 polygon(c(x0, x1, x1, x0), c(y0, y0, y1, y1), border= wholebraincolors[k], lwd=2 )
                 if(right==0){
                 polygon(c(h.x0, h.x1, h.x1, h.x0), c(h.y0, h.y0, h.y1, h.y1), col=rgb(0.2,0.2,0.2,0.1))
                 }
                 if(bottom==0){
                 polygon(c(v.x0, v.x1, v.x1, v.x0), c(v.y0, v.y0, v.y1,v.y1), col=rgb(0.9,0.2,0.2,0.1))     
                }
                if(bottom==0&right==0){
                 polygon(c(s.x0, s.x1, s.x1, s.x0), c(s.y0, s.y0, s.y1, s.y1), col=rgb(0.2,0.2,0.9,0.3))    
                }
                 
                 text((x0+x1)/2, (y0+y1)/2, k, col= wholebraincolors[k]  )
                 }
                 
                 k=k+1

            }
            
    }
    
    names(grid.layout)<-c('top' ,   'bottom', 'left',   'right',  'x0' ,    'x1'  ,   'y0'  ,   'y1')
  
    names(horizontal.overlaps)<-c('leftImage',   'rightImage',  'x0' ,    'x1'  ,   'y0'  ,   'y1')
    names(vertical.overlaps)<-c('topImage',   'bottomImage',  'x0' ,    'x1'  ,   'y0'  ,   'y1')
    names(small.overlaps)<-c('topleftImage', 'toprightImage',  'bottomleftImage', 'bottomrightImage',  'x0' ,    'x1'  ,   'y0'  ,   'y1')
    
    layout<-list(width=width, height=height, grid.layout = grid.layout, overlaps=list(horizontal = horizontal.overlaps, vertical = vertical.overlaps, small = small.overlaps) )
    
    return(layout)
    
}