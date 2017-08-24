#' Prune Cy3 spots from spatial transcriptomics (ST) array
#'
#' This takes segmented contours form Cy3 ST array and prunes it from aritfacts
#' and detects spots thar are missing.function loads a file as a matrix.
#'
#' @param segmentation a segmentation list obtained from segment()
#' @param array.dim the dimensions of the ST array. Default is c(33, 35)
#' @param corner
#' @param tol.dist distance that is tolerated in pixels between the center of two spots. Default is 200 px.
#' @return spot.id two column matrix giving the identity to each spot in row and column coordinates.
#' @return missing two column matrix giving the identity to each spot that couldn't be detected.
#' @examples
#' prune.spots(seg)
#' @export
prune.spots<-function(segmentation, array.dim=c(35, 33), corner=2, tol.dist=200){
  
  sorted.x<-sort(unique(round(segmentation$soma$x,-2)))
  sorted.y<-sort(unique(round(segmentation$soma$y,-2)))
  
  sorted.x<-c(sorted.x[which(diff(sorted.x)>200)],  max(sorted.x))
  sorted.y<-c(sorted.y[which(diff(sorted.y)>200)],  max(sorted.y))
  
  #get the theoretical position of each spot in a grid.
  theoretical.position<-expand.grid(sorted.x, sorted.y)
  #get the unique ID for each theoretical spot in the grid
  spot.index<-expand.grid(seq_along(sorted.x), seq_along(sorted.y))
  
  #match theoretical with actual position by computing the shortest distance returned a spot.id.
  spot.id<-cbind(seq_along(segmentation$soma$x), seq_along(segmentation$soma$x))
  for(i in seq_along(segmentation$soma$x)){
    distance<-sqrt((segmentation$soma$x[i]-theoretical.position$Var1)^2+(segmentation$soma$y[i]-theoretical.position$Var2)^2)
    spot.id[i,]<-as.integer(spot.index[which.min(distance),])
  }
  
  #check which spots are missing
  missing<-spot.index[which(!( paste(spot.index[,1], spot.index[,2])%in%paste(spot.id[,1], spot.id[,2]) )), ]
  tot<-prod(array.dim)
  n.missing<-nrow(missing)
  cat(paste('Number of missing spots:', n.missing, 'out of', tot, '(', round(n.missing/tot*100, 2), '%', ')\n' ) )
  
  #flip to the correct corner
  if(corner==2){
    spot.id[,1]<-abs(spot.id[,1]-c(array.dim[1]+1))
  }
  if(corner==3){
    spot.id[,1]<-abs(spot.id[,1]-c(array.dim[1]+1))
    spot.id[,2]<-abs(spot.id[,2]-c(array.dim[2]+1))
  }
  if(corner==4){
    spot.id[,2]<-abs(spot.id[,2]-c(array.dim[2]+1))
  }  
  
  return(list(spot.id=spot.id, missing=missing))
} 

#' Remove spot outliers
#'
#' This takes segmented cell nuclei from H&E staining and spot contours form Cy3 ST array 
#' and returns a parent vector for each cell nuclei.
#'
#' @param cells to be of H&E nucleu segmentation list obtained from segment()
#' @param array.dim the dimensions of the ST array. Default is c(33, 35)
#' @param corner
#' @param tol.dist distance that is tolerated in pixels between the center of two spots. Default is 200 px.
#' @return spot.id two column matrix giving the identity to each spot in row and column coordinates.
#' @return missing two column matrix giving the identity to each spot that couldn't be detected.
#' @examples
#' prune.spots(seg)
#' @export 
remove.spot.outliers<-function(segmentation, width=26000, height){
    if(missing(height)){
        height<-width
    }
    if(missing(width)){
        width<-height
    }
    
    
    poly<-list()
    poly$x<-c(median(segmentation$soma$x)+c(-width/2,+width/2,+width/2,-width/2))
    poly$y<-c(median(segmentation$soma$y)+c(-height/2,-height/2,+height/2,+height/2))
    
    remove<-which(point.in.polygon(segmentation$soma$x, segmentation$soma$y, poly$x, poly$y)==0)

    if(length(remove)>0){
    
    segmentation$soma$x<-segmentation$soma$x[-remove]
    segmentation$soma$y<-segmentation$soma$y[-remove]
    segmentation$soma$intensity<-segmentation$soma$intensity[-remove]
    segmentation$soma$area <-segmentation$soma$area[-remove]
    
    remove.contours<-which(segmentation$soma$contour.ID%in%(remove-1))
    
    segmentation$soma$contour.ID<-segmentation$soma$contour.ID[-remove.contours]
    segmentation$soma$contour.x<-segmentation$soma$contour.x[-remove.contours]
    segmentation$soma$contour.y<-segmentation$soma$contour.y[-remove.contours]
    }

    return(segmentation)
    
}

merge.st.data<-function(datasets){
  cat(paste0('Merging dataset: ', 1,'/', length(datasets), '\r'))
  load(datasets[1])
  master.data<-data.frame(cbind(dataset$spots, dataset$genes))
  for(i in seq_along(datasets)[-1]){
    load(datasets[i])
    master.data.tmp<-data.frame(cbind(dataset$spots, dataset$genes))
    master.data<-rbind.fill(master.data, master.data.tmp)
    cat(paste0('Merging dataset: ', i,'/', length(datasets), '\r'))
  }
  return(master.data)
}

combine.st.data<-function(spots, registration, stdata, nuclei, corner=1){
    spots<-remove.spot.outliers(spots)
    spot.id<-prune.spots(spots, corner=corner)

    cat("1/4 reading ST-data...\n")
    rnaseq<-read.table(stdata)

    rna.seq.id<-strsplit(row.names(rnaseq),'x')

    rna.seq.id<-do.call(rbind, rna.seq.id)
    rna.seq.id<-apply( rna.seq.id, 2, as.numeric )
    rna.seq.id<-rna.seq.id[,c(2,1)]

    temp.rna.seq.id<-apply(rna.seq.id, 1, function(x)paste(x, collapse="x"))       
    rna.seq.index<-integer()
    for(i in 1:nrow(spot.id$spot.id)){
        spot.ROI<-paste(spot.id$spot.id[i,],collapse="x")
        if(spot.ROI%in%temp.rna.seq.id){
            rna.seq.index<-append(rna.seq.index, which( temp.rna.seq.id == spot.ROI) )
        }else{
            rna.seq.index<-append(rna.seq.index, NA)
        }
    }
    cat('2/4 Transforming spots to atlas \n')
    dataset<-get.cell.ids(registration, spots, forward.warp = TRUE)
    cat('3/4 Transforming nuclei to atlas \n') 
    nuclei.dataset<-get.cell.ids(registration, nuclei, forward.warp = TRUE) 
    dataset$spot.id<-seq(1,nrow(dataset))
    nuclei.dataset$spot.id<-rep(NA,nrow(nuclei.dataset))
    dataset$nuclei<-rep(0,nrow(dataset))
    k<-1
    for(i in unique(spots$soma$contour.ID)){
        cat(paste('4/4 Counting nuclei in spot: ',i+1,'\r'))
        cell.in.spot<-point.in.polygon(nuclei$soma$x, nuclei$soma$y, spots$soma$contour.x[spots$soma$contour.ID==i], spots$soma$contour.y[spots$soma$contour.ID==i])
        dataset$nuclei[k]<-sum(cell.in.spot)
        if(dataset$nuclei[k]>0){
            nuclei.dataset$spot.id[which(cell.in.spot==1)]<-k
        }
        k<-k+1
    }
    dataset<-dataset[!is.na(rna.seq.index),]

    rna.seq.index<-na.omit(rna.seq.index)
    rna.seq.index<-as.integer(rna.seq.index)

  rnaseq<-rnaseq[rna.seq.index,]
  master.data<-list(spots=dataset, genes=rnaseq, nuclei=nuclei.dataset, corner=corner)
  return(master.data)
}

st.plot.quality<-function(){

}

spatial.transcriptomics<-function(anatomy, cDNA, feature.filter, cell.body.feature, plot=FALSE){
    anatomyTIF<-rgb2gray(anatomy)
    cDNATIF<-rgb2gray(cDNA, invert=FALSE)
    features<-segment(cDNATIF, filter= feature.filter, display=FALSE)
    print(paste('Features found: ', length(features$soma$x)))
    features$soma$x<-features$soma$x*4
    features$soma$y<-features$soma$y*4
    print('Segmenting out nuclei. this might take several minutes...')
    cell.bodies<-segment(anatomyTIF, filter= cell.body.feature, display=FALSE)
    print(paste('Cell nuclei found: ', length(cell.bodies$soma$x)))
    if(plot){
        par(yaxs='i', xaxs='i')
        plot(round(features$soma$x*(1/(0.5))), round(features$soma$y*(1/(0.5)), -2), asp=1, ylim=c(32768*2,0), xlim=c(0,2*29184), axes=F, xlab='', ylab='', bg=rgb(1,0,0,0.2), pch=21)
        points(cell.bodies$soma$x*2, cell.bodies$soma$y*2, pch=16, cex=0.25)
    }
    
    return(list(features =features, cell.bodies =cell.bodies, filename=anatomyTIF))
}