data(glassbrain, envir=environment())
data(EPSatlas, envir=environment())
data(atlasIndex, envir=environment())
#data(atlasOntology, envir=environment())
data(ontology, envir=environment())
data(schematics, envir=environment())


suggestions<-function(dataset, normalize.by=NULL, exclude.below=10, reduce.below=100, clusters=7){
dataset$acronym<-as.character(dataset$acronym)
tableCount<-table(dataset$acronym, dataset$animal)
tableCount <-tableCount[order(tableCount[,1], decreasing=TRUE),]
average.cells<-sort(rowMeans(tableCount), decreasing=TRUE)

if(length(which(dataset$acronym%in%c('fiber tracts', 'grey', names(which(average.cells<exclude.below)) ) ))){
dataset<-dataset[-which(dataset$acronym%in%c('fiber tracts', 'grey', names(which(average.cells<exclude.below)) ) ), ]
}
tableCount<-table(dataset$acronym, dataset$animal)
tableCount <-tableCount[order(tableCount[,1], decreasing=TRUE),]
average.cells<-sort(rowMeans(tableCount), decreasing=TRUE)

while(length(which(average.cells<reduce.below))>0){


to.be.replaced<-names(average.cells[average.cells<reduce.below])
for(i in to.be.replaced){
    dataset$acronym[which(dataset$acronym==i)]<- as.character( get.acronym.parent(i) ) 
    # dataset$acronym[which(dataset$acronym==i)]<-
}


tableCount<-table(dataset$acronym, dataset$animal)
tableCount <-tableCount[order(tableCount[,1], decreasing=TRUE),]
average.cells<-sort(rowMeans(tableCount), decreasing=TRUE)


}

if( length(which(row.names(tableCount)=='root')) > 0){
tableCount<-tableCount[-which(row.names(tableCount)=='root'),]
}

print(tableCount)
print(row.names(tableCount))

tableCount<-tableCount[order(row.names(tableCount)), ]

#group<-color.from.acronym(row.names(tableCount))

#d <- dist(t(col2rgb(group)), method = "euclidean") # distance matrix
#fit <- hclust(d, method="ward.D")

#tableCount <-tableCount[rev(fit$order),]

group<-unlist(lapply(row.names(tableCount), get.sup.structure))


new.order<-c('CTXpl','CTXsp','STR','PAL','TH','HY','MBsen','MBmot','MBsta','P','MY')

for(i in 1:length(new.order)){
    group[ which(group==new.order[i]) ]<-i
}
group[which(!group%in%as.character(c(1:length(new.order))))]<-(length(new.order)+1)
group<-as.numeric(group)
tableCount<-tableCount[order(group),]
group<-group[order(group)]

tableCount2<-tableCount
for(j in unique(group)){
    
    tableCount2[which(group==j), ]<-tableCount[which(group==j)[order( rowMeans(tableCount)[which(group==j)], decreasing=TRUE )],]
    row.names(tableCount2)[which(group==j)]<-row.names(tableCount)[which(group==j)[order( rowMeans(tableCount)[which(group==j)], decreasing=TRUE )]]
}
    return(tableCount2)

}    


se<-function (x, na.rm = TRUE) 
sqrt(var(x, na.rm = na.rm)/length(x[complete.cases(x)]))

plot.suggestions<-function(dataset, normalize.by=NULL, group = NULL, title='Groups:', color=gray(0.4), exclude.below=10, reduce.below=100, exclude.regions=NULL, include.regions=NULL, xlab='Cell count', log.scale=TRUE, bargraph=FALSE, fun = function(x) mean(x, 
        na.rm = TRUE), ci.fun = function(x) c(fun(x) - 1.96*se(x), 
        fun(x) + 1.96*se(x)), xlim=NULL, device=TRUE){

counts<-suggestions(dataset, exclude.below=exclude.below, reduce.below=reduce.below)

if(!is.null(exclude.regions)){
    remove<-which(row.names(counts)%in%exclude.regions)
    if(length(remove)>0){
        counts<-counts[-remove,]
    }
}

if(!is.null(include.regions)){
    dataset$acronym<-as.character(dataset$acronym)
    tableCount<-table(dataset$acronym, dataset$animal)

    to.be.included<-lapply(include.regions, function(x){ colSums(tableCount[which(row.names(tableCount)%in%get.sub.structure(x)), ]) } )

    to.be.included<-do.call(rbind, to.be.included)
    row.names(to.be.included)<-include.regions

    counts<-rbind(counts, to.be.included)
}


if(!is.null(normalize.by)){
    if(length(normalize.by)!=ncol(counts)){
        print(paste('Error: You have ', ncol(counts), ' animals, but entered values ', length(normalize.by),' for normalization.'))
        return()
    }
    
    for(j in 1:ncol(counts)){
        counts[,j]<-counts[,j]/normalize.by[j]
    }
}

if(log.scale){
counts <- log10(counts)
}

        if (device) {
            quartz(width = 7.036585, height = 0.2099039 * nrow(counts))
        }
        layout(matrix(c(1, 1, 1, 2, 2, 2, 2), nrow = 1))
        par(mar = c(4, 0, 4, 0))
        plot(rep(2.5, nrow(counts)), nrow(counts):1, col = 0, 
            axes = F, ylim = c(0.5, nrow(counts) + 0.5), ylab = "", 
            xlab = "", xlim = c(1, 5))
        mtext("Input region:", 3, cex = 0.9)
        for (i in 1:nrow(counts)) {
            regioncolor <- color.from.acronym(row.names(counts)[i])
            regioncolor <- adjustcolor(regioncolor, alpha.f = 0.2)
            y.lab <- (nrow(counts) + 1) - i
            polygon(c(1, 5, 5, 1), c(y.lab - 0.5, y.lab - 0.5, 
                y.lab + 0.5, y.lab + 0.5), col = regioncolor, 
                border = FALSE)
            text(3, y.lab, name.from.acronym(row.names(counts)[i]), 
                cex = 0.9)
        }
        
        
        if(!is.null(group)){
            par(xpd = TRUE)
        legend(3, -0.75, sort(unique(group)), pch = c(21), pt.bg = color, title = title, bg = "white", 
            horiz = TRUE, cex = 1.3, xjust = 0.5)
        par(xpd = FALSE)
        }
        
par(mar = c(4, 4, 4, 6))
        zeros <- floor(range(counts[is.finite(counts)])[1])
        no.zero.values<-FALSE
        print(which(!is.finite(counts), arr.ind=TRUE))
        print( zeros)
        if(length(which(!is.finite(counts)))==0){
            no.zero.values<-TRUE
        }
        counts[!is.finite(counts)] <- zeros
        if(is.null(xlim)){
                x.range <- c( floor(range(counts[is.finite(counts)])[1]), ceiling(range(counts[is.finite(counts)])[2]) )
        }else{
                x.range <- xlim
            }
        plot(apply(counts, 1, max), nrow(counts):1 - 0.125, pch = 21, 
            bg = "white", ylim = c(0.5, nrow(counts) + 0.5), 
            xlim = x.range, xlab = "", axes = F, ylab = "", col = 0)
        for (i in 1:nrow(counts)) {
            regioncolor <- color.from.acronym(row.names(counts)[i])
            regioncolor <- adjustcolor(regioncolor, alpha.f = 0.15)
            y.lab <- (nrow(counts) + 1) - i
            polygon(c(x.range[1]-1, x.range[2]+1, x.range[2]+1, x.range[1]-1), c(y.lab - 0.5, 
                y.lab - 0.5, y.lab + 0.5, y.lab + 0.5), col = regioncolor, 
                border = FALSE)
        }
        
        if(log.scale){
        
        log.range <- 10^seq(x.range[1], x.range[2])
        if(no.zero.values){
            axis(1, at = seq(x.range[1], x.range[2]), las = 1, labels = c( 
            log.range))
            axis(3, at = seq(x.range[1], x.range[2]), las = 1, labels = c( 
            log.range))
            }else{
        axis(1, at = seq(x.range[1], x.range[2]), las = 1, labels = c(0, 
            log.range[-1]))
        axis(3, at = seq(x.range[1], x.range[2]), las = 1, labels = c(0, 
            log.range[-1]))
        }
        log.range <- unlist(lapply(1:(length(log.range) - 1), 
            function(x) {
                seq(log.range[x], log.range[x + 1], by = log.range[x])
            }))
        axis(1, at = log10(log.range), labels = FALSE)
        axis(3, at = log10(log.range), labels = FALSE)
        axis(2, at = nrow(counts):1, labels = row.names(counts), 
            las = 1)
        axis(4, at = nrow(counts):1, labels = row.names(counts), 
            las = 1)
        abline(v = log10(log.range), col = "lightblue")
                    abline(h = 1:nrow(counts), lty = 2, col = "gray")

        }else{
            abline(h = 1:nrow(counts), lty = 2, col = "gray")
            axis(2, at = nrow(counts):1, labels = row.names(counts), 
            las = 1)
        axis(4, at = nrow(counts):1, labels = row.names(counts), 
            las = 1)
            abline(v = seq(x.range[1],x.range[2], length.out=6), col = "lightblue")

        }

    if(is.null(group)){

        if(bargraph){
            #barfunction
                lapply(1:nrow(counts), function(x) {
            polygon(c(x.range[1], rep(fun(counts[x, ]), 2), x.range[1]), c(rep(nrow(counts) - x + 1, 4 ) + c(-0.25,-0.25,0.25,0.25) ),  col=color );
           
              lines(ci.fun(counts[x, ]), rep(nrow(counts) - x + 1, 2), lwd=1.5);
               points(fun(counts[x, ]), nrow(counts) - x + 1, 
                pch = 21, bg = color, cex = 1.2)
            })

            }else{
                #pointfunction
lapply(1:nrow(counts), function(x) {
            points(counts[x, ], rep(nrow(counts) - x + 1, ncol(counts) ), 
                pch = 21, bg = color, cex = 1.2)})
            }
        }else{
            #GROUPS
            k<-1
            vertical<-seq(0.15,-0.15,length.out=length(unique(group)))
            for(j in sort(unique(group))){
                if(bargraph){
                    lapply(1:nrow(counts), function(x) {
                        polygon(c(x.range[1], rep(fun(counts[x,  which(group==j)]), 2), x.range[1]), c(rep(nrow(counts) - x + 1 + vertical[k], 4 ) + c(-0.25,-0.25,0.25,0.25)/length(unique(group)) ),  col=color[k] );
                        lines(ci.fun(counts[x,  which(group==j)]), rep(nrow(counts) - x + 1  + vertical[k], 2), lwd=1.5);
                        points(fun(counts[x, which(group==j) ]), nrow(counts) - x + 1  + vertical[k], pch = 21, bg = color[k], cex = 1.2)
                    })
                }else{
                    lapply(1:nrow(counts), function(x) {
                        points(counts[x, which(group==j)], rep(nrow(counts) - x + 1 + vertical[k], length(which(group==j)) ), pch = 21, bg = color[k], cex = 1.2)
                    })
                }
                k<-k+1
            }
        }
        
                box()
                if(log.scale){
        if(!no.zero.values){
        par(xpd = TRUE)
        polygon(c( mean(log10(log.range)[1:2]), log10(log.range)[2], log10(log.range)[2], mean(log10(log.range)[1:2])), c(-15, -15, nrow(counts) + 
            15, nrow(counts) + 15), col = "white", border = "white")
        
        #polygon(c(x.range[1] + 1/x.range[2], x.range[1] + 1/x.range[2]/2, x.range[1] + 
        #    1/x.range[2]/2, x.range[1] + 1/x.range[2]), c(-15, -15, nrow(counts) + 
        #    15, nrow(counts) + 15), col = "white", border = "white")
        par(xpd = FALSE)
        abline(v = c(mean(log10(log.range)[1:2]), log10(log.range)[2]))
        #abline(v = c(x.range[1] + 1/x.range[2], x.range[1] + 1/x.range[2]/2))
        }
        }
        mtext(xlab, 3, 2.2, cex = 0.8)
        mtext(xlab, 1, 2.2, cex = 0.8)
        
        

}


schematic.plot.sagittal<-function(regi, dataset, pch=21, cex=0.5, col='black'){
      EPSatlas<-SAGITTALatlas
      atlasIndex<-atlasIndex[atlasIndex$plane=="sagittal", ]
      plate.width<-1.159292
      width<-(528-214)*2
      coordinate<-regi$coordinate

  k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))

  scale.factor<-mean(c(dim(regi$transformationgrid$mx)[1]/regi$transformationgrid$height,    dim(regi$transformationgrid$mx)[2]/regi$transformationgrid$width) )
regi<-get.forward.warpRCPP(regi)

index<-round(scale.factor*cbind(dataset$y, dataset$x))
    #ensure that 0 indexes are 1
    index[index==0]<-1
    #check if point is outside image
    if( length( which(index[,1]>dim(regi$transformationgrid$mxF)[1]) ) ){
        index[which(index[,1]>dim(regi$transformationgrid$mxF)[1]) ,1]<-dim(regi$transformationgrid$mxF)[1]
    }
    if( length( which(index[,2]>dim(regi$transformationgrid$mxF)[2]) ) ){
      index[which(index[,2]>dim(regi$transformationgrid$mxF)[2]),2]<-dim(regi$transformationgrid$mxF)[2]
    }

    somaX<-regi$transformationgrid$mxF[index]/scale.factor
    somaY<-regi$transformationgrid$myF[index]/scale.factor
    
    atlas.fit.scale<-diff(range(regi$atlas$outlines[[1]]$yl/scale.factor))
    atlas.y.range<-range(EPSatlas$plates[[k]][[1]]@paths$path@y)
    atlas.scale<-diff(atlas.y.range)
    
    distance.to.atlas.center<-(atlas.y.range-68234.56/2)
    
    ypos<- -(range(regi$atlas$outlines[[1]]$yl/scale.factor) - (atlas.fit.scale/atlas.scale)*distance.to.atlas.center )[1]
    ypos<-scale.factor*ypos/4
    
    xpos<-min(-(EPSatlas$plates[[k]][[1]]@paths$path@x- (plate.width*97440)/2)+(plate.width*97440)/2)/(plate.width*97440)* plate.width*456 - scale.factor*min(regi$atlas$outlines[[1]]$xl/scale.factor)/4
    #xpos<-scale.factor*(xmin*(atlas.fit.scale/atlas.scale) )
    
    #somaY<-somaY-ypos
    
    x<-scale.factor*(somaX)/4
    y<-scale.factor*(somaY)/4
    y<-y+ypos
    x<-x+xpos
    x<-(x-(width)/2)*25/1000
    y<-(-y/1000*25)
           
    plot(x,y, col=0, xlim=c(-10,7), ylim=c(-4.5,4), axes=FALSE, ylab='', xlab='', asp = 1)
    
        
    style<-rep(gray(0.95), regi$atlas$numRegions)
                
                lapply(1:regi$atlas$numRegions, function(N) {
                    yPoly<- scale.factor*(regi$atlas$outlines[[N]]$yl/scale.factor)/4
    yPoly <-yPoly +ypos
    yPoly <-(-yPoly/1000*25)

    xPoly<- scale.factor*(regi$atlas$outlines[[N]]$xl/scale.factor)/4
    xPoly <-xPoly +xpos

    xPoly <-(xPoly-(width)/2)*25/1000

                    
                  polygon(xPoly, yPoly, 
                    col = style[N], border = "black", lwd = 1, 
                    lty = 3)
                                    })
                                    
                                    
  ventricles <- c(which(EPSatlas$plate.info[[k]]$style == 
                  "#aaaaaa")) 
   
   for (i in ventricles) {
                    yPoly<- scale.factor*(regi$atlas$outlines[[i]]$yl/scale.factor)/4
    yPoly <-yPoly +ypos
    yPoly <-(-yPoly/1000*25)

    xPoly<- scale.factor*(regi$atlas$outlines[[i]]$xl/scale.factor)/4
    xPoly <-xPoly +xpos

    xPoly <-(xPoly-(width)/2)*25/1000

                    
                  polygon(xPoly, yPoly, 
                      col = "black", border = "black", lwd = 4)

                  }
                
                fibertracts <- c(which(EPSatlas$plate.info[[k]]$style == 
                  "#cccccc"))
                for (i in fibertracts) {
                    yPoly<- scale.factor*(regi$atlas$outlines[[i]]$yl/scale.factor)/4
    yPoly <-yPoly +ypos
    yPoly <-(-yPoly/1000*25)

    xPoly<- scale.factor*(regi$atlas$outlines[[i]]$xl/scale.factor)/4
    xPoly <-xPoly +xpos

    xPoly <-(xPoly-(width)/2)*25/1000

                    
                  polygon(xPoly, yPoly, , 
                    col = gray(0.85), border = "black", lwd = 1)
                }                                 

polygon(c(5,6,6,5),  c(-4.2, -4.2, -4.1, -4.1), col='black' )

points(x,y, pch=pch, cex=cex, col=col, bg=dataset$color) }



schematic.plot<-function (dataset, coordinate = NULL, plane='coronal', title = TRUE, mm.grid = TRUE, 
    save.plots = FALSE, dev.size = c(5.4, 4.465), pch=21, cex=0.5, col='black', scale.bar=FALSE, region.colors=FALSE, device = TRUE) 
{


if(plane=="sagittal"){

      schematic.plot.sagittal(regi, dataset, pch = pch, cex=cex, col=col) #really really ugly hack
      return()
}

    if ( (!save.plots)&(device) ) {
        quartz(width = dev.size[1], height = dev.size[1])
    }
    #assign temporary animal label if this is not defined
    dataset$animal[is.na(dataset$animal)]<-'noname'

    if (is.null(coordinate)) {
        if (length(which(dataset$color == "#000000")) > 0) {
            dataset <- dataset[-which(dataset$color == "#000000"), 
                ]
        }
        datasets <- dataset
        for (animal.index in unique(datasets$animal)) {
            coordinates <- sort(unique(datasets$AP[which(datasets$animal == 
                animal.index)]), decreasing = TRUE)
            for (q in 1:length(coordinates)) {
                dataset <- datasets[which(datasets$AP == coordinates[q] & 
                  datasets$animal == animal.index), ]
                if (length(table(dataset$right.hemisphere)) > 
                  1) 
                  if (table(dataset$right.hemisphere)[2] < table(dataset$right.hemisphere)[1]) {
                    dataset$ML <- -dataset$ML
                  }
                coordinate <- unique(dataset$AP)
                k <- which(abs(coordinate - atlasIndex$mm.from.bregma) == 
                  min(abs(coordinate - atlasIndex$mm.from.bregma)))
                plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
                scale = 0.9579832
                if(length(dataset$image)>0){
                        bregmaX = 5640
                        bregmaY = 200
                    }else{
                      bregmaX = 5640
                      bregmaY = 200
                    }
                ventricles <- c(which(EPSatlas$plate.info[[k]]$style == 
                  "#aaaaaa"))
                if (atlasIndex$plate.id[k] %in% c(100960309, 
                  100960312, 100960316, 100960320)) {
                  fibertracts <- c(which(EPSatlas$plate.info[[k]]$style == 
                    "#cccccc"))
                  fibertracts <- fibertracts[which.min(sapply(fibertracts, 
                    function(x) {
                      min(EPSatlas$plates[[k]][[x]]@paths$path@y)
                    }))]
                  ventricles <- append(ventricles, fibertracts)
                }
                if (title) {
                  main.title <- paste(animal.index, "\n bregma: ", 
                    round(coordinate, 2), "mm", "\n image: ", 
                    unique(dataset$image))
                }
                else {
                  main.title <- ""
                }
                if (save.plots&device) {
                  quartz(width = dev.size[1], height = dev.size[1])
                }
                if (mm.grid) {
                }
                else {
                  par(mar = c(0, 0, 0, 0))
                }
                xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 
                  97440/2
                plot(EPSatlas$plates[[k]][[1]]@paths$path@x, 
                  EPSatlas$plates[[k]][[1]]@paths$path@y, col = 0, 
                  xlim = c(0, 97440), ylim = c(0, 68234.56), 
                  axes = F, ylab = "", xlab = "", asp = 1, main = "")
                if (mm.grid) {
                  abline(h = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y) + 
                    6 * (97440/456/25 * 1000), max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                    10 * (97440/456/25 * 1000), by = -(97440/456) * 
                    (1000/25)), col = "lightblue")
                  abline(v = seq(6040, 97440, by = (97440/456) * 
                    (1000/25)), col = "lightblue")
                  mtext("dorsoventral (mm)", 2, 3)
                  mtext("mediolateral (mm)", 1, 3)
                }
                polygon(EPSatlas$plates[[k]][[1]]@paths$path@x - 
                  xmin, EPSatlas$plates[[k]][[1]]@paths$path@y, 
                  col = gray(0.95), border = "black")
                polygon(-(EPSatlas$plates[[k]][[1]]@paths$path@x - 
                  xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[1]]@paths$path@y, 
                  col = gray(0.95), border = "black")
                if (title & (!mm.grid)) {
                  mtext(main.title, 1, -5.5, cex = 0.5, font = 2)
                }
                else {
                  title(main.title)
                }
                
                if (mm.grid) {
                  box()
                  axis(2, at = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y), 
                    max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                      6 * (97440/456/25 * 1000), by = -(97440/456) * 
                      (1000/25)/10), labels = FALSE, las = 1, 
                    col = "orange", tck = -0.0125)
                  axis(2, at = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y), 
                    max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                      6 * (97440/456/25 * 1000), by = -(97440/456) * 
                      (1000/25)/2), labels = FALSE, las = 1, 
                    col = "darkblue", tck = -0.025)
                  axis(2, at = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y), 
                    max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                      6 * (97440/456/25 * 1000), by = -(97440/456) * 
                      (1000/25)), labels = c(0:-6), las = 1)
                  axis(1, at = seq(6040, 91514, by = (97440/456) * 
                    (1000/25)/10), labels = FALSE, las = 1, col = "orange", 
                    tck = -0.0125)
                  axis(1, at = seq(6040, 91514, by = (97440/456) * 
                    (1000/25)/2), labels = FALSE, las = 1, col = "darkblue", 
                    tck = -0.025)
                  axis(1, at = seq(6040, 91514, by = (97440/456) * 
                    (1000/25)), labels = c(-5:5))
                }
                        if(scale.bar){
            right.corner<-c(82966.32, 91513.68)    
                polygon(c(right.corner[1], right.corner[1], right.corner[2], right.corner[2]),  c(0, 0-1000, 0-1000, 0), col='black' )
        }
                numPaths <- EPSatlas$plates[[k]]@summary@numPaths
                if(region.colors){
                    style<-as.character(EPSatlas$plate.info[[k]]$style)
                }else{
                    style<-rep(gray(0.95), numPaths)
                }
                lapply(1:numPaths, function(N) {
                  polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                    xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                    col = style[N], border = "black", lwd = 1, 
                    lty = 3)
                  polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                    xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                    col = style[N], border = "black", lwd = 1, 
                    lty = 3)
                })
                if (coordinate < 1.3) {
                  for (i in ventricles) {
                    polygon(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                      xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                      col = "black", border = "black", lwd = 4)
                    polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                      xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                      col = "black", border = "black", lwd = 4)
                  }
                }
                fibertracts <- c(which(EPSatlas$plate.info[[k]]$style == 
                  "#cccccc"))
                for (i in fibertracts) {
                  polygon(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                    xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                    col = gray(0.85), border = "black", lwd = 1)
                  polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                    xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                    col = gray(0.85), border = "black", lwd = 1)
                }
                if (!is.null(dataset)) {
                     if(length(dataset$image)!=0){
                     if(pch%in%c(21:25)){
                        points((dataset$ML * 1000) * 97440/456/25+ (97440/2), (8210 + dataset$DV * 1000 - bregmaY) * 97440/456/25, pch = pch, bg = as.character(dataset$color), col=col, cex = cex)

                    }
                    if(!(pch%in%c(21:25))){
                        points((dataset$ML * 1000) * 97440/456/25+ (97440/2), (8210 + dataset$DV * 1000 - bregmaY) * 97440/456/25, pch = pch, col = as.character(dataset$color), cex = cex)

                    }
                     }else{
                    if(pch%in%c(21:25)){
                  points((dataset$ML * scale * 1000 + bregmaX) * 
                    97440/456/25 + 1748.92, (8210 + dataset$DV * 
                    scale * 1000 - bregmaY) * 97440/456/25, pch = pch, 
                    bg = as.character(dataset$color), col=col, cex = cex)
                    }
                    if(!(pch%in%c(21:25))){
                    points((dataset$ML * scale * 1000 + bregmaX) * 
                    97440/456/25 + 1748.92, (8210 + dataset$DV * 
                    scale * 1000 - bregmaY) * 97440/456/25, pch = pch, 
                    col = as.character(dataset$color), cex = cex)
                    }
                    
                    }
                }
                if (save.plots) {
                  filename <- paste(animal.index, formatC(q, 
                    digits = 3, flag = "0"), round(unique(coordinate), 
                    3), ".pdf", sep = "_")
                  quartz.save(filename, type = "pdf")
                  dev.off()
                }
            }
        }
    }
    else {
        k <- which(abs(coordinate - atlasIndex$mm.from.bregma) == 
            min(abs(coordinate - atlasIndex$mm.from.bregma)))
        plate <- atlasIndex$plate.id[which(abs(coordinate - atlasIndex$mm.from.bregma) == 
            min(abs(coordinate - atlasIndex$mm.from.bregma)))]
        scale = 0.9579832
        bregmaX = 5640
        bregmaY = 0
        ventricles <- c(which(EPSatlas$plate.info[[k]]$style == 
            "#aaaaaa"))
        if (atlasIndex$plate.id[k] %in% c(100960309, 100960312, 
            100960316, 100960320)) {
            fibertracts <- c(which(EPSatlas$plate.info[[k]]$style == 
                "#cccccc"))
            fibertracts <- fibertracts[which.min(sapply(fibertracts, 
                function(x) {
                  min(EPSatlas$plates[[k]][[x]]@paths$path@y)
                }))]
            ventricles <- append(ventricles, fibertracts)
        }
        main.title <- ""
        if (save.plots) {
            quartz(width = dev.size[1], height = dev.size[1])
        }
        if (mm.grid) {
        }
        else {
            par(mar = c(0, 0, 0, 0))
        }
        xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 
            97440/2
        plot(EPSatlas$plates[[k]][[1]]@paths$path@x, EPSatlas$plates[[k]][[1]]@paths$path@y, 
            col = 0, xlim = c(0, 97440), ylim = c(0, 68234.56), 
            axes = F, ylab = "", xlab = "", asp = 1, main = "")
        if (mm.grid) {
            abline(h = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y) + 
                6 * (97440/456/25 * 1000), max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                10 * (97440/456/25 * 1000), by = -(97440/456) * 
                (1000/25)), col = "lightblue")
            abline(v = seq(6040, 97440, by = (97440/456) * (1000/25)), 
                col = "lightblue")
            mtext("dorsoventral (mm)", 2, 3)
            mtext("mediolateral (mm)", 1, 3)
        }
        polygon(EPSatlas$plates[[k]][[1]]@paths$path@x - xmin, 
            EPSatlas$plates[[k]][[1]]@paths$path@y, col = gray(0.95), 
            border = "black")
        polygon(-(EPSatlas$plates[[k]][[1]]@paths$path@x - xmin - 
            97440/2) + 97440/2, EPSatlas$plates[[k]][[1]]@paths$path@y, 
            col = gray(0.95), border = "black")
        if (title & (!mm.grid)) {
            mtext(main.title, 1, -5.5, cex = 0.5, font = 2)
        }
        else {
            title(main.title)
        }

        if (mm.grid) {
            box()
            axis(2, at = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y), 
                max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                  6 * (97440/456/25 * 1000), by = -(97440/456) * 
                  (1000/25)/10), labels = FALSE, las = 1, col = "orange", 
                tck = -0.0125)
            axis(2, at = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y), 
                max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                  6 * (97440/456/25 * 1000), by = -(97440/456) * 
                  (1000/25)/2), labels = FALSE, las = 1, col = "darkblue", 
                tck = -0.025)
            axis(2, at = seq(max(EPSatlas$plates[[k]][[1]]@paths$path@y), 
                max(EPSatlas$plates[[k]][[1]]@paths$path@y) - 
                  6 * (97440/456/25 * 1000), by = -(97440/456) * 
                  (1000/25)), labels = c(0:-6), las = 1)
            axis(1, at = seq(6040, 91514, by = (97440/456) * 
                (1000/25)/10), labels = FALSE, las = 1, col = "orange", 
                tck = -0.0125)
            axis(1, at = seq(6040, 91514, by = (97440/456) * 
                (1000/25)/2), labels = FALSE, las = 1, col = "darkblue", 
                tck = -0.025)
            axis(1, at = seq(6040, 91514, by = (97440/456) * 
                (1000/25)), labels = c(-5:5))
        }
        if(scale.bar){
            right.corner<-c(82966.32, 91513.68)    
            polygon(c(right.corner[1], right.corner[1], right.corner[2], right.corner[2]),  c(0, 0-1000, 0-1000, 0), col='black' )
        }
        numPaths <- EPSatlas$plates[[k]]@summary@numPaths
        if(region.colors){
            style<-as.character(EPSatlas$plate.info[[k]]$style)
        }else{
            style<-rep(gray(0.95), numPaths)
        }
        lapply(1:numPaths, function(N) {
            polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                col = style[N], border = "black", lwd = 1, 
                lty = 3)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                col = style[N], border = "black", lwd = 1, 
                lty = 3)
        })
        if (coordinate < 1.3) {
            for (i in ventricles) {
                polygon(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                  xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                  col = "black", border = "black", lwd = 4)
                polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                  xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                  col = "black", border = "black", lwd = 4)
            }
        }
        fibertracts <- c(which(EPSatlas$plate.info[[k]]$style == 
            "#cccccc"))
        for (i in fibertracts) {
            polygon(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                col = gray(0.85), border = "black", lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[i]]@paths$path@y, 
                col = gray(0.85), border = "black", lwd = 1)
        }
    }
    
    
}







paxTOallen<-function(paxinos){
 	round(214+(20-(paxinos*1000))/25)
}

glassbrain<-function(dataset, high.res=FALSE, dim=c(720,1080), device=TRUE, col='region', cex=0.5, hemisphere='right', spheres=FALSE, alpha=1, laterality=TRUE, plane='coronal'){
    if(sum(dataset$color=='#000000')>0){
	   dataset<-dataset[-which(dataset$color=='#000000'),]
    }
	if(device){
		open3d(windowRect = c(0,  0, 1280, 720))
	}
	par3d(persp)
	if(high.res){
		drawScene.rgl(list(VOLUME))
	}else{
		drawScene.rgl(list(VOLUMESMALL))
	}

	if(col=='region'){
		color<-as.character(dataset$color)
	}else{
		color<-col
	}
	if(hemisphere=='right'){
		hemisphere<-2
	}else{
		hemisphere<-1
	}

    if(plane=='coronal'){
        smp.AP<- rnorm(length(dataset$AP), 0,(320/9.75)*0.2 )
            smp.ML<-0
        }else{
            smp.ML<- rnorm(length(dataset$AP), 0,(320/9.75)*0.2 )
            smp.AP<-0
        }

    if(laterality){
	if(length(unique(dataset$AP))>1){
		laterality<-table(dataset$AP, dataset$right.hemisphere)
		for(i in 1:nrow(laterality)){
			if(hemisphere=='right'){
					if(laterality[i,1]>laterality[i,2]){
						 dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])]<- -dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])]
					} 
				}else{
					if(laterality[i,2]>laterality[i,1]){
						 dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])]<- -dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])]
					} 
				}
			
		}
	}
    }

    if(spheres){
    spheres3d(paxTOallen(dataset$AP)-530/2+smp.AP, -dataset$DV*1000/25-320/2, dataset$ML*1000/25+smp.ML, col=color, radius=cex, alpha=alpha )

        }else{
                points3d(paxTOallen(dataset$AP)-530/2+smp.AP, (-dataset$DV*1000/25*0.95)-320/2, dataset$ML*1000/25+smp.ML, col=color, size=cex, alpha=alpha )

        }

	
}

cortical.sideview<-function(dataset, right.hemisphere = TRUE, labels = TRUE,  col = NULL, pch = 16, pt.bg = NULL){
  hemisphere<-c('Left hemisphere', 'Right hemisphere')
  xlim<-rbind(c(5,-7), c(-7, 5))
  plot(region[,1:2], xlim=xlim[right.hemisphere + 1,], asp=1, ylim=c(-7,0), col=0, axes=F, ylab='Dorso-ventral [mm]', xlab='Anterior-posterior [mm]')
  mtext(hemisphere[right.hemisphere + 1], 3, 1, font=2)
  graylight<-0.9
  lapply(5:-8, function(x)abline(v=x, col='lightblue'))
  lapply(1:-8, function(x)abline(h=x, col='lightblue'))
  
  polygon(region[which(sideview$acronym=='brainstem'),1:2], col=gray(graylight))
  polygon(region[which(sideview$acronym=='cerebellum'),1:2], col=gray(graylight))
  polygon(region[which(sideview$acronym=='telencephalon'),1:2], col=gray(graylight))
  
  cortical<-unique(region$acronym[which(nchar(sideview$acronym)<=5)] )
  
  lapply(cortical, function(x){polygon(sideview[which(sideview$acronym==x),1:2], lty=3, col=gray(graylight+0.05))} )
  
  if(!is.null(col)){
    color<-col
  }else{
    color<-dataset$color
  }
    
  
  if(labels)
    lapply(cortical, function(x){coord<-apply(region[which(sideview$acronym==x),1:2], 2, median);text(coord[1], coord[2], x, cex=0.85, col=rgb(0.3,0,0.3))} )
  
  axis(2, at=c(0:-7), las=1)
  axis(1, at=c(5:-7))
}

cortical.plot<-function(dataset, type = c('left', 'right', 'top'), labels = TRUE,  col = NULL, pch = 16, pt.bg = NULL){
  #plot only both hemispheres
  if(c('left', 'right')%in%type){
    quartz(width=8.330275*2, height= 5.201835)
    par(mar=c(4,4,2,1), mfrow=c(1,2))
    cortical.sideview(right.hemisphere = FALSE, labels = labels)
    cortical.sideview(right.hemisphere = TRUE, labels = labels)
  }
  #plot all
  
  #plot only left or right
  if(xor('left'%in%type, 'right'%in%type)){
    cortical.sideview(right.hemisphere = (type == 'right'), labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
  }
  
  #plot  top
  if('top'%in%type){
    quartz(width=8.330275*2, height= 5.201835)
    par(mar=c(4,4,2,1), mfrow=c(1,2))
    cortical.sideview(right.hemisphere = FALSE, labels = labels)
    cortical.sideview(right.hemisphere = TRUE, labels = labels)
  }
  
}


bargraph <- function(dataset, device=TRUE, region.lab='Input region:') {
    print("bargraph() is deprecated! Use dot.plot() instead.")
    dot.plot(dataset, device=TRUE, region.lab='Input region:')
}

dot.plot<-function (dataset, device = TRUE, region.lab = "Input region:", xlab= 'Cell count') 
{
  if(ncol(dataset)<3){
    counts<-dataset
  }else{
    counts <- table(as.character(dataset$acronym), dataset$right.hemisphere)
  }
  hemisphere.to.sort <- which.max(colSums(counts*is.finite(counts),na.rm=TRUE ))
  counts <- counts[order(counts[, hemisphere.to.sort], decreasing = TRUE), ]
  counts <- log10(counts)

    if (device) {
      quartz(width = 7.036585, height = 0.2099039 * nrow(counts))
    }
    layout(matrix(c(1, 1, 1, 2, 2, 2, 2), nrow = 1))
    par(mar = c(4, 0, 4, 0))
    plot(rep(2.5, nrow(counts)), nrow(counts):1, col = 0, 
         axes = F, ylim = c(0.5, nrow(counts) + 0.5), ylab = "", 
         xlab = "", xlim = c(1, 5))
    mtext(region.lab, 3, cex = 0.9)
    for (i in 1:nrow(counts)) {
      regioncolor <- color.from.acronym(row.names(counts)[i])
      regioncolor <- adjustcolor(regioncolor, alpha.f = 0.1)
      y.lab <- (nrow(counts) + 1) - i
      polygon(c(1, 5, 5, 1), c(y.lab - 0.5, y.lab - 0.5, 
                               y.lab + 0.5, y.lab + 0.5), col = regioncolor, 
              border = FALSE)
      text(3, y.lab, name.from.acronym(row.names(counts)[i]), 
           cex = 0.9)
    }
    par(xpd = TRUE)
    legend(3, -0.75, c("Left", "Right"), pch = c(21), pt.bg = c("white", 
                                                                gray(0.2)), title = "Hemisphere:", bg = "white", 
           horiz = TRUE, cex = 1.3, xjust = 0.5)
    par(xpd = FALSE)
    par(mar = c(4, 4, 4, 6))
    zeros <- min(is.finite(counts)) - 1
    counts[!is.finite(counts)] <- zeros
    x.range <- ceiling(range(counts[is.finite(counts)]))
    plot(apply(counts, 1, max), nrow(counts):1 - 0.125, pch = 21, 
         bg = "white", ylim = c(0.5, nrow(counts) + 0.5), 
         xlim = x.range, xlab = "", axes = F, ylab = "", col = 0)
    for (i in 1:nrow(counts)) {
      regioncolor <- color.from.acronym(row.names(counts)[i])
      regioncolor <- adjustcolor(regioncolor, alpha.f = 0.05)
      y.lab <- (nrow(counts) + 1) - i
      polygon(c(x.range, rev(x.range)), c(y.lab - 0.5, 
                                          y.lab - 0.5, y.lab + 0.5, y.lab + 0.5), col = regioncolor, 
              border = FALSE)
    }
    log.range <- 10^seq(x.range[1], x.range[2])
    axis(1, at = seq(x.range[1], x.range[2]), las = 1, labels = c(0, 
                                                                  log.range[-1]))
    axis(3, at = seq(x.range[1], x.range[2]), las = 1, labels = c(0, 
                                                                  log.range[-1]))
    log.range <- unlist(lapply(1:(length(log.range) - 1), 
                               function(x) {
                                 seq(log.range[x], log.range[x + 1], by = log.range[x])
                               }))
    axis(1, at = log10(log.range), labels = FALSE)
    axis(3, at = log10(log.range), labels = FALSE)
    axis(2, at = nrow(counts):1, labels = row.names(counts), 
         las = 1)
    axis(4, at = nrow(counts):1, labels = row.names(counts), 
         las = 1)
    abline(h = 1:nrow(counts), lty = 2, col = "gray")
    abline(v = log10(log.range), col = "lightblue")
    lapply(1:nrow(counts), function(x) {
      lines(counts[x, ], rep(nrow(counts) - x + 1, 2), 
            lwd = 2)
    })
    points(counts[, 2], nrow(counts):1, 
           pch = 21, bg = gray(0.2), cex = 1.2)
    points(counts[, 1], nrow(counts):1, 
           pch = 21, bg = "white", cex = 1.2)
    box()
    par(xpd = TRUE)
    polygon(c(x.range[1] + 0.25, x.range[1] + 0.5, x.range[1] + 
                0.5, x.range[1] + 0.25), c(-2, -2, nrow(counts) + 
                                             3, nrow(counts) + 3), col = "white", border = "white")
    par(xpd = FALSE)
    abline(v = c(x.range[1] + 0.25, x.range[1] + 0.5))
    mtext(xlab, 3, 2.2, cex = 0.8)
    mtext(xlab, 1, 2.2, cex = 0.8)

}


plot.outlines<-function(regi, plot=FALSE){
scale.factor<-regi$transformationgrid$width/dim(regi$transformationgrid$mx)[2]
if(plot)
plot(c(regi$atlas$outlines[[1]]$xrT*scale.factor, regi$atlas$outlines[[1]]$xlT*scale.factor), c(regi$atlas$outlines[[1]]$yrT*scale.factor, regi$atlas$outlines[[1]]$ylT*scale.factor), type='l', lwd=2, xlim=c(regi$transformationgrid$width, 0), ylim=c(regi$transformationgrid$height, 0), axes=FALSE, asp=1, xlab='', ylab='')

for(i in 1:regi$atlas$numRegions){
  if(regi$atlas$col[i]=='#cccccc'){
    bg<-'#cccccc'
    lwd<-2
  }else{
    bg<-gray(0.95)
    lwd<-1
  }
  polygon(regi$atlas$outlines[[i]]$xlT*scale.factor, regi$atlas$outlines[[i]]$ylT*scale.factor, col=bg, lwd=lwd)
  polygon(regi$atlas$outlines[[i]]$xrT*scale.factor, regi$atlas$outlines[[i]]$yrT*scale.factor, col=bg, lwd=lwd)
}
}


#' Cortical topview plot
#'
#' Plots a schematic brain from the top showing only parts that are easily accessible with cranial windows without removing parts of the temporalis muscle.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param dv.cut numeric, value dorso-ventral plane where everything above dv.cut will be plotted. Default is -2.8 mm.
#' @param show.sections boolean, if lines should be drawn for each unique AP section. Default is TRUE.
#' @param labels if labels should be plotted with acronyms for each cortical region.
#' @param col color code or name.
#' @param pch plotting ‘character’, i.e., symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols.
#' @param pt.bg background (fill) color for the open plot symbols given by pch = 21:25.
#' @examples
#' cortical.topview(dataset, right.hemisphere = FALSE, labels = FALSE)
cortical.topview<-function(dataset, dv.cut = -2.8, show.sections = TRUE, labels = TRUE,  col = NULL, pch = 16, cex = 0.5, pt.bg = NULL){

brain.color<-rgb(0.25, 0.25, 0.25, 0.1)

par(pty='s')
plot(0,0, col=0, ylab='Anterior-Posterior (mm)', xlab='Medial-Lateral (mm)', xlim=c(-5,5), axes=F, ylim=c(-5,3.7))


abline(v=c(-5:5), col='lightblue')
abline(h=c(-5:5), col='lightblue')


polygon(topview$outline[,1]/1.95, topview$outline[,2]/1.95, lwd=2, col= brain.color )
polygon(-(topview$outline[,1]/1.95), topview$outline[,2]/1.95, lwd=2, col= brain.color)


axis(1, at=seq(-4,4, by=0.1), labels=F, col='gray', tck=-0.018)
axis(2, at=seq(-5,5, by=0.1), labels=F, col='gray', tck=-0.018)
axis(1, at=-4:4, labels=F)
axis(2, at=c(-5:5), labels=F)

axis(1, at=seq(-4,4, by=2))
axis(2, at=seq(-4,4, by=2))

abline(h=0)
points(0,0,pch=23, bg='white', cex=1.3, lwd=2)  
axis(4, at=0, labels='Bregma')
abline(h=-3.8)
points(0,-3.8,pch=23, bg='lightblue', cex=1.3, lwd=2)   
axis(4, at=-3.8, labels='lambda')

for(i in unique(topview$regions[,3])){
  polygon((topview$regions[which(topview$regions[,3]==i),1])/1.92, (topview$regions[which(topview$regions[,3]==i),2])/1.98, lty=3, border='black')
  polygon(-(topview$regions[which(topview$regions[,3]==i),1])/1.92, (topview$regions[which(topview$regions[,3]==i),2])/1.98, border='black', lty=3)
}

isocortex <- get.acronym.child(get.acronym.child(get.acronym.child("Isocortex")))
isocortex <- append(isocortex, get.acronym.child(get.acronym.child("Isocortex")))
isocortex <- na.omit(isocortex)
isocortex <- dataset$acronym %in% isocortex
isocortex <- (isocortex & (abs(dataset$DV) > dv.cut) )

if(!is.null(col)){
  color<-col
}else{
  color<-dataset$color[isocortex]
}

if(show.sections)
    abline(h = unique(dataset$AP), col="gray")

points(dataset$ML[isocortex], dataset$AP[isocortex], pch = pch, cex = cex, col = color, pt.bg = pt.bg)  

}

#' Cortical sideview plot
#'
#' Plots a schematic brain from the side with cortical regions outlined.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param right.hemisphere boolean, if true then right side of the hemisphere will be plotted, if false left hemisphere will be plotted.
#' @param ml.cut numeric, value medio-lateral plane where everything more medial than ml.cut will be excluded. Default is 2.0 mm. 
#' @param show.sections boolean, if lines should be drawn for each unique AP section. Default is TRUE.
#' @param labels if labels should be plotted with acronyms for each cortical region.
#' @param col color code or name.
#' @param pch plotting ‘character’, i.e., symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols.
#' @param pt.bg background (fill) color for the open plot symbols given by pch = 21:25.
#' @examples
#' cortical.sideview(dataset, right.hemisphere = FALSE, labels = FALSE)
cortical.sideview<-function(dataset, right.hemisphere = TRUE, ml.cut = 2.0, show.sections = TRUE, labels = TRUE,  col = NULL, pch = 16, pt.bg = NULL){
  hemisphere<-c('Left hemisphere', 'Right hemisphere')
  xlim<-rbind(c(5,-7), c(-7, 5))
  plot(region[,1:2], xlim=xlim[right.hemisphere + 1,], asp=1, ylim=c(-7,0), col=0, axes=F, ylab='Dorso-ventral [mm]', xlab='Anterior-posterior [mm]')
  mtext(hemisphere[right.hemisphere + 1], 3, 1, font=2)
  graylight<-0.9
  lapply(5:-8, function(x)abline(v=x, col='lightblue'))
  lapply(1:-8, function(x)abline(h=x, col='lightblue'))
  
  polygon(region[which(sideview$acronym=='brainstem'),1:2], col=gray(graylight))
  polygon(region[which(sideview$acronym=='cerebellum'),1:2], col=gray(graylight))
  polygon(region[which(sideview$acronym=='telencephalon'),1:2], col=gray(graylight))
  
  cortical<-unique(region$acronym[which(nchar(sideview$acronym)<=5)] )
  
  lapply(cortical, function(x){polygon(sideview[which(sideview$acronym==x),1:2], lty=3, col=gray(graylight+0.05))} )
  

  isocortex <- get.acronym.child(get.acronym.child(get.acronym.child("Isocortex")))
  isocortex <- append(isocortex, get.acronym.child(get.acronym.child("Isocortex")))
  isocortex <- na.omit(isocortex)
  isocortex <- dataset$acronym %in% isocortex
  isocortex <- (isocortex & (abs(dataset$ML) > ml.cut) )

  if(!is.null(col)){
    color<-col
  }else{
    color<-dataset$color[isocortex & dataset$right.hemisphere == right.hemisphere]
  }
  
  if(show.sections)
    abline(v = unique(dataset$AP), col="gray")
  
  points(dataset$AP[isocortex & dataset$right.hemisphere == right.hemisphere], dataset$DV[isocortex & dataset$right.hemisphere == right.hemisphere], pch = pch, col = color, pt.bg = pt.bg)  
  
  if(labels)
    lapply(cortical, function(x){coord<-apply(region[which(sideview$acronym==x),1:2], 2, median);text(coord[1], coord[2], x, cex=0.85, col=rgb(0.3,0,0.3))} )
  
  axis(2, at=c(0:-7), las=1)
  axis(1, at=c(5:-7))
}


#' Cortical plot
#'
#' Plots schematic brain(s) of cells in cortex from different views.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param type character vector, what type of plots should be made. Default is c('left', 'right', 'top').
#' @param show.sections boolean, if lines should be drawn for each unique AP section. Default is TRUE.
#' @param ml.cut numeric, value medio-lateral plane where everything more medial than ml.cut will be excluded. Default is 2.0 mm. 
#' @param dv.cut numeric, value dorso-ventral plane where everything above dv.cut will be plotted. Default is -2.8 mm.
#' @param labels if labels should be plotted with acronyms for each cortical region.
#' @param col color code or name.
#' @param pch plotting ‘character’, i.e., symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols.
#' @param pt.bg background (fill) color for the open plot symbols given by pch = 21:25.
#' @examples
#' ccortical.plot(dataset, type = c('left', 'right', 'top') )
cortical.plot<-function(dataset, type = c('left', 'right', 'top'), show.sections = TRUE, labels = TRUE, ml.cut = 2.0, dv.cut = -2.8, col = NULL, cex = 0.5, pch = 16, pt.bg = NULL){
  if(all(c('left', 'right')%in%type)  & (!'top'%in%type)){
    quartz(width=8.330275*2, height= 5.201835)
    par(mar=c(4,4,2,1), mfrow=c(1,2))
    cortical.sideview(dataset, right.hemisphere = FALSE, show.sections = show.sections, ml.cut = ml.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
    cortical.sideview(dataset, right.hemisphere = TRUE, show.sections = show.sections,ml.cut = ml.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
  }
  
  #plot only both hemispheres
  if(all(c('left', 'right', 'top')%in%type)){
    quartz(width=14.691892, height= 4.864865)
    layout(matrix(c(1, 1, 1, 1, 2, 2, 2,2,2,3, 3, 3, 3,3), 1, 14, byrow = TRUE))
    par(mar=c(4,4,4,2), yaxs="r")
    cortical.topview(dot, show.sections = FALSE)
    par(mar=c(4,2,3,0), yaxs="i" )
    cortical.sideview(dataset, right.hemisphere = FALSE, show.sections = FALSE)
    par(mar=c(4,2,3,4), yaxs="i" )
    cortical.sideview(dataset, right.hemisphere = TRUE, show.sections = FALSE)
  }
  #plot  top
  if(all(type == 'top')){
    cortical.topview(dataset, show.sections = show.sections, dv.cut = dv.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
  }
  
  #plot only left or right
  if(xor('left'%in%type, 'right'%in%type)){
    cortical.sideview(dataset, right.hemisphere = (type == 'right'), show.sections = show.sections, ml.cut = ml.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
  }
  
}
