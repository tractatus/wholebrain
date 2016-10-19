data(glasbrain, envir=environment())
data(EPSatlas, envir=environment())
data(atlasIndex, envir=environment())
#data(atlasOntology, envir=environment())
data(ontology, envir=environment())


schematic.plot<-function(dataset, title=TRUE, save.plot=FALSE, dev.size=c(5.4, 4.465)){
	

if(!save.plot){
	quartz(width=dev.size[1], height=dev.size[1])
}
if(length(which(dataset$color=='#000000'))>0){
dataset<-dataset[-which(dataset$color=='#000000'),]
}

datasets<-dataset


for(animal.index in unique(datasets$animal) ){
	
coordinates<-sort(unique(datasets$AP[which(datasets$animal==animal.index)]), decreasing=TRUE )	
for(q in 1:length(coordinates)){
dataset<-datasets[which(datasets$AP==coordinates[q] & datasets$animal==animal.index ),]


if( length(table(dataset$right.hemisphere)) > 1)
if(table(dataset$right.hemisphere)[2] < table(dataset$right.hemisphere)[1]){
	dataset$ML<- -dataset$ML
}


coordinate<-unique(dataset$AP)
#quartz(width=7.513513* 0.4171346 , height=4.540540* 0.4171346)
k<-which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))
plate<-atlasIndex$plate.id[which(abs(coordinate-atlasIndex$mm.from.bregma)==min(abs(coordinate-atlasIndex$mm.from.bregma)))]


scale = 0.9579832; 
bregmaX = 5640;
bregmaY = 0;#725
  
ventricles <-c(which(EPSatlas$plate.info[[k]]$style=='#aaaaaa')) #117 is optic chiasm och   which(EPSatlas$plate.info[[k]]$structure_id== 117)
if(atlasIndex$plate.id[k] %in%c(100960309, 100960312, 100960316, 100960320)){
  fibertracts<-c(which(EPSatlas$plate.info[[k]]$style=='#cccccc'))
  fibertracts <-fibertracts[which.min(sapply(fibertracts, function(x){min(EPSatlas$plates[[k]][[x]]@paths$path@y )}))]
  ventricles<-append(ventricles, fibertracts) 
}

if(title){
	main.title<-paste(animal.index,'\n bregma: ', round(coordinate,2), 'mm', '\n image: ', unique(dataset$image))
}else{
	main.title<-''
}

if(save.plots){
	quartz(width=dev.size[1], height=dev.size[1])
}
par(mar=c(0,0,0,0))
xmin<-min(EPSatlas$plates[[k]][[1]]@paths$path@x)-97440/2
plot(EPSatlas$plates[[k]][[1]]@paths$path@x, EPSatlas$plates[[k]][[1]]@paths$path@y, col=0, xlim=c(0,97440), ylim=c(0, 68234.56), axes=F, ylab='', xlab='', asp=1, main= '' )
polygon(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin, EPSatlas$plates[[k]][[1]]@paths$path@y, col=gray(0.95), border='black' )
polygon(-(EPSatlas$plates[[k]][[1]]@paths$path@x-xmin - 97440/2)+97440/2 , EPSatlas$plates[[k]][[1]]@paths$path@y, col=gray(0.95), border='black')	
mtext(main.title, 1,-5.5, cex=0.5, font=2)

numPaths<-EPSatlas$plates[[k]]@summary@numPaths

lapply(1:numPaths, function(N){polygon(EPSatlas$plates[[k]][[N]]@paths$path@x-xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, col=gray(0.95), border='black', lwd=1, lty=3);
 polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x-xmin- 97440/2)+97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, col=gray(0.95), border='black', lwd=1, lty=3)
 })
 
 if(coordinate <1.3){
for(i in ventricles){
polygon(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, col='black', border='black', lwd=4)
polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin- 97440/2)+97440/2, EPSatlas$plates[[k]][[i]]@paths$path@y, col='black', border='black', lwd=4)
}
}


  fibertracts<-c(which(EPSatlas$plate.info[[k]]$style=='#cccccc'))
for(i in fibertracts){
polygon(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin, EPSatlas$plates[[k]][[i]]@paths$path@y, col=gray(0.85), border='black', lwd=1)
polygon(-(EPSatlas$plates[[k]][[i]]@paths$path@x-xmin- 97440/2)+97440/2, EPSatlas$plates[[k]][[i]]@paths$path@y, col=gray(0.85), border='black', lwd=1)
}

points( (dataset$ML*scale*1000+bregmaX)*97440/11700+ 1748.92, (8210+dataset$DV*scale*1000-bregmaY)*97440/11700, pch=21, bg= as.character(dataset$color) , cex=0.5 )
if(save.plots){
filename <- paste(animal.index,formatC(q,digits=3,flag="0"), round(unique(coordinate),3),".pdf",sep="_")
quartz.save(filename, type = "pdf")
dev.off()
}

}
}


}

paxTOallen<-function(paxinos){
 	round(214+(20-(paxinos*1000))/25)
}

glassbrain<-function(dataset, high.res=FALSE, dim=c(720,1080), device=TRUE, col='region', cex=0.5, hemisphere='right'){
	dataset<-dataset[-which(dataset$color=='#000000'),]
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

	points3d(paxTOallen(dataset$AP)-530/2+rnorm(length(dataset$AP), 0,(320/9.75)*0.2 ), -dataset$DV*1000/25-320/2, dataset$ML*1000/25, col=color, size=cex )

	
}

