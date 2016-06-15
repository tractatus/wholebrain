stitch.animal<-function(folder, rotate=0, FFC=TRUE, web.map=TRUE){
	all.section.folder<-dir(folder)
	for(i in all.section.folder){
		section.folder<-paste(folder,i, sep='/')
		#get images 
		images<-get.images(section.folder)
		#order them
		index<-basename(images)
		index<-gsub("[A-z \\.\\(\\)]","",index)
		index<-as.numeric(index)
		index<-order(index)
		images<-images[index]
		#Flat-field correction
		if(FFC){
			flat.field.correction(images)
			FFC.folder<-paste('FFC',basename(section.folder), sep='_')
			FFC.folder<-paste(folder, FFC.folder, sep='/')
			images<-get.images(FFC.folder)
			images<-images[index]
			stitched<-paste('stitched_FFC',basename(section.folder), sep='_')
			stitched<-paste(folder, stitched, sep='/')
		}else{
			stitched<-paste('stitched_',basename(section.folder), sep='_')
			stitched<-paste(folder, stitched, sep='/')
		}		
		
		#stitch
		stitch(images, rotate=rotate)
		if(web.map){
			defaultwd<-getwd()
			setwd(folder)
			stitched<-get.images(stitched)
			makewebmap(stitched)
			setwd(defaultwd)
		}
	}
}


stitch.experiment<-function(folder){
	all.animals<-dir(folder)
	for(i in all.animals){
		section.folder<-paste(folder,i, sep='/')
		stitch.animal(section.folder)
	}
}