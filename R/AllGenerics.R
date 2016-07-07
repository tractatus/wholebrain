


acronym.from.id<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$id ==y))!=0){return(ontology$acronym[which(ontology$id ==y)])}else{return(NA)} }))
}

id.from.acronym <-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$acronym ==y))!=0){return(ontology$id[which(ontology$acronym ==y)])}else{return(NA)} }))
}

name.from.acronym<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$acronym ==y))!=0){return(name[which(ontology$acronym ==y)])}else{return(NA)} }))
}

name.from.id<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$id ==y))!=0){return(name[which(ontology$id ==y)])}else{return(NA)} }))
}

color.from.id<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$id ==y))!=0){return(ontology$allen.color[which(ontology$id ==y)])}else{return(NA)} }))
}

color.from.acronym<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$acronym ==y))!=0){return(as.character(ontology$allen.color[which(ontology$acronym ==y)]))}else{return(NA)} }))
}


check.progress<-function(barWidth, progress, stages, message, timing){
	cat("  [")
     	pos = round(barWidth * progress);
    	for(p in seq(0,barWidth-1)) {
        	if(p < pos){ 
        		cat("=")
   			}
        	else{if (p == pos){
        		cat(">")
    		}
    		else{
    			cat(" ")
    		}
    		}	
    	}
    	message.new<-paste("] ", round(progress * 100.0), "% | ",message," | time left: ", timing," |\r", sep='')
    	paste(message.new, paste(rep(' ', nchar(message.new)), collapse=' '), sep='')
    	cat(message.new)
    	Sys.sleep(.05)
    	progress <- progress + 1/(stages);
    	flush.console()
    	process.events()
    	return(progress)
}


stitch.animal<-function(folder, rotate=0, FFC=TRUE, web.map=TRUE, start.at=1, dont.run=NULL){
	all.section.folder<-list.dirs(folder, recursive=FALSE, full.names=FALSE)
	remove<-c(which(substr(all.section.folder, 1,6)%in%c('output','stitch')),which(substr(all.section.folder, 1,3)%in%c('FFC','Web')))
	if(length(remove)>0){
		all.section.folder<-all.section.folder[-remove]
	}
	if(!is.null(dont.run)){
		all.section.folder<-all.section.folder[-dont.run]
	}

	if(start.at>1){
		if(length(start.at)>1){
				all.section.folder<-all.section.folder[start.at[1]:start.at[2]]
			}else{
				all.section.folder<-all.section.folder[start.at:length(all.section.folder)]
			}
	}
	barWidth = 50;
   progress = 0.0;
   processing.steps<-length(all.section.folder)*(2+FFC+web.map)
   elapsed.time<-'?'
   elapsed.time.series<-numeric()

   if(length(rotate)==1){rotate<-rep(rotate, length(all.section.folder))}else{
   			rotate<-rotate[start.at:length(rotate)]
   }
   section.folder<-''
   images<-''
	for(i in all.section.folder){
		ptm <- proc.time()
		section.folder<-paste(folder,i, sep='/')
		progress<-check.progress(barWidth, progress, processing.steps, paste('Flat-field correction on:',i,'| image',which(all.section.folder==i),'of',length(all.section.folder)), elapsed.time )

		#get images 
		images<-get.images(section.folder)
		#order them
		index<-basename(images)
		index<-gsub("[A-z \\.\\(\\)]","",index)
		index<-gsub("[- \\.\\(\\)]","",index)
		index<-as.numeric(index)
		index<-order(index)
		images<-images[index]
		#Flat-field correction
		if(FFC){
			flat.field.correction(images, verbose=FALSE)
			progress<-check.progress(barWidth, progress, processing.steps, paste('Stitching:',i,'| image',which(all.section.folder==i),'of',length(all.section.folder)), elapsed.time )

			FFC.folder<-paste('FFC',basename(section.folder), sep='_')
			FFC.folder<-paste(folder, FFC.folder, sep='/')
			images<-get.images(FFC.folder)
			images<-images[index]
			stitched<-paste('stitched_FFC',basename(section.folder), sep='_')
			stitched<-paste(folder, stitched, sep='/')
		}else{
			stitched<-paste('stitched',basename(section.folder), sep='_')
			stitched<-paste(folder, stitched, sep='/')
		}		
		
		#stitch
		stitch(images, rotate=rotate[which(all.section.folder==i)], verbose=FALSE) 
		progress<-check.progress(barWidth, progress, processing.steps, paste('Web output:',i,'| image',which(all.section.folder==i),'of',length(all.section.folder)), elapsed.time )

		if(web.map){
			defaultwd<-getwd()
			setwd(folder)
			stitched<-get.images(stitched)
			makewebmap(stitched, verbose=FALSE)
			progress<-check.progress(barWidth, progress, processing.steps, paste('DONE:',i,'| image',which(all.section.folder==i),'of',length(all.section.folder)), elapsed.time )

			setwd(defaultwd)
		}
		elapsed.time.series<- append(elapsed.time.series , round( (proc.time()[3] - ptm[3])*(length(all.section.folder)-which(all.section.folder==i)) ) )
		elapsed.time<-mean(elapsed.time.series)
		elapsed.time<- format(.POSIXct(elapsed.time,tz="GMT"), "%H:%M:%S")
	}
	progress<-check.progress(barWidth, progress, processing.steps, ' FINISHED ', elapsed.time )

}


stitch.experiment<-function(folder){
	all.animals<-dir(folder)
	for(i in all.animals){
		section.folder<-paste(folder,i, sep='/')
		stitch.animal(section.folder)
	}
}


spreadsheet.animal<-function(folder, file, sep=','){
	filename<-list.dirs(folder, recursive=FALSE, full.names=FALSE)
	sectionname<-basename(filename)
	filename<-filename[substr(sectionname, 1, nchar('stitched_FFC'))=='stitched_FFC']
	sectionname<-sectionname[substr(sectionname, 1, nchar('stitched_FFC'))=='stitched_FFC']

	colNames<-c("filename", "section", "rotation", "Min", "Max", "soma.area", "eccentricity", "")

	data<-data.frame(matrix(rep(0, ncol(data)*length(sectionname)), ncol=length(colNames) ) )

	names(data)<-colNames

	data$filename<-filename
	data$section<-sectionname

	write.table(data, file=file, sep=sep, row.names=FALSE)
	return(data)
}
