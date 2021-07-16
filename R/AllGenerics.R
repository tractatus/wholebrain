#data(connectivity, envir=environment())

update.wholebrain<-function(){
	detach('package:wholebrain', unload=TRUE)
	remove.packages("wholebrain")
	if( .Platform$OS.type=="windows" ){
		devtools::install_github("tractatus/wholebrain", INSTALL_opts=c("--no-multiarch"))
	}else{
		devtools::install_github("tractatus/wholebrain")
	}
}

acronym.from.id<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$id ==y))!=0){return(as.character(ontology$acronym[which(ontology$id ==y)]))}else{return(NA)} }))
}

id.from.acronym <-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$acronym ==y))!=0){return(ontology$id[which(ontology$acronym ==y)])}else{return(NA)} }))
}

name.from.acronym<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$acronym ==y))!=0){return(ontology$name[which(ontology$acronym ==y)])}else{return(NA)} }))
}

name.from.id<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$id ==y))!=0){return(ontology$name[which(ontology$id ==y)])}else{return(NA)} }))
}

color.from.id<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$id ==y))!=0){return(ontology$allen.color[which(ontology$id ==y)])}else{return(NA)} }))
}

color.from.acronym<-function(x){
	unlist(lapply(x, function(y){if(length(which(ontology$acronym ==y))!=0){return(as.character(ontology$allen.color[which(ontology$acronym ==y)]))}else{return(NA)} }))
}

get.acronym.parent<-function(x){
	ids<-unlist(lapply(x, function(y){if(length(which(ontology$acronym ==y))!=0){if(y=='root'){return('997')}else{return(ontology$parent[which(ontology$acronym ==y)])}}else{return(NA)} }))
	return(acronym.from.id(ids))
}

get.acronym.child<-function(x){
	ids<-unlist(lapply(x, function(y){if(length(which(ontology$parent == ontology$id[which(ontology$acronym == y)]))!=0){if(y=='root'){return('997')}else{return(ontology$id[which(ontology$parent == ontology$id[which(ontology$acronym == y)])])}}else{return(NA)} }))
	return(acronym.from.id(ids))
}

get.sup.structure<-function(x, matching.string=c('CTX','CNU','IB','MB','HB','grey','root','VS','fiber tracts')){
	
	if(x%in%matching.string){
		return(x)
	}
	tmp<-get.acronym.parent(x)
	if((tmp%in%matching.string)){
    	tmp<-x  
	}
	tmp2<-tmp
	while(!(tmp%in% matching.string) ){
    	tmp2<-tmp
    	tmp<-get.acronym.parent(tmp)
    }
	if(tmp=='root'|tmp=='grey'|tmp=='CH'){
		return(tmp2)
	}else{
		return(tmp)
	}
}

get.sub.structure<-function(x){
	tmp<-get.acronym.child(x)
	if(sum(is.na(tmp)/length(tmp))!=0){
    	tmp<-x
    	return(tmp)  
	}
	tmp2<-tmp
	for(i in tmp){
		tmp2<-append(tmp2, get.sub.structure(i) )
	}
	return(tmp2)
}

remove.whitematter<-function(dataset){
	dataset<-dataset[-which(dataset$acronym%in%c('grey', 'fibertracts')),]
	return(dataset)
}

get.common.regions<-function(dataset1, dataset2){
	regions.to.reduce<-names(table(dataset1$acronym))[which(!names(table(dataset1$acronym))%in%names(table(dataset2$acronym)))]

	new.names<-get.acronym.parent(regions.to.reduce)

}


agg.datasets<-function(folder){
	files<-dir(folder, full.names=TRUE)
	files<-files[tools::file_ext(files)=='RData']
	load(files[1])
	dataset.master<-dataset
	for(i in seq_along(files)[-1]){
		load(files[i])
		dataset.master<-rbind(dataset.master, dataset)
	}
	return(dataset.master)
}


get.vector.intensity<-function(input, x, y){
	file <- as.character(input)
    ## check for existence
    if(!file.exists(file)){
      stop(file, ", file not found")
      return(NULL)
  	}
    file <- path.expand(file)

	intensity<-.Call("getVectorIntensity", file, as.integer(x), as.integer(y))
	return(intensity)
}

#' Get region volume
#'
#' Providing an acronym, id or name this function returns the volume in cubic millimeters of the region. Option to only output for single hemisphere.
#' Returns a numeric vector of values in cubic millimeters.
#' @param regions vector of either acronyms, names or integer ids of regions that you want volume of.
#' @param bilateral boolean. If volume should be calculated for a single hemisphere or both. Default is TRUE.
#' @export
#' @examples
#' volume <- get.region.volume(dataset$acronym, bilteral = FALSE)
get.region.volume <- function(regions, bilateral = FALSE){
  volume<-ontology$volume[match(regions,ontology$acronym)]
  return(volume*(1-0.5*bilateral))
}

#' Get pixel intensity
#'
#' This function extracts the pixel intensity for a given set of x and y coordinates.
#' Returns a numeric vector of values.
#' @param input file path to a valid monochrome 16-bit TIFF.
#' @param x pixel x coordinate(s).
#' @param y pixel y coordinate(s).
#' @param type what type of intensity value default is 'SNR' set to 'intensity' to get raw pixel intensities.
#' @param roi .
#' @param background .
#' @keywords pixel intensity
#' @export
#' @examples
#' snr <- get.pixel.intensity(image, dataset$x, dataset$y, type='SNR', roi=10, background=40)
#' intensity <- get.pixel.intensity(image, dataset$x, dataset$y, type='intensity', roi=10)
get.pixel.intensity<-function(input, x, y, type='SNR', roi= 9,  background=40){
	file <- as.character(input)
    ## check for existence
    if(!file.exists(file)){
      stop(file, ", file not found")
      return(NULL)
  	}
    file <- path.expand(file)
    if(type=='SNR'){
    	type=2
    }else{
    	type=1
    }
	intensity<-.Call("getPixelIntensity", file, as.integer(x), as.integer(y), type, roi, background)
	return(intensity$intensity)
}

get.range<-function(input){
	file <- as.character(input)
    ## check for existence
    if(!file.exists(file)){
      stop(file, ", file not found")
      return(NULL)
  	}
    file <- path.expand(file)
 
	maxmin<-.Call("getMaxMin", file)
	return(maxmin)
}

imstats<-function(input){
	file <- as.character(input)
    ## check for existence
    if(!file.exists(file)){
      stop(file, ", file not found")
      return(NULL)
  	}
    file <- path.expand(file)
 
	maxmin<-.Call("getImgStats", file)
	return(maxmin)
}

#' Get multi channel
#'
#' Gets intensity or signal-to-noise measurement for segmented features in multiple channels.
#' @param seg a segmentation object.
#' @param images a character vector with file paths for images.
#' @param images a character vector with file paths for images.
#' @examples
#' dataset<-reflect.cells(dataset)
#get.multi.channel<-function(image = c('/Volumes/General/Daniel/C1-MAX_05_N1_2017-09-12_006_006.tif', '/Volumes/General/Daniel/C2-MAX_05_N1_2017-09-12_006_006.tif','/Volumes/General/Daniel/C3-MAX_05_N1_2017-09-12_006_006.tif',  '/Volumes/General/Daniel/C4-MAX_05_N1_2017-09-12_006_006.tif'), seg = seg, roi=9, background = 40){
#  mylist<-list()
  
#  for(i in )
#  assign(paste(), i) 
  
#  mylist$Cy5<-get.pixel.intensity(image[1], seg$soma$x, seg$soma$y, roi=roi, background = background)
#  mylist$TxRd<-get.pixel.intensity(image[2], seg$soma$x, seg$soma$y, roi=roi, background = background)
#  mylist$Cy3<-get.pixel.intensity(image[3], seg$soma$x, seg$soma$y, roi=roi, background = background)
#  mylist$FITC<-get.pixel.intensity(image[4], seg$soma$x, seg$soma$y, roi=roi, background = background)
#  mylist$x<-seg$soma$x
#  mylist$y<-seg$soma$y
#  basecall<-data.frame(do.call("cbind", mylist))
#  return(basecall)
#}



#' Flip cells to the ipsilateral side
#'
#' Flips the cells to the right hemisphere.
#' @param dataset a data frame generated by inspect.registration() or get.cell.ids().
#' @param right.hemisphere boolean, if TRUE then cells will be flipped to right side.
#' @examples
#' dataset<-reflect.cells(dataset)
reflect.cells<-function(dataset){
	if (length(table(dataset$right.hemisphere)) > 
                  1){ 
                  if (table(dataset$right.hemisphere)[2] < table(dataset$right.hemisphere)[1]) {
                    dataset$ML <- -dataset$ML
                  }
                  }
                  return(dataset)
}


#' Get cell counts for a specific Region of Interest (ROI).
#'
#' Get cell counts for a specific Region of Interest (ROI) in tidy format.
#' @param dataset a data frame generated by inspect.registration() or get.cell.ids().
#' @param rois a character vector with acronyms for brain regions you want.
#' @examples
#' #get cell position and identity
#' dataset<-get.cell.ids(regi, seg)
#' #extract roi counts
#' data<-roi.cell.count(dataset, rois=c('MO','TH'))
roi.cell.count<-function(dataset, rois=c('MO','TH')){
	cell.counts<-table(dataset$acronym, dataset$animal)
	regions<-get.sub.structure(rois[1])
	regions<-which(row.names(cell.counts)%in%regions)
	if(length(dim(cell.counts[regions,]))>1){
		roi.data<-colSums(cell.counts[regions,])
	}else{
		roi.data<-cell.counts[regions,]
	}
	for(i in rois[-1]){
		cell.counts<-table(dataset$acronym, dataset$animal)
		regions<-get.sub.structure(i)
		regions<-which(row.names(cell.counts)%in%regions)
		if(length(dim(cell.counts[regions,]))>1){
			roi.data<-rbind(roi.data, colSums(cell.counts[regions,]))
		}else{
			roi.data<-rbind(roi.data, cell.counts[regions,])
		}
		
	}
	
	roi.data <-data.frame(acronym=rep(rois, ncol(roi.data)), animal= rep(colnames(roi.data), each=nrow(roi.data)), cell.count = as.numeric(roi.data) )

	return(roi.data)	
}

summary.func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      err = se(x[[col]], na.rm=TRUE))
  }


#' Summary mean and dispersion
#'
#' Function to calculate the mean and the dispersion for each group
#' @param a data frame
#' @param varname the name of a column containing the variable to be summariezed
#' @param groupnamesvector of column names to be used as grouping variables
#' @examples
#' df <- data.summary(data, varname="cell.count", groupnames=c("genotype", "dose"))
data.summary <- function(data, varname, groupnames){
  
  data_sum<-ddply(data, groupnames, .fun=summary.func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

#' Add group variable
#'
#' Function to add a group variable such as genotype or experimental group based on animal ID.
#' @param data a data frame
#' @param subjectID a character vector 
#' @param group a character vector ordered accordin to subjectID which indicates what group.
#' @examples
#' df <- add.group(data, subjectID=c('R0057A317','R0057A319','R0057A322','R0057A335','R0057A341','R0057A342','R0057A345','R0057A346' 'R0057A347' 'R0057A349','R0057A350','R0057A355','R0057A358','R0057A360','R0057A361', 'R0057A365','R0057A371','R0057A375'), group=c("genotype", "dose"))
add.group <- function(data, subjectID, group){

  data$group<-as.character(data$animal)
  for(i in 1:length(subjectID)){
  	data$group[which(as.character(data$group)==subjectID[i])]<-group[i]
  }	

 return(data)
}

#' Convert to grayscale
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
rgb2gray<-function(input,verbose=TRUE, savefilename=TRUE, invert=TRUE, rotate=0){
    file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)
    
     if(savefilename==TRUE){
      pos <- regexpr("\\.([[:alnum:]]+)$", basename(file))	
      filename<-ifelse(pos > -1L, substring(basename(file), 1, pos - 1L), "")
      savefilename<-paste(getwd(),'/grayscale_',filename,'.tif', sep='')
     }

    a<-.Call("rgbTogray", file, as.integer(verbose), savefilename, as.integer(invert), rotate)
    
    return(savefilename)
}

invert.img<-function(input, verbose=TRUE, savefilename=TRUE){
    file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)
    
     if(savefilename==TRUE){
      pos <- regexpr("\\.([[:alnum:]]+)$", basename(file))	
      filename<-ifelse(pos > -1L, substring(basename(file), 1, pos - 1L), "")
      savefilename<-paste(getwd(),'/inverted_',filename,'.tif', sep='')
     }

    a<-.Call("invertImg", file, as.integer(verbose), savefilename)
    
    return(savefilename)
}

morph<-function(input, element = 0, size = 0, operator = 'tophat', save.as.8bit = TRUE, verbose = TRUE, savefilename=TRUE) {
	file <- as.character(input)
    ## check for existence
    if(!file.exists(file))
      stop(file, ", file not found")
    file <- path.expand(file)
    #"Operator:\n 0: Opening - 1: Closing \n 2: Gradient - 3: Top Hat \n 4: Black Hat"
    if(operator=='tophat'){
    	morphOperator<-3
    }
    if(operator=='open'){
    	morphOperator<-0
    }
    if(operator=='close'){
    	morphOperator<-1
    }
    if(operator=='gradient'){
    	morphOperator<-2
    }
    if(operator=='blackhat'){
    	morphOperator<-4
    }

    if(savefilename==TRUE){
      pos <- regexpr("\\.([[:alnum:]]+)$", basename(file))	
      filename<-ifelse(pos > -1L, substring(basename(file), 1, pos - 1L), "")
      savefilename<-paste(getwd(),'/',operator,'_',filename,'.tif', sep='')
     }

     a<-.Call("morphologyEx", input, element, size, morphOperator, as.integer(save.as.8bit), as.integer(verbose), savefilename)
     return(savefilename)
}

make.movie<-function(directory){
	
}

legacy<-function(dataset){

}

pax.to.allen<-function(paxinos){
 	round(214+(20-(paxinos*1000))/25)
}

read.brains<-function(filenames, animalID=NULL, type = '.RData', dataset.object = 'dataset'){

data<-read.table(filenames[1], sep=',', header=TRUE, fill=TRUE)
if(!is.null(animalID)){
data$animal<-animalID[1]
}
for(i in 1:length(filenames)){
	data.tmp<-read.table(filenames[i], sep=',', header=TRUE, fill=TRUE)
	if(!is.null(animalID)){
		data.tmp$animal<-animalID[i]
	}
	data<-rbind(data, data.tmp)
}
return(data)
}

read.brains2<-function(filenames, animalID=NULL, type = '.RData', dataset.object = 'dataset'){

if(type == '.RData'){
  read.section<-function(x){
    return(read.table(x, sep=',', header=TRUE, fill=TRUE))
  }
}
if(type == '.csv'){
  read.section<-function(x){
    return(read.table(x, sep=',', header=TRUE, fill=TRUE))
  }
}

data<-read.table(filenames[1], sep=',', header=TRUE, fill=TRUE)
if(!is.null(animalID)){
data$animal<-animalID[1]
}
for(i in 1:length(filenames)){
  data.tmp<-read.table(filenames[i], sep=',', header=TRUE, fill=TRUE)
  if(!is.null(animalID)){
    data.tmp$animal<-animalID[i]
  }
  data<-rbind(data, data.tmp)
}
return(data)
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


stitch.experiment<-function(folder, rotate=0, FFC=TRUE, web.map=TRUE){
	all.animals<-dir(folder)
	for(i in all.animals){
		section.folder<-paste(folder,i, sep='/')
		stitch.animal(section.folder, rotate = rotate, FFC= FFC, web.map=web.map)
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

#' Normalize volume
#'
#' Returns a count table of cell counts per region normalized by cubic millimeters.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param bilateral boolean, if true then cells will be normalized by total volume regardless of hemisphere.
#' @param aggregate compute volume for subregins based on parent region. Add which parent region you want to aggregate.
#' @examples
#'densities<-normalize.volume(dataset, aggregate = c('MO', 'SSp'))
normalize.volume  <- function(dataset, bilateral = TRUE, aggregate = NULL){
  if(!bilateral){
    counts<-table(dataset$acronym, dataset$right.hemisphere)
  }else{
    counts<-table(dataset$acronym)
  }
  
  parent.regions<-row.names(counts)
  
  if(!is.null(aggregate)){
    parent.regions<-row.names(counts)
    parent.regions[parent.regions%in%get.sub.structure(aggregate)]<-aggregate
  }
  
  volume<-ontology$volume[match(parent.regions,ontology$acronym)]
  
  while(any(volume==0)){
    parent.regions[volume==0]<-get.acronym.parent(get.acronym.parent(parent.regions[volume==0]))
    volume<-ontology$volume[match(parent.regions,ontology$acronym)]
    
  }
  
  
  counts<-sweep(counts, 1, volume*(1-0.5*bilateral), FUN = "/")
  
  return(counts)
  
}

#' Get index over cells in specific cortical layers
#'
#' Returns a two-dimensional data frame with the first column, layer, indicating which layer the cell is in and the second column, index, indicating which row in the original data frame the cell can be found.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param layer integer vector with unique layers of interest c(1:6) will give you all layers, default is c(1:6).
#' @examples
#'only.layers<-get.cortex.layer(dataset, layer = c(1:6))
#'dataset$acronym[only.layers$index] <- get.acronym.parent(dataset$acronym[only.layers$index])
#'dot.plot(dataset)'
get.cortex.layer<-function(dataset=dataset, layer=c(1:6)){
	layer[(layer==2)|(layer==3)]<-23
	remove.non.numbers<-abs(as.numeric(gsub("[A-z]|/","", dataset$acronym)))
	index<-which(remove.non.numbers%in%layer)
	index<-data.frame(layer= remove.non.numbers[index], index= index)
	return(index)
}

flag.outlier<-function(dataset){
	
}

#' Get Allen 25 micron voxel coordinates from stereotactic
#'
#' Returns integer of anterior-posterior voxel index.
#' @param paxinos AP form bregma.
pax.to.allen<-function(paxinos){
  round(214+(20-(paxinos*1000))/25)
}


inside.blobs<-function(dataset, blobs){
  if(class(dataset)=='list'){
  }else if(class(dataset)=='data.frame'){
    test<-(lapply(unique(blobs$soma$contour.ID), function(x)point.in.polygon(dataset$x, dataset$y, blobs$soma$contour.x[which(blobs$soma$contour.ID == x)], blobs$soma$contour.y[which(blobs$soma$contour.ID == x)])))
    return(apply(do.call("cbind", test), 1, sum)>0)
  }
}

get.region.from.coordinate<-function(AP, ML, DV){
    url<-paste0('http://mouse.brain-map.org/agea/data/P56/voxel_info?seed=', round(-981.1*AP+5244.1,0), '%2C', round(abs(DV)*1000,0), '%2C', round(1000*(ML+11.256/2), 0) )
    data<-data.frame(id = integer(), abbreviation = character(),  name = character(),  color = character()  )
    for(i in seq_along(url)){
        doc <- xmlParse(url[i])
        doc<-xmlToDataFrame(xpathApply(doc, "//voxelInfo/voxel/label"))
        data<-rbind(data, doc)
    }
    
    return(data)
}




get.projection.strength<-function(target = cbind(AP, ML, DV), source = cbind(AP, ML, DV), figure = 4){
    cat('Fetching from Allen atlas...\n')
    if(figure == 4){
        figX<-fig4
    }else{
        figX<-fig3
    }
    target.acronym<-get.region.from.coordinate(target[,1], target[,2], target[,3])$abbreviation
    cat('DONE\n')
    target.acronym<-as.character(target.acronym)
    
    simplify<-function(target.acronym, figX){
        need.to.fetch <- target.acronym %in% row.names(figX$ipsi)
        for(i in which(!need.to.fetch)){
            parent <- get.acronym.parent(target.acronym[i])
            if( parent %in% row.names(figX$ipsi) ){
                target.acronym[i] <- parent
            }
        }
        return(target.acronym)
    }
    
    INSERT_NA_ROW <- function(indice, tabla) {
        new_Row <- NA
        colName<-names(tabla)
        rowName<-row.names(tabla)
        long <- NROW(tabla)
        if(long>1){
            new_Data<- rbind(tabla[1:indice,], new_Row ,tabla[(indice + 1):(long),])
            new_Data<-t(new_Data)
            tabla<- rbind(new_Data[1:indice,], new_Row ,new_Data[(indice + 1):(long),])
            tabla<-data.frame(t(tabla))
        }else{
            new_Data<- tabla[1:(indice + 1),]
            new_Data<-t(new_Data)
            tabla<- rbind(new_Data[1:indice,], NA ,new_Data[(indice + 1):nrow(new_Data),])
            tabla<-data.frame(t(tabla))
        }
        
        names(tabla)[-(indice+1)]<-colName
        row.names(tabla)[-(indice+1)]<-rowName
        index<-which(row.names(tabla) == 'NA')
        if(length(index)>0){
            row.names(tabla)[index]<-paste('NA', seq_along(index), sep='.')
            names(tabla)[index]<-paste('NA', seq_along(index), sep='.')
        }
        
        return(tabla)
    }
    
    target.acronym<-simplify(target.acronym, figX)
    
    
    source.acronym<-get.region.from.coordinate(source[,1], source[,2], source[,3])$abbreviation
    source.acronym<-as.character(source.acronym)
    source.acronym<-simplify(source.acronym, figX)
    
    ipsi<-figX$ipsi[which(row.names(figX$ipsi) %in% source.acronym ), which(row.names(figX$ipsi) %in% target.acronym )]
    for(i in which(source.acronym == 'NA') - 1){
        ipsi<-INSERT_NA_ROW(i, ipsi)
    }
    
    
    contra<-figX$contra[which(row.names(figX$contra) %in% source.acronym ), which(row.names(figX$contra) %in% target.acronym )]
    for(i in which(source.acronym == 'NA') - 1){
        contra<-INSERT_NA_ROW(i, contra)
    }
    
    connections<-list(ipsi = ipsi, contra = contra)
    
    labels.dir <- expand.grid(row.names(connections$ipsi), names(connections$ipsi))
    
    labels.dir <- apply( labels.dir , 1 , paste , collapse = " -> " )
    
    projections <- data.frame(projection = labels.dir, weight = matrix(data.matrix(connections$ipsi), ncol=1), contra = FALSE )
    
    labels.dir <- expand.grid(row.names(connections$contra), names(connections$contra))
    
    labels.dir <- apply( labels.dir , 1 , paste , collapse = " -> " )
    
    projections.tmp <- data.frame(projection = labels.dir, weight = matrix(data.matrix(connections$ipsi), ncol=1), contra = TRUE )
    
    connections$projections <- rbind(projections, projections.tmp)
    
    return(connections)
}


get.projection<-function(injection = c(AP, DV, ML), AP = numeric(), DV = numeric(), ML = numeric(), display = FALSE ){
    
    matched.scale<- round(c(paxTOallen(injection[1]), -injection[2]*1000/25, injection[3]*1000/25+456/2) *c(13150/528, 7900/320, 11200/456), 0)
    url.download<-'http://connectivity.brain-map.org/projection/csv?criteria=service::mouse_connectivity_injection_coordinate[injection_structures$eq8,304325711][seed_point$eq%d,%d,%d][primary_structure_only$eqtrue][injection_distance_threshold$eq500]'
    cat("Downloading from allen...\n")
    allen.download<-read.table(url(sprintf(url.download, matched.scale[1], matched.scale[2], matched.scale[3])), sep=',', header=TRUE)
    allen.download<-allen.download[which(allen.download$distance<500),]
    experiment.id<-allen.download$id[which.max(allen.download$injection.volume)]
    
    file.to.download <- sprintf("http://api.brain-map.org/grid_data/download_file/%d??image=projection_density&resolution=100", experiment.id)
    file.to.open<-download.file(file.to.download, 'trashThis.nrrd', method = 'curl')
    
    scale.factor<-8
    src<-list(AP = integer(), DV= integer(), ML= integer())
    src$AP<- as.integer(paxTOallen(injection[1])/scale.factor)
    src$DV<- as.integer((-injection[2] * 1000/25)/scale.factor)
    src$ML<- as.integer( (injection[3]* 1000/25 + 456/2)/scale.factor )
    trgt<-list(AP = integer(), DV= integer(), ML= integer())
    
    trgt$AP <- as.integer(paxTOallen(AP)/scale.factor)
    trgt$DV <- as.integer(-(DV) * 1000/25/scale.factor ) #as.integer(-(DV-1.75) * 1000/25/scale.factor )
    trgt$ML <- as.integer( (ML* 1000/25 + 456/2)/scale.factor)
    
    VOL<-read.nrrd('trashThis.nrrd')
    VOL<-VOL>0+0
    VOL<-VOL[seq(1, dim(VOL)[1], by=2), seq(1, dim(VOL)[2], by=2), seq(1, dim(VOL)[3], by=2)]
    
    if(display){
        plot(c(0, dim(VOL)[3]), c(dim(VOL)[2], 0), ylim=c(40, 0), asp=1, axes=FALSE, ylab='', xlab='', type='n')
        rasterImage(VOL[src$AP, , ]>0, xleft = 0, ybottom = 40, xright = 57, ytop = 0)
        points(src$ML, src$DV, pch=16, col='orange')
        plot(c(0, dim(VOL)[3]), c(dim(VOL)[2], 0), ylim=c(40, 0), asp=1, axes=FALSE, ylab='', xlab='', type='n')
        rasterImage(VOL[trgt$AP[3], , ]>0, xleft = 0, ybottom = 40, xright = 57, ytop = 0)
        points(trgt$ML[1], trgt$DV[1], pch=16, col='purple')
    }
    
    
    paths<-list()
    for(i in seq_along(trgt$AP)){
        paths[[i]]<-searchPath(VOL, c(src$AP, src$DV, src$ML), c(trgt$AP[i], trgt$DV[i], trgt$ML[i]))
    }

    return(paths)
}
 
