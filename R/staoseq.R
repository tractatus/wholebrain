staoseq.prune.features<-function(features, width=26000, height){
 if(missing(height)){
 	height<-width
 }
 if(missing(width)){
 	width<-height
 }


poly<-list()
poly$x<-c(median(features$soma$x)+c(-width/2,+width/2,+width/2,-width/2))
poly$y<-c(median(features$soma$y)+c(-height/2,-height/2,+height/2,+height/2))

remove<-which(point.in.polygon(features$soma$x, features$soma$y, poly$x, poly$y)==0)

features$soma$x<-features$soma$x[-remove]
features$soma$y<-features$soma$y[-remove]
features$soma$intensity<-features$soma$area[-remove]
features$soma$area <-features$soma$area[-remove]


return(features)

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
plot(round(features$soma$x*(1/(0.5))), round(features$soma$y*(1/(0.5)), -2), asp=1, ylim=c(32768*2, 0), xlim=c(0,2*29184), axes=F, xlab='', ylab='', bg=rgb(1,0,0,0.2), pch=21)
points(cell.bodies$soma$x*2, cell.bodies$soma$y*2, pch=16, cex=0.25)
}

	return(list(features =features, cell.bodies =cell.bodies, filename=anatomyTIF))
}