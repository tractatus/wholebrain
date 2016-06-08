EPSILON<-0.000001

crossProduct<-function(a, b){
	return(a[1]*b[2]-b[1]*a[2])
}


#
# Checks if a Point is on a line
# @param a line (interpreted as line, although given as line
#                segment)
# @param b point
# @return <code>true</code> if point is on line, otherwise
#         <code>false</code>
#

isPointOnLine<-function(LineSegment.a, Point.b) {
    # Move the image, so that a.first is on (0|0)
    aTmp = rbind( c(0, 0),
    						  c(
            				LineSegment.a[2,1] - LineSegment.a[1,1], LineSegment.a[2,2] - LineSegment.a[1,2])
            			)
    bTmp = c(Point.b[1]  - LineSegment.a[1,1], Point.b[2] - LineSegment.a[1,2])
    r = crossProduct(aTmp[2,], bTmp);
    return(abs(r)<EPSILON)
}

##
# Checks if a point is right of a line. If the point is on the
# line, it is not right of the line.
# @param a line segment interpreted as a line
# @param b the point
# @return <code>true</code> if the point is right of the line,
#         <code>false</code> otherwise
##
isPointRightOfLine<-function(LineSegment.a, Point.b) {
    # Move the image, so that a.first is on (0|0)
 	aTmp = rbind( c(0, 0),
    						  c(
            				LineSegment.a[2,1] - LineSegment.a[1,1], LineSegment.a[2,2] - LineSegment.a[1,2])
            			)
    bTmp = c(Point.b[1]  - LineSegment.a[1,1], Point.b[2] - LineSegment.a[1,2])
    return( crossProduct(aTmp[2,], bTmp) < 0 )
}

##
# Check if line segment first touches or crosses the line that is 
# defined by line segment second.
#
# @param first line segment interpreted as line
# @param second line segment
# @return <code>true</code> if line segment first touches or
#                           crosses line second,
#         <code>false</code> otherwise.
##
lineSegmentTouchesOrCrossesLine<-function(LineSegment.a, LineSegment.b) {
    return( isPointOnLine(LineSegment.a, LineSegment.b[1,])
            || isPointOnLine(LineSegment.a, LineSegment.b[2,])
            || xor(isPointRightOfLine(LineSegment.a, LineSegment.b[1,]), 
                isPointRightOfLine(LineSegment.a, LineSegment.b[2,]))
                )
}


# GET BOUNDING BOX
getBoundingBox <-function(LineSegment){
	min.x<-min(LineSegment[,1])
	min.y<-min(LineSegment[,2])
	max.x<-max(LineSegment[,1])
	max.y<-max(LineSegment[,2])	
	
	return(matrix(c(min.x, min.y, max.x-min.x, max.y-min.y), ncol=2))
}

##
# Check if bounding boxes do intersect. If one bounding box
# touches the other, they do intersect.
# @param a first bounding box
# @param b second bounding box
# @return <code>true</code> if they intersect,
#         <code>false</code> otherwise. Where first row is coordinates and second row are width and height
##
doBoundingBoxesIntersect<-function(Box.a, Box.b) {
   a.min.x = Box.a[1,1]
   b.min.x = Box.b[1,1]
   a.max.x = Box.a[1,1] + Box.a[1,2]
   b.max.x = Box.b[1,1] + Box.b[1,2]
   
   a.min.y = Box.a[2,1]
   b.min.y = Box.b[2,1]
   a.max.y = Box.a[2,1] + Box.a[2,2]
   b.max.y = Box.b[2,1] + Box.b[2,2]
       
   return( (a.max.x >= b.min.x)  # a is right of b
   && (a.min.x <= b.max.x)  # a is left of b
   && (a.max.y >= b.min.y)  # a is below b
   && (a.min.y <= b.max.y)  # a is above b 
   )  
}


##
# Check if line segments intersect
# @param a first line segment
# @param b second line segment
# @return <code>true</code> if lines do intersect,
#         <code>false</code> otherwise
##
doLinesIntersect<-function(LineSegment.a, LineSegment.b) {
   box1 = getBoundingBox(LineSegment.a);
   box2 = getBoundingBox(LineSegment.b);
   return( doBoundingBoxesIntersect(box1, box2)
                && lineSegmentTouchesOrCrossesLine(LineSegment.a, LineSegment.b)
                && lineSegmentTouchesOrCrossesLine(LineSegment.b, LineSegment.a)
          )
}


getIntersection<-function(a, b){
    # the intersection [(x1,y1), (x2, y2)]
    #   it might be a line or a single point. If it is a point,
    #   then x1 = x2 and y1 = y2.  */


   if ( isTRUE(all.equal(a[1,1], a[2,1]) )) {
        # Case (A)
        # As a is a perfect vertical line, it cannot be represented
        # nicely in a mathematical way. But we directly know that
        #
        x1 <- a[1,1]
        x2 <- x1
        if ( isTRUE(all.equal(b[1,1] , b[2,1]) ))  {
          # Case (AA): all x are the same!
            # Normalize
            if(a[1,2] > a[1,2]) {
                tmp<-a[1,] ; a[1,] <- a[2,]; a[2,]<-tmp; #a={"first": a[2,, "second": a[1,};
            }
            if(b[1,2] > b[2,2]) {
                tmp<-b[1,] ; b[1,] <- b[2,]; b[2,]<-tmp; #b = {"first": b[2,, "second": b[1,};
            }
            if(a[1,2] > b[1,2]) {
                tmp<-a
                a<-b
                b<-tmp
            }

            # Now we know that the y-value of a[1, is the 
            # lowest of all 4 y values
           # this means, we are either in case (AAA):
            #   a: x--------------x
            #   b:    x---------------x
            # or in case (AAB)
            #   a: x--------------x
            #   b:    x-------x
            # in both cases:
            # get the relavant y intervall
            y1 = b[1,2];
            y2 = min(c(a[2,2], b[2,2]) );
        } else {
           # Case (AB)
            # we can mathematically represent line b as
            #     y = m*x + t <=> t = y - m*x
            # m = (y1-y2)/(x1-x2)
            m = (b[1,2] - b[2,2])/
                (b[1,1] - b[2,1]);
            t = b[1,2] - m*b[1,1];
            y1 = m*x1 + t;
            y2 = y1
        }
    } else{ if( isTRUE(all.equal(b[1,1],b[2,1])  )) {
       # Case (B)
        # essentially the same as Case (AB), but with
        # a and b switched
        x1 = b[1,1];
        x2 = x1;

        tmp = a;
        a = b;
        b = tmp;

        m = (b[1,2] - b[2,2])/
            (b[1,1] - b[2,1]);
        t = b[1,2] - m*b[1,1];
        y1 = m*x1 + t;
        y2 = y1
    } else {
       # Case (C)
        # Both lines can be represented mathematically
        ma = (a[1,2] - a[2,2])/
             (a[1,1] - a[2,1]);
        mb = (b[1,2] - b[2,2])/
             (b[1,1] - b[2,1]);
        ta = a[1,2] - ma*a[1,1];
        tb = b[1,2] - mb*b[1,1];
        if ( isTRUE(all.equal(ma , mb) )) {
            # Case (CA)
            # both lines are in parallel. As we know that they 
            # intersect, the intersection could be a line
            # when we rotated this, it would be the same situation 
            # as in case (AA)

            # Normalize
            if(a[1,1] > a[2,1]) {
                 tmp<-a[1,] ; a[1,] <- a[2,]; a[2,]<-tmp;# a = {"first": a["second"], "second": a["first"]};
            }
            if(b[1,1] > b[2,1]) {
                 tmp<-b[1,] ; b[1,] <- b[2,]; b[2,]<-tmp; #b = {"first": b["second"], "second": b["first"]};
            }
            if(a[1,1] > b[1,1]) {
                tmp = a;
                a = b;
                b = tmp;
            }

            # get the relavant x intervall
            x1 = b[1,1];
            x2 = min(a[2,1], b[2,1]);
            y1 = ma*x1+ta;
            y2 = ma*x2+ta;
        } else {
           # Case (CB): only a point as intersection:
            # y = ma*x+ta
            # y = mb*x+tb
            # ma*x + ta = mb*x + tb
            # (ma-mb)*x = tb - ta
            # x = (tb - ta)/(ma-mb)
            x1 = (tb-ta)/(ma-mb);
            y1 = ma*x1+ta;
            x2 = x1;
            y2 = y1;
        }
    }
    }

    if( (x1==x2)&&(y1==y2) ){
    			return(c(x1,y1))
    		}else{
    			return(c(x1,y1,x2,y2))
    		}
}


##########
# getPrincipalAxes
############


getPrincipalAxes<-function(contour, plot=F){
	load<-princomp(contour)$loadings
	
	slope <- load[2, ]/load[1, ]
	mn <- apply(contour, 2, mean)
	intcpt <- mn[2] - (slope * mn[1])
	height<-c(min(contour[,2]), max(contour[,2]))
	width<-c(min(contour[,1]), max(contour[,1]))
	y1<-height[1]-diff(height)*0.04
	y2<-height[2]+diff(height)*0.04
	x1<-width[1]-diff(width)*0.04
	x2<-width[2]+diff(width)*0.04
	#first components (usually vertical) get end points out of the contour
	PC_1_x<-( (c(y1, y2)-intcpt[1])/slope[1] )
	#second component (usually horizontal) get end points out of the contour
	PC_2_y <-(intcpt[2]+slope[2]*c(x1, x2) )
	
	
	intersect<-getIntersection(matrix(c(x1, x2, PC_2_y[1], PC_2_y[2]), ncol=2), matrix(c(PC_1_x[1], PC_1_x[2], y1, y2), ncol=2))
	
	if(plot){
		plot(contour, axes=F, type='l', ylab='', xlab='', asp=1)
		polygon(contour, col=gray(0.9))
		
		#first principal components
		arrows(PC_1_x[1], y1, PC_1_x[2], y2, length=0.15, code=3, lwd=2, col='darkred')	
		#second principal components
		arrows(x1, PC_2_y[1], x2, PC_2_y[2], length=0.15, code=3, lwd=2, col='darkred')
		
		points(intersect[1], intersect[2], pch=21, bg='red', cex=1.5)
	}
	
	principalaxes<-list( PC2=round(matrix(c(x1, x2, PC_2_y[1], PC_2_y[2]), ncol=2), 3),  PC1= round(matrix(c(PC_1_x[1], PC_1_x[2], y1, y2), ncol=2), 3), intersect=round(intersect,3) )
	
	return(principalaxes)
}

isLinesOverlapping<-function(line, contour){
    if(isTRUE(all.equal(line[1,2],line[2,2]))&isTRUE(all.equal(contour[1,2],contour[2,2]))){ #horizontal
        return(TRUE)
    }
    if(isTRUE(all.equal(line[1,1],line[2,1]))&isTRUE(all.equal(contour[1,1],contour[2,1]))){ #vertical
        return(TRUE)  
    }
}


getIntersectwithContour<-function(contour, line, plot=F){
	index<-c(1, rep(seq(2, nrow(contour)-1), each=2),  nrow(contour) )
	index<-matrix(index, ncol=2, byrow=T)
	intersectloop<-function(x){ doLinesIntersect(line, contour[index[x,],]) }
    overlaploop<-function(x){ isLinesOverlapping(line, contour[index[x,],]) }
	criteria<-lapply(1:nrow(index),  intersectloop )
	
	indexExtrapolate<-index[which(unlist(criteria)==TRUE),]


    if(length(indexExtrapolate)==0){
       criteria<-lapply(1:nrow(index),  overlaploop )
        indexExtrapolate<-index[which(unlist(criteria)==TRUE),]
    }

	coordinates<-c('x', 'y') #wierdness going on here
    if(length(indexExtrapolate)<3 ){
        coordinates <-rbind(coordinates, getIntersection(line, contour[indexExtrapolate,] ) )
        }else{
	for(i in 1:nrow(indexExtrapolate)){
		coordinates <-rbind(coordinates, getIntersection(line, contour[indexExtrapolate[i,],] ) )
	}
    }
	coordinates <- matrix(as.numeric(coordinates[-1,]), ncol=2, byrow=F)
	coordinates <-coordinates[!duplicated(coordinates), ]
	return(coordinates)
}	

reduceCoordinates<-function(coordinates){
	maxY<-which.max(coordinates[,2])
	minY<-which.min(coordinates[,2])
	maxX<-which.max(coordinates[,1])
	minX<-which.min(coordinates[,1])
	
	coordinates <-coordinates[c(minY,minX,maxY,maxX),]
	return(coordinates)	
}

minimumDistance<-function(midpoint, p.new){
	distance<-numeric()

	for(i in 1:nrow(p.new)){
		distance<-append(distance, sqrt((midpoint[1]-p.new[i,1])^2+(midpoint[2]-p.new[i,2])^2) )
	}#returns the popint with minimum distance
    return(p.new[which.min(distance),])	
}

#############
#  Generate homology points
############
automatic.correspondences<-function(contour, R, plot=F, cex=2){
	if(plot){
		par(mfrow=c(1,R))
	}
	
	PC<-getPrincipalAxes(contour, plot=FALSE)
	PC1<-PC$PC1
	PC2<-PC$PC2
	q<-PC$intersect
	
	#get bounaries for contour with 0.04 marginal
	height<-c(min(contour[,2]), max(contour[,2]))
	width<-c(min(contour[,1]), max(contour[,1]))
	y1<-height[1]-diff(height)*0.04
	y2<-height[2]+diff(height)*0.04
	x1<-width[1]-diff(width)*0.04
	x2<-width[2]+diff(width)*0.04
	
	p<-matrix(rep(0,8), ncol=2)
	p<-getIntersectwithContour(contour, PC$PC1)
	p<-rbind(p, getIntersectwithContour(contour, PC$PC2) )
    INDEX<-c(which.min(p[,2]), which.max(p[,1]), which.max(p[,2]), which.min(p[,1]) )
    #INDEX<-c(1,3,2,4)
    #p<-p[chull(p), ]
	p<-p[INDEX, ]
	
	r<-1
	n<-nrow(p)
	while(r<=R){
	p.out<-p
	p<-rbind(p, p[1,])
	p.tmp<-matrix(c(0,0), ncol=2)
	for(i in 1:(nrow(p)-1) ){
		midpoint<-colMeans(p[i:(i+1),])
		slope<-apply(p[i:(i+1),], 2, diff)
        if(slope[1]==0){ #vertical line 
         line<-cbind(p[i:(i+1),1], c(y1, y2) )
        }else if(slope[2]==0){ #horizontal line 
         line<-cbind( c(x1, x2), p[i:(i+1),2] )
        }else{
		  slope<-slope[2]/slope[1]
		  slope<- (-1/slope) #perpendicular -1
          intercpt<- midpoint[2] - slope*midpoint[1]
          line<-cbind(c(x1,x2), c(slope*c(x1, x2)+ intercpt) )
        }


		p.new<-getIntersectwithContour(contour, line) #get the points on the line that intersect the contour

        if(!is.null(nrow(p.new))){
		  p.new <-minimumDistance(midpoint, p.new) #get the point with the minimum distance to qi
        }
		p.tmp <-rbind(p.tmp, p.new)
		
		q <-rbind(q, midpoint)
	}
	p.tmp<-p.tmp[-1,]
	n<-2*n
	p<-matrix(rep(0, 2*n), ncol=2)
	p[seq(1,nrow(p)-1, by=2),]<-p.out
	p[seq(2,nrow(p), by=2),]<-p.tmp

	
	if(plot){
		plot(contour, type='l', axes=F, main=paste('level ',r), ylab='', xlab='', asp=1)
		polygon(contour, col=gray(0.9))
		lines(PC1, lwd=2)
		lines(PC2, lwd=2, lty=2)
		box()
		points(q, pch=21, bg='white', cex=cex)
		points(p.out, pch=16, cex=cex)
	}
	
	r<-r+1
	
	}
	
	return(list(p=p.out, q=matrix(as.numeric(q), ncol=2)))
}

getWholeSection<-function(contour){
	xmin<-min(contour[,1])
	contour[,1]<-(contour[,1]-xmin)
	contour<-cbind(c(contour[,1], rev(-contour[,1])),
					c(contour[,2], rev(contour[,2]) ) )
	return(contour)
}