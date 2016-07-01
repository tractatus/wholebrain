
xprod <- function(a,b) {
    if(length(a)==3 && length(b)==3){
        getDet <- function(w,x,y,z) {
            return(det(matrix(c(w,x,y,z),ncol=2)))
        }
        return(c(getDet(a[2],b[2],a[3],b[3]),-getDet(a[1],b[1],a[3],b[3]),getDet(a[1],b[1],a[2],b[2]))) 
    }else{
        return(FALSE)
    }
}


signed.triangle.area<-function(A, B, C){
	u <- B - A;
	v <- C - A;
	cp <- xprod(c(u, 0), c(v, 0) )
	a <- sign(cp[1,3]) * 0.5 * sqrt(sum(cp^2))
	return(a)
}

contour_area<-function(c){
n<-nrow(c)
# get centroid of points
cog <- apply(c, mean, 2)

# get total area
total.contour.area <- signed.triangle.area(cog, c[n,], c[1,])
for(i in 1:(length(c)-1) ){
  total.contour.area <- total.contour.area + signed.triangle.area(cog, c[i,], c[i+1,])
}
return(total.contour.area)
}

area.normalize<-function(old.c){
	# get the centroid of the vertices of the contour
	cog <- mean(old.c)
	
	# get total area
	total.contour.area <- contour.area(old.c);

	# normalize by dividing the square root of the area
	ratio <- sqrt(abs(total.contour.area));

	for(i in seq_along(old.c)){
  		u <- old_c[i,] - cog
  		new.c[i,] <- cog + u/ratio - cog
	}
	return(new.c)
}

dist2<-function(x, c){
#DIST2	Calculates squared distance between two sets of points.
#
#	Description
#	D = DIST2(X, C) takes two matrices of vectors and calculates the
#	squared Euclidean distance between them.  Both matrices must be of
#	the same column dimension.  If X has M rows and N columns, and C has
#	L rows and N columns, then the result has M rows and L columns.  The
#	I, Jth entry is the  squared distance from the Ith row of X to the
#	Jth row of C.
#
#	See also
#	GMMACTIV, KMEANS, RBFFWD
#

#	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

ndata<-dim(x)[1]
ncentres<-dim(c)[1]
if(dim(x)[2] != dim(c)[2]){
	stop('Data dimension does not match dimension of centres')
}

n2 <- t(rep(1, ncentres) %*% t(apply(t(x^2), 2, sum))) + rep(1, ndata) %*% t(apply(t(c^2),2, sum)) - 2*(x%*%t(c))
return(n2)
}

linspace <- function(a, b, n = 100) {
    if (!is.numeric(a)) {
        stop(sprintf("argument %s must be numeric", sQuote("a")))
    } else if (!(length(a) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("a")))
    } 

    if (!is.numeric(b)) {
        stop(sprintf("argument %s must be numeric", sQuote("b")))
    } else if (!(length(b) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("b")))
    }

    if (!is.numeric(n)) {
        stop(sprintf("argument %s must be numeric", sQuote("n")))
    } else if (!(length(n) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("n")))
    }   

    n <- floor(n)    ## Undocumented but required
    if (n < 2) {
        b
    } else {
        seq(a, b, length = n)
    }
}

logspace <- function(a, b, n = 50) {
    if (b == pi) {
        b <- log10(pi)
    }

    10^linspace(a, b, n)
}

mod <- function(x, y) {
    ans <- x %% y
    ## Substitute x[off] in answer anywhere y[off] equals zero
    if (length(zero.off <- which(y == 0))) {
        ans[zero.off] <- if (length(x) == 1) x else x[zero.off]
    }
    return(ans)
}

rem <- function(x, y) {
    ans <- mod(x, y)
    if (!((x > 0 && y > 0) ||
          (x < 0 && y < 0))) {
        ans <- ans - y
    }

    return(ans)
}  		

sc_compute<-function(Bsamp,Tsamp,mean_dist,nbins_theta,nbins_r,r_inner,r_outer,out_vec){
# [BH,mean_dist]=sc_compute(Bsamp,Tsamp,mean_dist,nbins_theta,nbins_r,r_inner,r_outer,out_vec);
#
# compute (r,theta) histograms for points along boundary 
#
# Bsamp is 2 x nsamp (x and y coords.)
# Tsamp is 1 x nsamp (tangent theta)
# out_vec is 1 x nsamp (0 for inlier, 1 for outlier)
#
# mean_dist is the mean distance, used for length normalization
# if it is not supplied, then it is computed from the data
#
# outliers are not counted in the histograms, but they do get
# assigned a histogram
#

nsamp<-ncol(Bsamp)
in_vec<-out_vec==0

# compute r,theta arrays
r_array<-Re(sqrt(dist2(t(Bsamp), t(Bsamp))) ; # real is needed to
                                          # prevent bug in Unix version
theta_array_abs<-t(atan2(as.matrix(Bsamp[2,])%*%rep(1,nsamp)-as.matrix(rep(1,nsamp))%*%Bsamp[2,],as.matrix(Bsamp[1,:])%*%rep(1,nsamp)-as.matrix(rep(1,nsamp))%*%Bsamp[1,:]) );
theta_array<-theta_array_abs-as.matrix(Tsamp)%*%rep(1,nsamp)

# create joint (r,theta) histogram by binning r_array and
# theta_array

# normalize distance by mean, ignoring outliers
tmp<-r_array[in_vec,];
tmp<-tmp[,in_vec];
mean_dist<-mean(as.vector(tmp)); 

r_array_n<-kronecker(r_array,mean_dist)

# use a log. scale for binning the distances
r_bin_edges<-logspace(log10(r_inner),log10(r_outer),5)
r_array_q<-matrix(rep(0, nsamp*nsamp), ncol=nsamp)
for(m in 1:nbins_r){
   r_array_q<-r_array_q+(r_array_n<r_bin_edges[m])
}
fz<-r_array_q>0 # flag all points inside outer boundary

# put all angles in [0,2pi) range
theta_array_2 <- rem(rem(theta_array,2*pi)+2*pi,2*pi)
# quantize to a fixed set of angles (bin edges lie on 0,(2*pi)/k,...2*pi
theta_array_q <- 1+floor(theta_array_2/(2*pi/nbins_theta));

nbins<-nbins_theta*nbins_r;
BH<-matrix(rep(0, nsamp*nbins),nrow=nsamp,ncol=nbins);
for(n in 1:nsamp){
   fzn<-fz[n,]&in_vec
   Sn=sparse(theta_array_q(n,fzn),r_array_q(n,fzn),1,nbins_theta,nbins_r);
   BH[n,]<-t(Sn(:));
}

}

extract_shape_context<-function(Y){
# Get contour size
n <- nrow(Y);

# Set shape context parameters
nbins_theta<-12;
nbins_r<-5;
r_inner<-1/8;
r_outer<-2;
out_vec <- rep(0, n)

# Extract shape context
[D, mean_dist] = sc_compute(Y, zeros(1, n), [], nbins_theta, nbins_r, r_inner, r_outer, out_vec);
return(D)
}

shape.matching<-function(Y1, Y2, algorithm = 'aco', distance='chisquare'){
	Y1<-area.normalize(Y1)
	Y2<-area.normalize(Y2)

	# Get contour sizes
	n1 <- nrow(Y1)
	n2 <- nrow(Y2)

	# Switch contours if n1 > n2
	if(n1 > n2){
    	temp_n <- n1
    	n1 <- n2
    	n2 <- temp_n
    	temp_Y <- Y1
    	Y1 <- Y2
    	Y2 <- temp_Y
	}

	# Compute shape context for contours
	D1 <- extract_shape_context(Y1)
	D2 <- extract_shape_context(Y2)


}