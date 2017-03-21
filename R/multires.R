mrd <- function(input, scales=6, family='db2', output='waveletoutput',  scale.of.interest = 3, energy.trace=0, coherency=0, orientation=0, max=0, min=0, maskoriginal=FALSE, out.SOI.only=FALSE) {
    file <- as.character(input)[1]
    family <- as.character(family)[1]
    outputfile <- as.character(output)[1]
    scales<-as.integer(scales)
    
    ## check for existence
    if(!file.exists(file))
    stop(file, "not found")
    
    create.output.directory('approximations_coefficents')
    lapply(1:scales-1, function(x){create.output.directory( paste('d',x, sep='') ) } )
    create.output.directory('output')
    if(energy.trace>0){
        create.output.directory('trace')
        if(energy.trace==TRUE){
            energy.trace<-10
        }
    }
    if(coherency){
        create.output.directory('coherency')
    }
    if(orientation){
        create.output.directory('orientation')
    }

    ## expand path
    file <- path.expand(file)
    #############################
    #           CALL            #
    #############################
    maskoriginal<-as.integer(maskoriginal)
    a <- .Call("multiresolutiondecomposition", file, scales, family, outputfile, scale.of.interest, energy.trace, coherency, orientation, ceiling(max), floor(min), maskoriginal, out.SOI.only)
    cat("IMAGE PROCESSED\n")

  }