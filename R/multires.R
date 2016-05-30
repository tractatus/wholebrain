mrd <- function(input, scales=6, family='db2', output='waveletoutput') {
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
    
    

    ## expand path
    file <- path.expand(file)
    #############################
    #           CALL            #
    #############################
    a <- .Call("multiresolutiondecomposition", file, scales, family, outputfile)
    cat("IMAGE PROCESSED\n")

  }