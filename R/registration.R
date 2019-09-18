registration <- function (input,
                          coordinate = NULL,
                          plane = "coronal",
                          right.hemisphere = NULL, 
                          interpolation = "tps",
                          intrp.param = NULL,
                          brain.threshold = 200, 
                          blurring = c(4, 15),
                          pixel.resolution = 0.64,
                          resize = (1/8)/4, 
                          correspondance = NULL,
                          resolutionLevel = c(4, 2),
                          num.nested.objects = 0, 
                          display = TRUE,
                          plateimage = FALSE,
                          forward.warp = FALSE, 
                          filter = NULL,
                          output.folder = "../",
                          batch.mode = FALSE, 
                          channel = 0,
                          verbose = TRUE){

  # We need to access this data that is part of the package
  data("EPSatlas", package = "wholebrain")
  # Input check #####
  # check platform and change batch.mode
  if (.Platform$OS.type == "windows" | grepl("linux-gnu", R.version$os)) {
    batch.mode = TRUE
  }
  # If not given, make coordinate = 0 
  if (is.null(coordinate)) {
    if (is.null(correspondance)) {
      coordinate <- 0
    }
    else {
      # Look for the coordinate on correspondance object
      coordinate <- correspondance$coordinate
    }
  }
  
  # File checking #####
  file <- as.character(input)
  if (!file.exists(file)) 
    stop(file, ", file not found")
  # Generate filenames #####
  file <- path.expand(file)
  # Get the name of file
  outputfile <- basename(file)
  # remove extension
  outputfile <- tools::file_path_sans_ext(outputfile)
  
  # Set the output folder path
  defaultwd <- getwd()
  if (file.exists(output.folder)) {
    parentpath <- output.folder
  }
  else if (output.folder == "./") {
    parentpath <- dirname(input)[1]
  }
  else if (output.folder == "../") {
    parentpath <- dirname(dirname(input))[1]
  }
  # Generate final output folder name
  outfolder <- paste("output", outputfile, sep = "_")
  # directory change needed to feed create.output.directory function
  setwd(parentpath)
  # create new directory for output files
  # verbose wrapper of dir.create (wholebrain::create.outupt.directory)
  create.output.directory(outfolder, verbose = verbose)
  # Change directory back
  setwd(defaultwd)
  outputfile <- paste(parentpath, outfolder,
                      paste("Registration", 
                            outputfile, sep = "_"), sep = "/")
  
  # Get Contours #####
  # atlasIndex object contains metadata for sections
  # subset according to plane of view
  # try View(atlasIndex)
  plate.width <- 1
  SAGITTAL <- TRUE
  if (plane == "sagittal") {
    EPSatlas <- SAGITTALatlas
    # subset sagittal plate from atlasIndex
    atlasIndex <- atlasIndex[atlasIndex$plane == "sagittal", 
                             ]
    plate.width <- 1.159292
    SAGITTAL <- !SAGITTAL
  }
  else {
    # subset coronal plate from atlasIndex
    atlasIndex <- atlasIndex[atlasIndex$plane == "coronal", 
                             ]
  }
  # Try to look for plate image (only if provided in file)
  if (plateimage != FALSE) {
    plateimage <- as.character(plateimage)
    if (!file.exists(plateimage)) 
      stop(plateimage, ", file not found")
    plateimage <- path.expand(plateimage)
  }
  # If filter was not provided set intensities for display 
  if (is.null(filter)) {
    MaxDisp <- 0
    MinDisp <- 0
  }
  # Otherwise get the info from the provided filter
  else {
    MaxDisp <- filter$Max
    MinDisp <- filter$Min
    brain.threshold <- filter$brain.threshold
    resize <- filter$resize
    # blurring[1] is used for image
    # blurring[2] is used for atlas
    blurring[1] <- filter$blur
  }
  if (is.null(MinDisp)) {
    MinDisp <- 0
  }
  # If no correspondance provided
  if (is.null(correspondance)) {
    # Try to look for pre-processed contour
    if(!is.null(filter$biggest)){
      contourInput <- filter$biggest
      # This call to unique might be problematic
      # The original code looks for contour.ID of the which.min(x)
      contoursI <- unique(filter$biggest$contour.ID)
      
      # We need to resize, otherwise the transformationgrid$mx[index] command
      # will give out of bounds errors
      # let's try something like
      
      ##### THIS IS KEY #########
      
      contourInput$x <- resize * contourInput$x
      contourInput$y <- resize * contourInput$y
      
      
    } else {
      # Get contour from C .Call to "getcontour"
      # Returns list with xy coordinates and contour ID
      contourInput <- get.contour(file, 
                                  thresh = brain.threshold,
                                  resize = resize,
                                  blur = blurring[1],
                                  num.nested.objects = num.nested.objects, 
                                  channel = channel, display = FALSE)
      
      # get the contour.ID number of the min x coordinate
      # This is done for each id
      
      contoursI <- as.numeric(names(sort(tapply(contourInput$x, 
                                                contourInput$contour.ID, min))))
      
      # Slower but consider readability
      # contoursI <- contourInput %>% 
      #	as.data.frame() %>%
      #	group_by(contour.ID) %>%
      #	summarise(x_min = min(x)) %>% 
      #	arrange(contour.ID) %>%
      #	pull(contour.ID)
    }
    
    
    # Get atlas binary image with helper function
    # see ?wholebrain::get.atlas.image
    filename <- get.atlas.image(coordinate,
                                right.hemisphere = right.hemisphere, 
                                plane = plane)
    # Get contour for the atlas image
    contourAtlas <- get.contour(filename,
                                resize = 1,
                                # blurring[2] is used for atlas
                                blur = blurring[2], 
                                num.nested.objects = num.nested.objects,
                                display = FALSE)
    
    contours <- as.numeric(names(sort(tapply(contourAtlas$x, 
                                             contourAtlas$contour.ID, min))))
    
    # Generate empty lists 
    cor.pointsInput <- list()
    cor.pointsAtlas <- list()
    
    # Enter the void :)
    for (i in 1:length(contours)) {
      # Toggle resLevel value depending on whether contours is 0
      if (contours[i] == 0) {
        # use low level
        resLevel <- resolutionLevel[1]
      }
      else {
        # use high level
        resLevel <- resolutionLevel[2]
      }
      # Perform correspondences for atlas image
      # This makes use of wholebrain::automatic.correspondences and other helper functions
      # Helper functions need annotation!
      # We could potentially perform
      # contourAtlas_frame <- contourAtlas %>% as.data.frame() # put outside of the loop
      # and feed this instead of cbind (readability improvement ?)
      # as.matrix(contourAtlas_frame)[which(contourAtlas_frame$contour.ID == contours[i]), ]
      cor.pointsAtlas[[i]] <- automatic.correspondences(cbind(contourAtlas$x[which(contourAtlas$contour.ID == 
                                                                                     contours[i])], contourAtlas$y[which(contourAtlas$contour.ID == 
                                                                                                                           contours[i])]), resLevel, plot = FALSE)
      # Perform correspondences for Input image
      cor.pointsInput[[i]] <- automatic.correspondences(cbind(contourInput$x[which(contourInput$contour.ID == 
                                                                                     contoursI[i])], contourInput$y[which(contourInput$contour.ID == 
                                                                                                                            contoursI[i])]), resLevel, plot = FALSE)
    }
    # Get the correspondance coordinates into a list 
    cor.points <- list(atlas = cor.pointsAtlas, input = cor.pointsInput)
    
    
    centroidAtlas <- cor.points$atlas[[1]]$q[1, ]
    # Where does this 456/2 come from ? how flexible is it ?       
    offsetAtlas <- centroidAtlas - 456/2
    targetP.x <- numeric()
    referenceP.x <- numeric()
    targetP.y <- numeric()
    referenceP.y <- numeric()
    shape <- numeric()
    centroidNorm <- centroidAtlas - cor.points$input[[1]]$q[1, 
                                                            ]
    for (i in 1:length(contours)) {
      referenceP.x <- append(referenceP.x, 4 * (cor.points$atlas[[i]]$p[, 
                                                                        1] - centroidNorm[1]))
      referenceP.y <- append(referenceP.y, 4 * (cor.points$atlas[[i]]$p[, 
                                                                        2] - centroidNorm[2]))
      targetP.x <- append(targetP.x, 4 * cor.points$input[[i]]$p[, 
                                                                 1])
      targetP.y <- append(targetP.y, 4 * cor.points$input[[i]]$p[, 
                                                                 2])
      shape <- append(shape, rep(i, length(cor.points$input[[i]]$p[, 
                                                                   2])))
    }
  }
  else {
    correspondance$correspondance <- na.omit(correspondance$correspondance)
    targetP.x <- correspondance$correspondance[, 1]
    targetP.y <- correspondance$correspondance[, 2]
    referenceP.x <- correspondance$correspondance[, 3]
    referenceP.y <- correspondance$correspondance[, 4]
    shape <- correspondance$correspondance[, 5]
    centroidNorm <- correspondance$centroidNorm
    centroidAtlas <- correspondance$centroidAtlas
  }
  resizeP <- resize * 4
  if (interpolation == "cpd") {
    if (is.null(intrp.param)) {
      intrp.param <- list(beta = 3, lambda = 3, gamma = 0.7, 
                          sigma = 0, max.iter = 150)
    }
    transformationgrid <- cpdNonrigid(file, targetP.x, targetP.y, 
                                      referenceP.x, referenceP.y, resizeP, MaxDisp, MinDisp, 
                                      outputfile, intrp.param$beta, intrp.param$lambda, 
                                      intrp.param$gamma, intrp.param$sigma, intrp.param$max.iter)
  }
  else {
    # Call C "ThinPlateRegistration"
    # Returns list: mx, my, width, height
    # List of 4
    # $ mx    : num [1:1928, 1:2400] 1 1 0 0 0 0 -1 -1 -1 -1 ...
    # $ my    : num [1:1928, 1:2400] 0 1 1 2 3 4 5 6 7 8 ...
    # $ width : int 1200
    # $ height: int 964
    
    transformationgrid <- .Call("ThinPlateRegistration", 
                                file, targetP.x, targetP.y, referenceP.x, referenceP.y, 
                                resizeP, MaxDisp, MinDisp, outputfile, channel)
  }
  
  # find the closest plate to the coordinate
  min_dist <- min(abs(coordinate - atlasIndex$mm.from.bregma))
  k <- which(abs(coordinate - atlasIndex$mm.from.bregma) == min_dist)
  
  
  # EPSatlas #####
  # Where does this scaling factor come from, how flexible is it ? 
  xmin <- (plane != "sagittal") * (min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 
                                     97440/2)
  numPaths <- EPSatlas$plates[[k]]@summary@numPaths
  # Where does this scaling factor come from, how flexible is it ? 
  scale.factor <- 456/97440
  style <- EPSatlas$plate.info[[k]]$style
  outlines <- list()
  for (i in 1:numPaths) {
    if (is.null(right.hemisphere)) {
      # Do the scaling math 
      xr <- 4 * ((EPSatlas$plates[[k]][[i]]@paths$path@x - 
                    xmin) * scale.factor - centroidNorm[1])
      yr <- 4 * (((-EPSatlas$plates[[k]][[i]]@paths$path@y) * 
                    scale.factor + 320) - centroidNorm[2])
      xl <- 4 * ((-(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                      xmin - (plate.width * 97440)/2) + (plate.width * 
                                                           97440)/2) * scale.factor - centroidNorm[1])
      yl <- 4 * (((-EPSatlas$plates[[k]][[i]]@paths$path@y) * 
                    scale.factor + 320) - centroidNorm[2])
      index <- cbind(as.integer(round(yr)), as.integer(round(xr)))
      # TODO: These lines throw uninformative subset errors
      # Code should check whether is even possible to do subset
      # The code relies on the resize parameter...should that be improved ???
      # Unless the resize is within bound, it will break here with no info

      xrT <- (xr + (transformationgrid$mx[index] - xr))
      yrT <- (yr + (transformationgrid$my[index] - yr))
      index <- cbind(as.integer(round(yl)), as.integer(round(xl)))
      xlT <- (xl + (transformationgrid$mx[index] - xl))
      ylT <- (yl + (transformationgrid$my[index] - yl))
    }
    else {
      if (right.hemisphere == SAGITTAL) {
        xr <- 4 * ((EPSatlas$plates[[k]][[i]]@paths$path@x - 
                      xmin) * scale.factor - centroidNorm[1])
        yr <- 4 * (((-EPSatlas$plates[[k]][[i]]@paths$path@y) * 
                      scale.factor + 320) - centroidNorm[2])
        xl <- NA
        yl <- NA
        index <- cbind(as.integer(round(yr)), as.integer(round(xr)))
        xrT <- (xr + (transformationgrid$mx[index] - 
                        xr))
        yrT <- (yr + (transformationgrid$my[index] - 
                        yr))
        xlT <- NA
        ylT <- NA
      }
      else {
        xr <- NA
        yr <- NA
        xl <- 4 * ((-(EPSatlas$plates[[k]][[i]]@paths$path@x - 
                        xmin - (plate.width * 97440)/2) + (plate.width * 
                                                             97440)/2) * scale.factor - centroidNorm[1])
        yl <- 4 * (((-EPSatlas$plates[[k]][[i]]@paths$path@y) * 
                      scale.factor + 320) - centroidNorm[2])
        xrT <- NA
        yrT <- NA
        index <- cbind(as.integer(round(yl)), as.integer(round(xl)))
        xlT <- (xl + (transformationgrid$mx[index] - 
                        xl))
        ylT <- (yl + (transformationgrid$my[index] - 
                        yl))
      }
    }
    outlines[[i]] <- list(xr = xr, yr = yr, xl = xl, yl = yl, 
                          xrT = xrT, yrT = yrT, xlT = xlT, ylT = ylT)
  }
  # Display results #####
  if (display) {
    par(yaxs = "i", xaxs = "i", bg = "black", mar = c(0, 
                                                      0, 0, 0))
    img <- paste(outputfile, "_undistorted.png", sep = "")
    img <- readPNG(img)
    if (length(dim(img)) > 2) {
      img = as.raster(img[, , ])
    }
    else {
      img = as.raster(img[, ])
    }
    if (batch.mode) {
      img <- apply(img, 2, rev)
    }
    par(xaxs = "r", yaxs = "r")
    plot(c(0, dim(img)[2] * 2), c(0, dim(img)[1]), axes = F, 
         asp = 1, col = 0, xlab = "", ylab = "", ylim = c(dim(img)[1], 
                                                          0))
    polygon(c(-5, dim(img)[2] + 5, dim(img)[2] + 5, -5), 
            c(-5, -5, dim(img)[1] + 5, dim(img)[1] + 5), col = "orange")
    polygon(c(-5 + dim(img)[2], 2 * dim(img)[2] + 5, 2 * 
                dim(img)[2] + 5, -5 + dim(img)[2]), c(-5, -5, dim(img)[1] + 
                                                        5, dim(img)[1] + 5), col = "purple")
    rasterImage(img, 0, 0, dim(img)[2], dim(img)[1])
    rasterImage(img, dim(img)[2], 0, 2 * dim(img)[2], dim(img)[1])
    abline(v = dim(img)[2], lwd = 2, col = "white")
    if (is.null(right.hemisphere)) {
      lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xr, outlines[[x]]$yr, border = "orange")
      })
      lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xl, outlines[[x]]$yl, border = "orange")
      })
      lapply(1:numPaths, function(x) {
        polygon(dim(img)[2] + outlines[[x]]$xrT, outlines[[x]]$yrT, 
                border = "purple")
      })
      lapply(1:numPaths, function(x) {
        polygon(dim(img)[2] + outlines[[x]]$xlT, outlines[[x]]$ylT, 
                border = "purple")
      })
    }
    else {
      if (right.hemisphere == SAGITTAL) {
        lapply(1:numPaths, function(x) {
          polygon(outlines[[x]]$xr, outlines[[x]]$yr, 
                  border = "orange")
        })
        lapply(1:numPaths, function(x) {
          polygon(dim(img)[2] + outlines[[x]]$xrT, outlines[[x]]$yrT, 
                  border = "purple")
        })
      }
      else {
        lapply(1:numPaths, function(x) {
          polygon(outlines[[x]]$xl, outlines[[x]]$yl, 
                  border = "orange")
        })
        lapply(1:numPaths, function(x) {
          polygon(dim(img)[2] + outlines[[x]]$xlT, outlines[[x]]$ylT, 
                  border = "purple")
        })
      }
    }
    lapply(1:length(targetP.x), function(x) {
      points(c(targetP.x[x] + dim(img)[2], referenceP.x[x]), 
             c(targetP.y[x], referenceP.y[x]), pch = c(19), 
             col = "black", cex = 1.8)
      text(c(targetP.x[x] + dim(img)[2], referenceP.x[x]), 
           c(targetP.y[x], referenceP.y[x]), label = x, 
           col = "white", cex = 0.7)
    })
  }
  style <- EPSatlas$plate.info[[k]]$style
  returnlist <- list(atlas = list(outlines = outlines, numRegions = numPaths, 
                                  col = style), transformationgrid = transformationgrid, 
                     correspondance = data.frame(targetP.x, targetP.y, referenceP.x, 
                                                 referenceP.y, shape), centroidAtlas = centroidAtlas, 
                     centroidNorm = centroidNorm, coordinate = coordinate, 
                     resize = resize, outputfile = outputfile, plane = plane)
  if (forward.warp) {
    returnlist <- get.forward.warp(returnlist)
  }
  return(returnlist)
}
