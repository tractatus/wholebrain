
file_ext<-function(x){
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
}

#' Get images from directory
#'
#' Find dimensions of stitched image from number of tiles.
#' @param input a list of file names.
#' @keywords grid
#' @export
#' @examples
#' get.images('/Volumes/Seagate Slim Drive/D159/FFC_section001_1_ownFFC')
get.images<-function(x, type='tif'){
  image.files<-which(file_ext(dir(x))== type)
  
  image.files <-paste(path.expand(x), dir(x)[image.files], sep='/' )
  return(image.files)
}
