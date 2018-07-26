.onAttach<-function(...){
    packageStartupMessage("WholeBrain (version 0.1.1) \"Thurstone\" \n by Daniel Fürth, 2018, by using this software you agree to the EULA")
    
 	get_os <- function(){
  		sysinf <- Sys.info()
  		if (!is.null(sysinf)){
    		os <- sysinf['sysname']
    	if (os == 'Darwin')
      		os <- "osx"
  		} else { ## mystery machine
    		os <- .Platform$OS.type
    	if (grepl("^darwin", R.version$os))
      		os <- "osx"
    	if (grepl("linux-gnu", R.version$os))
     		 os <- "linux"
  		}
  		tolower(os)
	}   

    if(get_os()=="windows") {
  		quartz<-function(width,height){windows(width, height)}
  		assign("quartz", quartz, envir = .GlobalEnv) 
	}

	if(get_os()=="linux") {
  		quartz<-function(width,height){X11(width, height)}
  		assign("quartz", quartz, envir = .GlobalEnv) 
	}

}