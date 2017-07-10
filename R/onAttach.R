.onAttach<-function(...){
    packageStartupMessage("WholeBrain (version 0.0.9) \"Thurstone\" \n by Daniel Fürth, 2017, by using this software you agree to the EULA")
    if(.Platform$OS.type=="windows") {
  		quartz<-function(width,height){windows(width, height)}
  		assign("quartz", quartz, envir = .GlobalEnv) 
	}
}