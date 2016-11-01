.onAttach<-function(...){
    packageStartupMessage("WholeBrain (version 0.0.7) \"Thurstone\" \n by Daniel Fürth, 2016, by using this software you agree to the EULA")
    if(.Platform$OS.type=="windows") {
  		quartz<-function() windows(width, height)
	}
}