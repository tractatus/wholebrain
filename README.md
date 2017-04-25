WholeBrain
=========

`WholeBrain` is intended as a software suite, meaning a collection of computer programs sharing a common interface and ability to integrate and exchange data with each other.

The purpose of WholeBrain is to provide a user-friendly and efficient way for scientist with minimal knowledge of computers to create anatomical maps and integrate this information with behavioral and physiological data for sharing on the web.

WholeBrain is conceived and created by Daniel Fürth, a PhD student in Konstantinos Meletis lab, at Department of Neuroscience, Karolinska Institutet.

### Example processing a single section

```R
#load package
library(wholebrain)
```
Get the image(s) you want to process.
```R
#set folder with 16-bit raw single-channel TIFF images
folder<-"~/Users/Documents/myexperiment/"

#get images
images<-get.images(folder)
```
Segment out neurons and the brain outline from autofluorescence.
```R
#segment neurons and brain outline
seg<-segment(images[1])
# to access parameters programmatically tryout edit(seg$filter) on the output
```
Register to coordinate 0.38 mm anterior-posterior from bregma.
```R
#register to atlas coordinate 0.38 mm from bregma anterior-posterior
regi<-registration(images[1], coordinate= 0.38, filter=seg$filter)
```

Get cell counts as well as stereotactic coordinates of each cell and display the results.
```R
#get cell counts in regions as well as stereotactic coordinates while inspecting registration results
dataset<-inspect.registration(regi, seg, forward.warps = TRUE)
# try get.cell.ids() to just get the cell dataset object without plotting registration results.
```
Write web-based interactive map of your tissue section.
```R
#set pixel resolution in microns 
pixel.resolution<-0.64
#name of channel imaged
protein <- "EGFP"
#make a web map output of your result
makewebmap(images[1], 
		seg$filter, 
		registration = regi, 
		dataset = dataset, 
		scale = pixel.resolution, 
		fluorophore = protein
	)
```

### For developers

Sorry for the poor documentation. I’m working on it. DM me if interested in contributing or specific dev questions.

### Installation instructions

http://www.wholebrainsoftware.org/cms/install/

### How to get started

http://www.wholebrainsoftware.org/cms/tutorials/

#### Author(s)

Daniel Fürth - [@wholebrainsuite](https://twitter.com/wholebrainsuite) - <daniel.furth@ki.se>



