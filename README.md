WholeBrain
=========

`WholeBrain` is intended as a software suite, meaning a collection of computer programs sharing a common interface and ability to integrate and exchange data with each other.

The purpose of WholeBrain is to provide a user-friendly and efficient way for scientist with minimal knowledge of computers to create anatomical maps and integrate this information with behavioral and physiological data for sharing on the web.

WholeBrain is conceived and created by Daniel Fürth, a PhD student in Konstantinos Meletis lab, at Department of Neuroscience, Karolinska Institutet.

### Example processing a single section

```
library(wholebrain)
folder<-’~/Users/Documents/myexperiment/’
images<-get.images(folder)

seg<-segment(images[1])
regi<-registration(images[1], coordinate= 0.38, filter=seg$filter)

dataset<-inspect.registration(regi, seg, forward.warps = TRUE)

pixel.resolution<-0.64
protein <- "EGFP"
makewebmap(images[1], seg$filter, registration = regi, dataset = dataset, scale = pixel.resolution, fluorophore = protein)
```

### For developers

Sorry for the poor documentation. I’m working on it. DM me if interested in contributing or specific dev questions.

### Installation instructions

http://www.wholebrainsoftware.org/cms/install/

### How to get started

http://www.wholebrainsoftware.org/cms/tutorials/

#### Author(s)

Daniel Fürth - [@wholebrainsuite](https://twitter.com/wholebrainsuite) - <daniel.furth@ki.se>



