# *QCMI*
# An R package for easy modeling, filtering, and quantifying putative biotic interactions of microbial communities.

<a href="https://github.com/joshualiuxu/qcmi/blob/main/qcmi.test.r/"><img src="https://github.com/joshualiuxu/qcmi/blob/main/data/fig.jpg" width=150 align="right" ></a>


## Overview

qcmi calculates the strength of biotic associations and quantifies the contributions on microbial diversity. 

qcmi provides some convenient verbs to make it easy to process data and results:


## Pipeline

  + Step 1. Construct ecological networks for microbial communities.  📜 

  + Step 2. Assign the ecological assembly processes to each significantly pair ASVs. 📈

  + Step 3. Quantify the strength of biotic associations to each local site. 📊

  + Step 4. Calculate the effects of biotic associations on alpha and beta diversity of microbial communities. ❤️


<img src="https://github.com/joshualiuxu/qcmi/blob/main/data/Figure1.jpg" width="90%" />


## Function

  + trans_ps() converts the data to phyloseq format.

  + filter_ps() filters OTU table by occurrence and abundance.

  + cal_network() infers ecological networks.

  + test_link_env() classifies the ecological associations to environmental filtering

  + test_link_dl() classifies the ecological associations to dispersal limitation

  + assigned_process() identifies ecological associations as environmental filtering and dispersal limitation

  + qcmi() quantifies the local intensity of microbial biotic interaction.

  + cal_alphacon() calculates the contributions of microbial interaction on alpha diversity

  + cal_betacon() calculates the contributions of microbial interaction on beta diversity

For a detailed introduction, please see qcmi.test.r.



## Installation

to get the development version from GitHub:
```r
# If devtools package is not installed, first install it
install.packages("devtools")
devtools::install_github("joshualiuxu/qcmi")
```

load the package:
```r
library(“qcmi”)
```

If you find a bug, please file a minimal reproducible example in the issues



## Usage

Please see the document of qcmi.test.r



