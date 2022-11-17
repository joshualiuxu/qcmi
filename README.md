# QCMI*: Quantify Community-level Microbial Interactions*
# An R package for easy modeling, filtering, and quantifying putative biotic associations of microbial communities.

<a href="https://github.com/joshualiuxu/qcmi/blob/main/qcmi.test.r/"><img src="https://github.com/joshualiuxu/qcmi/blob/main/data/fig.jpg" width=150 align="right" ></a>


## Overview

qcmi calculates the strength of putative biotic associations and quantifies the contributions on microbial diversity. 

qcmi provides some convenient verbs to make it easy to process data and results:


## Pipeline

  + Step 1. Construct ecological networks for microbial communities.  📜 

  + Step 2. Assign the ecological assembly processes to each significantly pair ASVs. 📈

  + Step 3. Quantify the strength of putative biotic associations to each local site. 📊

  + Step 4. Calculate the effects of putative biotic associations on alpha and beta diversity of microbial communities. ❤️


<img src="https://github.com/joshualiuxu/qcmi/blob/main/data/Figure1.jpg" width="80%" />


## Function

  + trans_ps() converts the data to phyloseq format.

  + filter_ps() filters OTU table by occurrence and abundance.

  + cal_network() infers ecological networks.

  + rmt() filters correlation coefficient.

  + test_link_env() classifies the ecological associations to environmental filtering

  + test_link_dl() classifies the ecological associations to dispersal limitation

  + assigned_process() identifies ecological associations as environmental filtering and dispersal limitation to dig out putative biotic associations from complex ecological networks

  + qcmi() quantifies the strength of microbial biotic associations at the community level.

  + cal_alphacon() calculates the contributions of microbial associations on alpha diversity

  + cal_betacon() calculates the contributions of microbial associations on beta diversity

For a detailed introduction, please see Tutorial.r.



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

Please see the document of qcmi.tutorial.r


## Contributing

I’m happy to receive bug reports, suggestions, questions, and (most of
all) contributions to fix problems and add features. I personally prefer
using the `GitHub` issues system over trying to reach out to me in other
ways (personal e-mail, Twitter, etc.). Pull Requests for contributions
are encouraged.

Here are some simple ways in which you can contribute (in the increasing
order of commitment):

-   Read and correct any inconsistencies in the
    [documentation](https://github.com/joshualiuxu/qcmi/blob/main/Tutorial.r)

-   Raise issues about bugs or wanted features

-   Review code

-   Add new functionality (in the form of new plotting functions or
    helpers for preparing subtitles)


## Citation

Not yet QAQ 👻👻👻

On the way 🌱🌱🌱

