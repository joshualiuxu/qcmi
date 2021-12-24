#**Overview**

qcmi calculate the local intensity of microbial biotic interaction and quantify the contributions on microbial diversity patterns. 

qcmi provides some convenient verbs to make it easy to process data and results:

trans_ps() converts the data to phyloseq format.

filter_ps() filters OTU table by occurrence and abundance.

cal_network() calculates correlation network.

test_link_env() classifies ecological associations of environmental filtering

test_link_dl() classifies ecological associations of dispersal limitation

assigned_process() identifies ecological associations as environmental filtering and dispersal limitation

qcmi() quantifies the local intensity of microbial biotic interaction.

cal_alphacon() calculates the contributions of microbial interaction on alpha diversity

cal_betacon() calculates the contributions of microbial interaction on beta diversity

For a detailed introduction, please see qcmi.test.r.

****

#**Installation**

~to get the development version from GitHub:
install.packages("devtools")
devtools::install_github("joshualiuxu/qcmi")

~load the package:
Library(“qcmi”)

If you find a bug, please file a minimal reproducible example in the issues.

****

#**Usage**

Please see the document of qcmi.test.r

****

#**Pipeline**

![image]( https://github.com/joshualiuxu/qcmi/blob/main/data/FigS4.png)
