# -----------------------------------------------------------------------------#
# Reducing consumer pressure increases community dissimilarity leading to 
# coral- or algal- dominance on a coral reef
#  
# 0_download_data
# -----------------------------------------------------------------------------#

# this script downloads data from EDI data repository
# creates necessary folders to store data and results 

# Packages --------------------------------------------------------------------#

library(EDIutils)
library(tidyverse)

# create folders for storing data, figures, and model outputs ------------------#

dir.create("data")
dir.create("data/data_summaries")
dir.create("figures")
dir.create("model outputs")
dir.create("model outputs/adonis")
dir.create("model outputs/betadisper")
dir.create("model outputs/betadisper_tukey")
dir.create("model outputs/simper")

# Download data ----------------------------------------------------------------#

# EDI data packageID is: "knb-lter-mcr.5056.1"
read_data_package(packageId="knb-lter-mcr.5056.1")

# read data package, download package zip to folder
read_data_package_archive(packageId="knb-lter-mcr.5056.1", path = "data/")

# unzip the file
unzip("data/knb-lter-mcr.5056.1.zip", exdir ="data/")




