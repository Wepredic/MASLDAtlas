chooseCRANmirror(ind=1)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#for Shiny app
reqPkg = c("shiny","Rtools", "shinyhelper", "data.table", "Matrix", "hdf5r", "shinyWidgets","xlsx","DT",
           "ggplot2", "gridExtra", "magrittr", "ggdendro", "shinyhelper", "shinymanager","data.table", "Matrix", "hdf5r", "ggplot2", "gridExtra",
           "glue", "readr", "RColorBrewer", "R.utils", "reticulate","bslib","dplyr","shinydisconnect","shinycssloaders","shinyjs","ggpubr","shinyBS","stringr","fenr")
for(lib in reqPkg){
    if(!lib %in% installed.packages()){
        if(lib %in% available.packages()[,1]){
            install.packages(lib,dependencies=TRUE)
        }else {(BiocManager::install(lib))
    }}
}