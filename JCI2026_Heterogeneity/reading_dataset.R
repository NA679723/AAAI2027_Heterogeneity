install.packages("readr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_0.3-21.tar.gz",force = TRUE)
install.packages("Matrix", repos = "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.2-17.tar.gz")
BiocManager::install("np",force = TRUE)
library(devtools); install_github('JeffreyRacine/R-Package-np')

library(readr)
library(cubature)
library(np)

data<-read_csv("insurance.csv")


bw<-npcdist(formula=data$charges~data$bmi,data=data)

          