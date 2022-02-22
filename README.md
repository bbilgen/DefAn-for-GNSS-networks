# DefAn-for-GNSS-networks
DefAn-for-GNSS-networks is a new scientific software encoded in R programming language for carrying out geodetic deformation analysis in GNSS networks using Bernese v5.2-based output files. The software determines absolute deformations in the GNSS network established for detecting deformations on engineering structures and earth's crust, illustrates input data and output results, and provides results as reports in two different formats HTML and PDF. The developed software uses two different types of ASCII input files with .SNX and .OUT extensions produced by Bernese v5.2 scientific GNSS software. In the first file, there are degrees of freedom for the adjustments of the observations, the root mean square error (RMSE) of the unit weighted observation, the observation date and the coordinates of all points in the geocentric Cartesian coordinate system. In the second file, there are the lower triangular elements of the covariance matrix regarding the adjustments of the observations. Since at least two measurement campaigns are required for geodetic deformation analysis, at least two of these files must be available with different observation dates. When the number of this file pair is more than two, the program determines the reference measurements according to the observation dates and performs all calculations based on the reference measurement campaign. 


# Instructions to download the software
Click on the green "Code" button and "Download.ZIP". Unzip the downloaded zip archive in "C:/" directory on your computer.

# Instructions to run
To run the software, first download the R v4.1.1 or newer version from https://cran.r-project.org and RStudio Desktop v1.4 or newer version from https://www.rstudio.com for your operating system and install them. It should be noted that RStudio requires R v3.0.1 or newer versions.

Open the R studio desktop and type the following code via RStudio Console for installing the required packages:
```
install.packages(‘dplyr’ , ‘ISwR’ , ‘matlib’  , ‘readr’, ‘rmarkdown’ , ‘rstudioapi’ , ‘stringr’ , dependencies = TRUE)
```

The following ‘tinytex’ package must also be installed so that the program generates PDF and HTML result reports. For this, type the following codes on the RStudio console, respectively:
```
install.packages(‘tinytex’ , dependencies = TRUE)
```
```
tinytex :: install_tinytex()
```

Afterwards, you must open the "DefAnforGNSS.Rproj" file, which exists in the downloaded and unzipped directory in "C:\", on your Rstudio desktop software by clicking on "File --> Open project". Type the following small block of code on RStudio console;
```
source("C:/DefAn-for-GNSS-networks-main/DefAnforGNSS.R")
```
Then, file selection dialog will start on your computer for selecting the .OUT extensioned input files. You must select both F1_162080.OUT and F1_162090.OUT files from the downloaded "Sample data" directory and click "Open" button. After that, file selection dialog will start again on your computer for selecting the .SNX extensioned input files.  You must select both F1_162080.SNX and F1_162090.SNX files from the downloaded "Sample data" directory and click "Open" button. Then, the software will start processing files, perform geodetic deformation analysis and write two output reports which are named as "Reporting.html" and "Reporting.pdf" in the downloaded and unzipped directory in "C:\". You will be able to access geodetic deformation analysis results by opening one of the reports produced by the software.

To use the software for your own GNSS deformation network, two types of Bernese v5.2 files (.SNX and .OUT extensioned files) of at least two GNSS campaign processing results must be copied in the "Sample data" directory. Then, the software should be run and these files should be selected to perform geodetic deformation analysis.

# Contents of the archived directory
In DefAn-for-GNSS-networks-main folder:
  * 'DefAnforGNSS.R' (download, and open this file on your RStudio Desktop to see the source code of the software) 
  * 'DefAnforGNSS.Rproj' (R project file to automatically open RStudio Desktop)
  * 'Reporting.Rmd' (R markdown file used automatically by the software's source codes to report analysis results)


In Sample data folder: 
  * contains sample input data files for users to test the application 


_Last Modified February 22, 2022, Burhaneddin Bilgen_

Please contact Burhaneddin Bilgen (bbilgen@ktun.edu.tr) with any problems, questions, or concerns.
