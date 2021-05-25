# EFA Tutorial
 Accompanying scripts and material to conduct a temporal EFA on ERP data as explained in Scharf, Widmann, Bonmassar & Wetzel (under review, url: TBA).

The repository is structured as follows:
- The main scripts "01*.R" to "05*.R" conduct the EFA-based ERP analysis step by step as explained in the article.
- "data" contains the participant average data after pre-processing for both groups and conditions as MATLAB files.
- The subdirectories in "results" contain the output files of the respective script.
- All scripts with prefix "99" contain the code for custom functions which are needed in the main scripts. 

To use the code, please follow these steps:
1. Make sure that current versions of R (https://cran.r-project.org) and RStudio (https://www.rstudio.com) are installed 
2. Download the whole repository (Code -> Download ZIP)
3. Unzip to a local folder on your computer
4. Open the file EFAtutorial.Rproj. 
5. This should open a separate instance of Rstudio with all necessary paths.
6. You can open and run all code from this point.

In case you encounter bugs or unexpected behavior, please, open an issue or report via mail to 

florian.scharf dot uni-muenster dot de 
