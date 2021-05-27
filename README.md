# EFA Tutorial
 Accompanying scripts and material to conduct a temporal EFA on ERP data as explained in Scharf, Widmann, Bonmassar & Wetzel (under review, url: TBA).

The repository is structured as follows:
- The main scripts "01*.R" to "05*.R" conduct the EFA-based ERP analysis step by step as explained in the article.
- "tools" contains some custom functions which are called in the main scripts
- "data" contains the participant average data after pre-processing for both groups and conditions as MATLAB files<sup>1</sup>. 
- The subdirectories in "results" contain the output files of the respective script.
- All scripts with prefix "99" contain the code for custom functions which are needed in the main scripts. 

<sup>1</sup> <small>Note: The epochs in the set-files are participant averages, **not** trials. We use the electrode positions from these files to make plotting the factor tropographies easier. See *03b_topoplot_allFactors.R* for further details.</small> 

To use the code, please follow these steps:
1. Make sure that current versions of R (https://cran.r-project.org) and RStudio (https://www.rstudio.com) are installed 
2. Download the whole repository (Code -> Download ZIP)
3. Unzip to a local folder on your computer
4. Open the file "00_EFAtutorial.Rproj". 
5. This should open a separate instance of Rstudio with all necessary paths.
6. You can open and run all code from this point.

In case you encounter bugs or unexpected behavior, please, open an issue or report via mail to 

florian.scharf dot uni-muenster dot de 
