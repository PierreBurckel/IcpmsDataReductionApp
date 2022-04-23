# ICP-MS data reduction app

## Description

This is an R Shiny application that facilitates data reduction for ICP-MS analyses (currently for Agilent format files).

## Installation

### Installing R and RStudio

Go to [this page](https://cran.rstudio.com/) to download the latest R version (chose the link adapted to your operating system).
Follow the instructions on the website/installer to properly install R.

Go to [this page](https://www.rstudio.com/products/rstudio/download/#download) to download the latest RStudio version. RStudio is the main R Integrated Development Environment (IDE) and  facilitates writing, reading and running R code.


### Preparing the R scripts

On [this page](https://github.com/PierreBu/R_ICPMS_Process/tree/master) (master branch, this is the stable version of the application) click on the green **Code** button and select **Download ZIP**. Extract the ZIP folder where you want to run the ICP-MS data reduction software. Keep all files together in the same folder.

Open **IcpmsDataReductionApp_main.R** with RStudio (open RStudio, go to **File**, **Open File** then navigate to **IcpmsDataReductionApp_main.R**). On top of the R Script, there should be a warning symbol with a message stating that some packages are not installed (see image below). You can click on **Install** to install the required packages (this will take several minutes, an internet connexion is required).


![Rstudio screenshot](https://github.com/PierreBurckel/IcpmsDataReductionApp/blob/master/RStudioCaptureEcran_mod.png)

If the message is not displayed, either the packages are already installed and you can continue, or you'll need to manually install the packages. To do that, uncomment all the `install.packages` lines (see image) by selecting them then go to the **Code** menu and click on **Comment/Uncomment Lines**. With the lines still selected, go back to **Code** then click on **Run Selected Line(s)**. This will install the packages. Once they are installed (this takes several minutes), you can comment the lines back (same procedure as to uncomment). Save the script by going to the **File** menu then clicking on **Save**. Saving is not necessary to keep the packages installed, but since the code was changed RStudio wants the script to be saved prior to execution. Alternatively, if you don't want to save unexpected changes, you can close the script, click on **Don't Save** to discard changes, then open it back up. The packages will still be installed and the script will be the original.

## Run

To run the application, open **IcpmsDataReductionApp_main.R** with RStudio and click on the **Run App** button (next to a green "Play" icon, see image above) on top of the R Script. This will run the application with your default web browser.

## Tutorials

Here are 4 tutorial videos for processing data on the Agilent 7900 ICP-MS. Video n°4 covers the use of the R-Shiny application. I highly encourage you to download the videos prior to watching as they are several hundreds of megabytes large.
- Part1, [Opening and checking data](https://drive.google.com/file/d/1mFeMndzGmAtN5Qt_tsqvyM2BrRVbyo06/view?usp=sharing)
- Part2, [Data correction](https://drive.google.com/file/d/1mOb3AqzBAstOO8Fs1d2B7wsbLvItuS8a/view?usp=sharing)
- Part3, [Data processing (Excel)](https://drive.google.com/file/d/1YT1coEkk1zpPVjIeC4uTmfB8GYZMdzC7/view?usp=sharing)
- Part4, [Data processing (R)](https://drive.google.com/file/d/1-ycANNSknwSEfpGpWF8atFMBAEwdLh2z/view?usp=sharing)

**Attention, some information from the videos are not up to date anymore, please refer to the Patch information section for more information**

## Patch information

- In the current release, the name of the analytes in the standard file should be unique and refer exactly to the names in the data file. For instance, if lithium 7 was measured in no gas mode, its name will be " 7  Li  [ No Gas ] " in the data file and it should be also " 7  Li  [ No Gas ] " in the standard file, not just "Li" anymore.

## Questions

If you have issues running the application, please send me an e-mail with a copy of the error message displayed in the console. For any questions, my e-mail address is burckel@ipgp.fr 
