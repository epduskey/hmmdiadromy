# Weclome to hmmdiadromy

The HMM folder contains all code necessary to replicate the simulation study found in Duskey and Sullivan 2022, and referenced and altered in Duskey, Higgs, Fox, Breece, and Sullivan 2021.  We will update with the full citations when the articles have been accepted for publication.  Briefly, the code can simulate migratory paths for diadromous fish, and then run a Bayesian Hidden markov model (HMM) with nested Markov movement processes describing upriver and downriver migration.  The model also includes a logistic model of habitat preference, and a basical model of search and detection.

## Requirements

The Bayesian model contained in HMM run in JAGS via the R package jagsUI.  See the following website for instructions on how to install JAGS:

https://mcmc-jags.sourceforge.io/

## Usage

Follow these steps prior to first time usage to ensure the code runs properly on your machine:

1. Download the "hmmdiadromy-main" zipped folder from Github in its entirety
2. Unzip and place "hmmdiadromy-main" in your preferred directory
3. Navigate to the Code folder and open the file "depend.R" in R or RStudio
4. Confirm that all listed packages are installed on your machine; uncomment and run lines corresponding to packages you have not yet installed, then re-comment
5. Change the "mypath" variable to the file path which contains the hmmdiadromy folder
6. Save and close "depend.R"

Upon all subsequent uses, we recommend running the .R files in the following order:

1. depend.R -- to set mypath; alternatively, change setwd() within other R files
2. RUN_hmm.R -- to simulate data and run all models for HMM data
3. RUN_smm.R -- to simulate data and run all models for Standard Markov data
4. allplot.R -- to recreate all figures and supplemental figures in Duskey and Sullivan 2022

All files not appearing here are sourced in one or more of the files listed above, and typically contain functions necessary to run the code in each script.

NOTE: There are several nearly empty output folders contained in HMM, containing only a textfile called ".keep".  These are meant to contain and organize your own output.  Feel free to delete the ".keep" files after downloading.  The original output is currently contained in all folders, but may be archived or replaced as you see fit.

## Contact us

If you have questions or concerns regarding this code, or would like help in re-formatting it for your own use, please do not hesitate to contant the corresponding author at:

elizabeth.duskey@slu.se

## License

Copyright 2022 Elizabeth Duskey

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.