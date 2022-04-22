# Weclome to hmmdiadromy

The HMM folder contains all code necessary to replicate the simulation study found in Duskey and Sullivan 2021, and referenced and altered in Duskey, Higgs, Fox, Breece, and Sullivan 2021.  We will update with the full citations when the articles have been accepted for publication.  Briefly, the code can simulate migratory paths for diadromous fish, and then run a Bayesian Hidden markov model (HMM) with nested Markov movement processes describing upriver and downriver migration.  The model also includes a logistic model of habitat preference, and a basical model of search and detection.

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

1. depend.R -- to set mypath
2. RUN_hmm.R -- to simulate data and run both models for HMM data
3. RUN_smm.R -- to simulate data and run both models for Standard Markov data
4. allplot.R -- to recreate all simulation figures in Duskey and Sullivan 2021

All files not appearing here are sourced in one or more of the files listed above, and typically contain functions necessary to run the code in each script.

NOTE: There are several nearly empty output folders contained in HMM, containing only a textfile called ".keep".  These are meant to contain and organize your own output.  
Feel free to delete the ".keep" files after downloading.  If you would like our original model output, it is available upon request.

## Contact us

If you have questions or concerns regarding this code, or would like help in re-formatting it for your own use, please do not hesitate to contant the corresponding author at:

elizabeth.duskey@slu.se

## License

None.