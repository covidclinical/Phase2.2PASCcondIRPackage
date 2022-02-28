PASC Conditional Independence Testing
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

This package implements conditional independence testing procedures to
identify codified concepts related to prior exposure to SARS-CoV-2 in
4CE Phase2.2 Data. Running this package necessitates refreshed Phase2.2 
data to Feb 1, 2022 inclusive of all patient cohorts and membership in 
the 4CE Consortium.

4CE analysts at each participating healthcare system can follow the
steps described bellow to run the analysis locally, using their
healthcare system’s 2.2 data. The summarized results are then shared and
combined with those obtained from other healthcare systems. Those
phenotype conditional testing results shared include PheCode, beta coefficient,
standard errors, p-values, and prevalence estimates.

### 1. Install the package

Install the development version of **Phase2.2PASCcondIRPackage** from
GitHub using remotes by running the code in your R console below:

``` r
remotes::install_github('covidclinical/Phase2.2PASCcondIRPackage',
                        upgrade = FALSE)
library(Phase2.2PASCcondIRPackage)
```

### 2. Run the Analysis

Please note that the analysis may take 2-3 days to
complete. To start the analysis, define the following 3 variables below
in your R console. "siteid" must be in all capital letters. 

``` r
siteid=""     # specify the site ID (capital letters)
dir.data=""   # specify the 2.2 data directory
dir.repo=""   # specify the directory to save the results 

# read the data
obs = fread(paste0(dir.data,"Phase22all_LocalPatientObservations.csv"),stringsAsFactors = F)
summary = fread(paste0(dir.data,"Phase22all_LocalPatientSummary.csv"),stringsAsFactors = F)

# run the analysis
runAnalysis()
```

### 3. Submit the results

<!-- #### 1. Submit the results via GitHub -->
<!-- Finally, please submit the results to [Phase2.2PASCcondISummariesPublic](https://github.com/covidclinical/Phase2.2PASCcondISummariesPublic): -->
<!-- - Share your GitHub handle with @Clara-Lea via direct message so you can be added as contributor to the repository. -->
<!-- - Note that you will need to use a token to access **private** repositories. Please see the details [here](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token). -->
<!-- To generate a new token, go to your GitHub settings -> Developer settings -> Personal access tokens -> Generate. -->
<!-- Then, run the following: -->
<!-- ```{r, eval=FALSE} -->
<!-- submitAnalysis() -->
<!-- ``` -->
<!-- #### 2. Alternatively, you can submit the results via Slack.  -->

Your results will be saved under the path that you specify in
`dir.repo`, as an .Rdata file. Please share the results with @Harrison
and @Clara-Lea via the Slack channel #post-acute-sequelae.

Thank you very much for your contribution!
