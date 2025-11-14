# ecofowlcasting <img src="https://www.alexbaecher.com/project/disease-ecology/featured_hu3854238940016619547.webp" alt="ecofowlcasting tool for predicting the transmission of vector borne diseases in a sentinel chicken system" align="right" height="200"/>

> A repository for predicting the transmission of vector borne diseases in a sentinel chicken system

### Author: [J. Alex Baecher](https://www.alexbaecher.com/)

### Correspondance: please contact the first author, J. Alex Baecher, for questions about the code or data: 
  - e-mail: jbaecher@gmail.com
  - github: https://github.com/slamander
  - personal website: https://www.alexbaecher.com/

#### Citation:
[J. Alex Baecher](https://www.alexbaecher.com/), V. A. Akshay, [Robert P. Guralnick](https://www.gurlab.net/), Amy M. Bauer, Yasmin N. Tavares, Yesenia Sánchez, [James T. Thorson](https://sites.google.com/site/thorsonresearch/), [Lindsay P. Campbell](https://lcampbelllab.wixsite.com/campbell-lab/). [Toward Ecological Forecasting of West Nile Virus in Florida: Insights from Two Decades of Sentinel Chicken Surveillance](https://doi.org/10.1016/j.scitotenv.2025.180308). *Science of the Total Environment* 
__________________________________________________________________________________________________________________________________________

## Repository Directory

### [`./functions`](./functions): Contains helper functions for executing code contained in [scripts](./scripts)
  - [`buildmermod_to_glmmtmb.R`](./functions/buildmermod_to_glmmtmb.R): Converting `buildmermod` objects to `glmmTMB` objects
  - [`compare_models.R`](./functions/compare_models.R): Automatically fits `glmmTMB` models and performs AIC selection
  - [`extract_model_data.R`](./functions/extract_model_data.R): Custom methods for extracting parameter estimates (excluding unwanted parameters)
  - [`stepwise_vif.R`](./functions/stepwise_vif.R): Iteratively performs VIF calculations while removing multicollinear variables
  - [`truncate.R`](./functions/truncate.R): Truncates raster values based on supplied threshold values to improve visualization

### [`./scripts`](./scripts): Contains code for modeling spatiotemporal transmission dynamics of West Nile virus
  - [`/data_processing`](./scripts/data_processing): Contains scripts for processing WNV monitoring data
  - [`/environmental_data`](./scripts/environmental_data): Contains scripts for assembling environmental predictor data
  - [`/models`](./scripts/models): Contains scripts for executing temporal and spatiotemporal models
    - [`/glmmtmb`](./scripts/models/glmmtmb): Contains script for temporal modeling using [glmmTMB](https://github.com/glmmTMB/glmmTMB)
      - [`base_model_development.R`](./scripts/models/glmmtmb/base_model_development.R): Script for model selection and variable reduction
    - [`/sdmtmb`](./scripts/models/sdmtmb): Contains scripts for spatiotemporal modeling with [sdmTMB](https://pbs-assess.github.io/sdmTMB/)
      - [`/fitting`](./scripts/models/sdmtmb/fitting): Contains scripts for model calibration and selection
        - [`model_calibration.R`](./scripts/models/sdmtmb/fitting/model_calibration.R): Fitting `sdmTMB` model
        - [`model_selection.R`](./scripts/models/sdmtmb/fitting/model_selection.R): Comparing model metrics for **SI Table 1** and **SI Table 2**
      - [`/prediction`](./scripts/models/sdmtmb/prediction): Contains scripts for obtaining predictions
        - [`calculate_rmse.R`](./scripts/models/sdmtmb/fitting/calculate_rmse.R): Calculate rmse values for manuscript **Table 1**
        - [`point_based.R`](./scripts/models/sdmtmb/fitting/point_based.R): Comparing epsilon-correction and plotting for **Figure 2** 
        - [`response_curves.R`](./scripts/models/sdmtmb/fitting/response_curves.R): Obtaining predictions for response curves in **Figure 3**
        - [`statewide_preds.R`](./scripts/models/sdmtmb/fitting/statewide_preds.R): Obtaining predictions for statewide plots in **Figure 4** and **Figure 5**

## Data
Georeferenced sentinel chicken seroconversion data is available upon request through the Florida Department of Health Arbovirus Surveillance program upon agreement from participating Florida mosquito control programs through a memorandum of understanding. The authors did not receive special privileges in accessing the data that other researchers would not have. Contact information for data requests are available through the [FDOH website](https://www.floridahealth.gov/diseases-and-conditions/mosquito-borne-diseases/surveillance.html). 
