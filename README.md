## Extended Simple Stock Synthesis 

Extended Simple Stock Synthesis (XSSS) is an assessment method for application to data-limited stocks.  XSSS uses adaptive importance sampling to update parameter priors based upon index data.  

*Code available on this website comes with no warranty or guarantee of accuracy.*

## Installation

```S
install.packages("devtools")
devtools::install_github("XSSS")
```

Note: devtools may give this message: "*WARNING: Rtools is required to build R packages, but is not currently installed.*" However, Rtools is NOT required for installing xsss via devtools, so ignore the warning.

Once you have installed the xsss package, it can be loaded using:

```S
library(xsss)
```

## Setting up a model

XSSS is data-limited assessment approach for application when only catches and indices of abundance are available allowing for multiple fishing and survey fleets. XSSS requires four input files used by Stock Synthesis: starter, forecast, data, and control file.  Biological relationships (e.g., weight-at-length, fecundity, maturity) and selectivity forms must be defined in the control file.  Example model files are located within the "example" folder.       

## Running XSSS

Example model files and calls to the main function are included in the "example" folder. The user must define distributions for three parameters: natural mortality, steepness, and relative stock status is a specific year.

## Description of package

Wetzel, C.R., and Punt, A.E. 2015. Evaluating the performance of data-moderate and catch-only assessment methods for U.S. west coast groundfish. Fisheries Research, 171: 170-187.  https://www.sciencedirect.com/science/article/pii/S016578361500185X