## Extended Simple Stock Synthesis 

Extended Simple Stock Synthesis (XSSS) is an assessment method for application to data-limited stocks.  XSSS uses adaptive importance sampling to update parameter priors based upon index data.  

*Code available on this website comes with no warranty or guarantee of accuracy. It merely represents an ongoing attempt to integrate output plotting, statistics and diagnostics. It is absolutely necessary that prior to use with a new application, the user checks the output manually to verify that there are no plotting or statistical bugs which could incorrectly represent the output files being analyzed.*

## Installation

```S
install.packages("devtools")
devtools::install_github("XSSS")
```

Note: devtools may give this message: "*WARNING: Rtools is required to build R packages, but is not currently installed.*" However, Rtools is NOT required for installing r4ss via devtools, so ignore the warning.

Once you have installed the r4ss package, it can be loaded using:

```S
library(xsss)
```

## Setting up a model

## Running XSSS