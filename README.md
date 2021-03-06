# fwerapprox

R package for

(1) computing score test statistics for testing whether each of a large number of coefficients (typically corresponding to genetic markers) in a GLM is zero in presence of a smaller number of covariates (typically environmental covariates), and estimation of correlations between the test statistics, and

(2) computing Glaz–Johnson-type intersection approximations of the asymptotic multivariate normal distribution of the test statistics, which can be used for FWER control.

To install:

```r
library(devtools)
install_github("oyvind-bakke/fwerapprox")
```

To load:

```r
library(fwerapprox)
```

Basic information:

```r
package?fwerapprox
```
