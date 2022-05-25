# ciCalibrate

`ciCalibrate` is an R package for computing support intervals for unknown
univariate parameters, based on either a parameter estimate and standard error
or from a confidence interval for the respective parameter. The main function
for doing so is `ciCalibrate`, see the manual for the available options.

## Installation

```r
remotes::install_github(repo = "SamCH93/ciCalibrate")
```

## Usage

``` r
library(ciCalibrate)

## data from RECOVERY trial
logHR <- -0.19 # estimate
se <- 0.05 # standard error of estimate
ci95 <- logHR + c(-1, 1) * qnorm(p = 0.975) * se # 95% Wald-CI

## default normal prior for logHR under the alternative
pm <- 0 # center around value of no effect
psd <- 2 # unit-information standard deviation for a logHR

## compute a support interval with support level = 3
si3 <- ciCalibrate(estimate = logHR, se = se, siLevel = 3, method = "SI-normal",
                   priorMean = pm, priorSD = psd)
si3

## compute instead with confidence interval as input
si3 <- ciCalibrate(ci = ci95, ciLevel = 0.95, siLevel = 3, method = "SI-normal",
                   priorMean = pm, priorSD = psd)
si3

## plot Bayes factor function and support interval
plot(si3)
```

## References

* Pawel, S., Ly, A., and Wagenmakers, E.-J. (2022). Evidential calibration of
  confidence intervals. to appear on arXiv soon

* Wagenmakers, E.-J., Gronau, Q. F., Dablander, F., and Etz, A. (2020). The
  support interval. Erkenntnis.
  [10.1007/s10670-019-00209-z](https://doi.org/10.1007/s10670-019-00209-z)
