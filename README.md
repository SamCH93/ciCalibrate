# ciCalibrate

**ciCalibrate** is an R package for computing support intervals for unknown
univariate parameters. A support interval can either be computed based on a
parameter estimate and standard error or based on a confidence interval for the
respective parameter. The main function for doing so is `ciCalibrate`, see the
manual with `?ciCalibrate` for the available options. Theoretical background on
support intervals is provided in the accompanying preprint 
[Pawel et al. (2022)](https://doi.org/10.48550/arXiv.2206.12290) and also 
[Wagenmakers et al.(2020)](https://doi.org/10.1007/s10670-019-00209-z).

## Installation

```r
remotes::install_github(repo = "SamCH93/ciCalibrate")
```

## Usage

``` r
library("ciCalibrate")

## data from RECOVERY trial
logHR <- -0.19 # estimate
se <- 0.05 # standard error of estimate
ci95 <- logHR + c(-1, 1) * qnorm(p = 0.975) * se # 95% Wald-CI

## default normal prior for logHR under the alternative H1
pm <- 0 # center around value of no effect
psd <- 2 # unit-information standard deviation for a logHR

## compute a support interval with support level = 10
si10 <- ciCalibrate(estimate = logHR, se = se, siLevel = 10, method = "SI-normal",
                    priorMean = pm, priorSD = psd)

## compute instead with confidence interval as input
si10 <- ciCalibrate(ci = ci95, ciLevel = 0.95, siLevel = 10, method = "SI-normal",
                    priorMean = pm, priorSD = psd)
si10

#> Point Estimate [95% Confidence Interval] 
#> -0.19 [-0.29,-0.092]
#> 
#> Calibration Method
#> Normal prior for parameter under alternative
#> with mean m = 0 and standard deviation sd = 2
#> 
#> k = 10 Support Interval
#> [-0.27,-0.11]

## plot Bayes factor function and support interval
plot(si10)
```
![Output of the command plot(si10): the Bayes factor function and the 10 support
interval](SIexample.png)

## References

* Pawel, S., Ly, A., and Wagenmakers, E.-J. (2022). Evidential Calibration of
  Confidence Intervals.
  [doi:10.48550/arXiv.2206.12290](https://doi.org/10.48550/arXiv.2206.12290)

* Wagenmakers, E.-J., Gronau, Q. F., Dablander, F., and Etz, A. (2020). The
  support interval. Erkenntnis.
  [10.1007/s10670-019-00209-z](https://doi.org/10.1007/s10670-019-00209-z)
