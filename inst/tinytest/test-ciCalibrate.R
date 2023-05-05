library(tinytest)
library(ciCalibrate)

## set up testing grid
k <- c(1/10, 1, 10)
x <- 1
se <- 0.2
mp <- 0.5
sep <- 2
siMethods <- c("SI-normal", "SI-normal-local", "SI-normal-nonlocal")
msiMethods <- c("mSI-all", "mSI-normal-local", "mSI-eplogp")
grid1 <- expand.grid(method = siMethods, k = k, stringsAsFactors = FALSE)
grid2 <- expand.grid(method = msiMethods, k = k[k <= 1], stringsAsFactors = FALSE)
grid <- rbind(grid1, grid2)

## recompute bf from support interval
resultsDF <- do.call("rbind", lapply(X = seq(1, nrow(grid)), FUN = function(i) {
    method <- grid$method[i]
    siLevel <- grid$k[i]
    siObject <- ciCalibrate(estimate = x, se = se, method = method,
                            siLevel = siLevel, priorMean = mp, priorSD = sep)
    si <- siObject$si
    bf <- siObject$bfFun(x = si)
    out <- data.frame(method = method, k = siLevel, lower = si[1], upper = si[2],
                      klower = bf[1], kupper = bf[2])
    return(out)
}))

## check that recomputed bf corresponds to input support level
expect_equal(cbind(resultsDF$k, resultsDF$k),
             cbind(resultsDF$klower, resultsDF$kupper),
             info = "Checking that BFs and (minimum) SIs correspond")
