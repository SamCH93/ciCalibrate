#' @title Calibrate confidence intervals to support intervals
#'
#' @description This function TODO write documentation
#'
#' @param ci Confidence interval
#' @param ciLevel Confidence level
#' @param estimate Parameter estimate
#' @param se Standard error of the estiamte
#' @param siLevel Support level
#' @param method Calibration method
#' @param priorMean Prior mean
#' @param priorSD Prior standard deviation
#'
#' @return A supInt object
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## confidence interval of hazard ratio needs to be transformed to log-scale
#' ciHR <- c(0.75, 0.93)
#' ci <- log(ciHR)
#'
#' ## need prior under alternative hypothesis H1
#' m <- log(0.8)
#' s <- 2
#'
#' ## compute 3 support interval
#' si <- ciCalibrate(ci = ci, method = "SI-normal", priorMean = m, priorSD = s,
#'                   siLevel = 10)
#' si # on logHR scale
#' exp(si$si) # on HR scale
#'
#' ## plot Bayes factor function and support interval
#' plot(si)
#'
#' ## minimum support interval
#' msi <- ciCalibrate(ci = ci, method = "mSI-normal-local")
#' plot(si)
#'
#' @export
ciCalibrate <- function(ci = NULL,
                        ciLevel = 0.95,
                        estimate = mean(ci),
                        se = diff(ci) * 0.5 / stats::qnorm(p = 0.5*(1 + ciLevel)),
                        siLevel = 1,
                        method = c("SI-normal", "SI-normal-local", "mSI-all",
                                   "mSI-normal-local", "mSI-eplogp"),
                        priorMean,
                        priorSD) {
    ## input checks
    stopifnot(
        ## TODO write input checks for ci abd ciLevel

        length(ciLevel) == 1,
        is.numeric(ciLevel),
        is.finite(ciLevel),
        0.5 < ciLevel, ciLevel < 1,

        length(estimate) == 1,
        is.numeric(estimate),
        is.finite(estimate),

        length(se) == 1,
        is.numeric(se),
        is.finite(se),
        0 < se,

        length(siLevel) == 1,
        is.numeric(siLevel),
        is.finite(siLevel),
        0 < siLevel

    )
    method <- match.arg(method)

    if (method == "SI-normal") {
        ## input checks
        stopifnot(
            length(priorMean) == 1,
            is.numeric(priorMean),
            is.finite(priorMean),

            length(priorSD) == 1,
            is.numeric(priorSD),
            is.finite(priorSD),
            0 <= priorSD
        )
        ## standard error multiplier mSE
        mSE <- sqrt(log(1 + priorSD^2/se^2) +
                    (estimate - priorMean)^2/(se^2 + priorSD^2) -
                    2*log(siLevel))
        ## Bayes factor function bfFun
        bfFun <- function(x) {
            sqrt(1 + priorSD^2/se^2)*exp(-0.5*((estimate - x)^2/se^2 -(estimate - priorMean)^2/
                                          (se^2 + priorSD^2)))
        }
    }

    ## global normal prior under the alternative
    if (method == "SI-normal-local") {
        ## input checks
        stopifnot(
            length(priorSD) == 1,
            is.numeric(priorSD),
            is.finite(priorSD),
            0 < priorSD
        )
        ## standard error multiplier mSE
        mSE <- sqrt((log(1 + priorSD^2/se^2) - 2*log(siLevel))*(1 + se^2/priorSD^2))
        ## Bayes factor function bfFun
        bfFun <- function(x) {
            sqrt(1 + priorSD^2/se^2)*exp(-0.5*((estimate - x)^2/(1 + se^2/priorSD^2)))
        }
    }

    ## class of all priors under the alternative
    if (method == "mSI-all") {
        ## standard error multiplier mSE
        mSE <- sqrt(-2*log(siLevel))
        ## Bayes factor function bfFun
        bfFun <- function(x) {
            exp(-0.5*(estimate - x)^2/se^2)
        }
    }

    ## class of local normal priors under the alternative
    if (method == "mSI-normal-local") {
        ## standard error multiplier mSE
        mSE <- sqrt(-lamW::lambertWm1(x = -siLevel^2/exp(1)))
        ## Bayes factor function bfFun
        bfFun <- function(x) {
            z <- (estimate - x)/se
            ifelse(abs(z) < 1, 1, abs(z)*exp(-0.5*(z^2 - 1)))
        }
    }

    ## class of Beta(a, 1), a >= 1 priors for the p-value under the alternative
    if (method == "mSI-eplogp") {
        ## standard error multiplier mSE
        mSE <- stats::qnorm(p = 1 - 0.5*exp(lamW::lambertWm1(x = -siLevel/exp(1))))
        ## Bayes factor function bfFun
        bfFun <- function(x) {
            p <- 2*(1 - stats::pnorm(q = abs(estimate - x)/se))
            ifelse(p > exp(-1), 1, -exp(1)*p*log(p))
        }
    }

    ## compute support interval
    si <- estimate + c(-1, 1)*se*mSE
    if (is.nan(mSE)) {
        warning("Support interval does not exist for specified support level")
    }
    res <- list(si = si, bfFun = bfFun, estimate = estimate, se = se,
                level = siLevel, method = method)
    class(res) <- "supInt"
    return(res)
}


#' Print method for supInt object
#' @method print supInt
#' @param x A supInt object
#' @param ... Other arguments
#' @export
print.supInt <- function(x, ...) {
    print(x$si)
    invisible(x)
}

#' Plot method for supInt object
#' @method plot supInt
#' @param x A supInt object
#' @param xlim Limits of x-axis
#' @param ... Other arguments
#' @export
plot.supInt <- function(x,
                        xlim = x$estimate + c(-1, 1)*3*x$se,
                        ...) {
    ## define suitable x-axis
    xseq <- seq(xlim[1], xlim[2], length.out = 1000)

    ## plot Bayes factor function
    bf <- x$bfFun(x = xseq)
    plot(x = xseq, y = bf, log = "y", type = "l", las = 1,
         xlab = bquote("Parameter value" ~ theta[scriptstyle("0")]),
         ylab = bquote("BF"["01"]),
         main = bquote({italic(H)[0] * ":" ~ theta == theta[scriptstyle("0")]} ~
                           "vs." ~
                           italic(H)[1] * ":" ~ theta != theta[scriptstyle("0")]),
         ...)
    graphics::points(x = x$estimate, y = x$level,
                     # y = x$bfFun(x = x$estimate),
                     pch = 20)
    graphics::abline(h = 1, col = "#00000033")
    graphics::arrows(x0 = x$si[1], x1 = x$si[2], y0 = x$level, angle = 90,
                     code = 3, length = 0.1)

}
