#' @title Calibrate confidence intervals to support intervals
#'
#' @description This function computes a support interval for an unknown
#'     parameter based on either a confidence interval for the parameter or a
#'     parameter estimate with standard error.
#'
#' @details A *support interval* with support level \eqn{k} is defined by the
#'     parameter values \eqn{\theta_0}{theta0} for which the Bayes factor
#'     contrasting \eqn{H_0\colon \theta = \theta_0}{H0: theta = theta0} to
#'     \eqn{H_1\colon \theta \neq \theta_0}{H1: theta != theta0} is larger or
#'     equal than \eqn{k}, i.e., the parameter values for which the data are at
#'     least \eqn{k} times more likely than under the alternative. Different
#'     prior distributions for the parameter \eqn{\theta}{} under the
#'     alternative \eqn{H_1}{H1} are available:
#'
#' * \code{method = "SI-normal"}: a normal prior centered around
#' \code{priorMean} with standard deviation \code{priorSD}, i.e., \eqn{\theta
#' \,|\, H_1 \sim N(\code{priorMean}, \code{priorSD}^2)}{theta | H1 ~
#' N(priorMean, priorSD^2)}
#'
#' * \code{method = "SI-normal-local"}: a local normal prior with standard
#' deviation \code{priorSD}, i.e., \eqn{\theta \,|\, H_1 \sim N(\theta_0,
#' \code{priorSD}^2)}{theta | H1 ~ N(theta0, priorSD^2)}
#'
#' * \code{method = "SI-normal-nonlocal"}: a nonlocal normal moment prior with
#' spread parameter \code{priorSD}, i.e., a prior with density \eqn{f(\theta
#' \,|\, H_1) = N(\theta \,|\, \theta_0, \code{priorSD}^2) \times (\theta -
#' \theta_0)^2/\code{priorSD}^2}{f(theta | H1) = N(theta | theta0, priorSD^2)*
#' (theta - theta0)^2/\code{priorSD}^2}
#'
#' The function also allows to compute *minimum support intervals* which require
#' to only specify a class of priors for the parameter under the alternative and
#' then compute the minimum Bayes factor over the class of alternatives. The
#' following classes of prior distribution are available:
#'
#' * \code{method = "mSI-all"}: the class of all prior distributions under the
#' alternative, this leads to the narrowest support interval possible
#'
#' * \code{method = "mSI-normal-local"}: the class of local normal prior
#' distributions under the alternative, i.e., \eqn{\theta \,|\, H_1 \sim
#' N(\theta_0, v)}{theta | H1 ~ N(theta0, v)} with \eqn{v \geq 0}{v >= 0}
#'
#' * \code{method = "mSI-eplogp"}: the class of monotonically decreasing beta
#' prior distributions on the p-value of the data \eqn{p = 2(1 -
#' \Phi(|\code{estimate} - \theta_0|/\code{se}))}{p = 2*(1 - pnorm(abs(estimate
#' - theta0)/se))}, i.e. \eqn{p \,|\, H_1 \sim \mbox{Be}(\xi, 1)}{p | H1 ~
#' Be(xi, 1)} with \eqn{\xi \geq 1}{xi >= 1}
#'
#' @md
#'
#' @param ci Confidence interval given as a numeric vector of length two.
#' @param ciLevel Confidence level. Defaults to 0.95.
#' @param estimate Parameter estimate. Only required if no confidence interval
#'     and confidence level are specified.
#' @param se Standard error of the parameter estimate. Only required if no
#'     confidence interval and confidence level are specified.
#' @param siLevel Support level. Defaults to 1.
#' @param method Calibration method. Can either be \code{"SI-normal"},
#'     \code{"SI-normal-local"}, \code{"SI-normal-nonlocal"}, \code{"mSI-all"},
#'     \code{"mSI-normal-local"}, or \code{"mSI-eplogp"}. Defaults to
#'     \code{"SI-normal"}. See details for more information.
#' @param priorMean Prior mean, only required for \code{"SI-normal"}.
#' @param priorSD Prior standard deviation / spread, only required for
#'     \code{"SI-normal"}, \code{"SI-normal-local"}, \code{"SI-normal-nonlocal"}.
#'
#' @return Returns an object of class \code{"supInt"} which is a list
#'  containing:\tabular{ll}{
#'    \code{si} \tab The computed support interval. \cr
#'    \tab \cr
#'    \code{bfFun} \tab The computed Bayes factor function. \cr
#'    \tab \cr
#'    \code{estimate} \tab The specified parameter estimate. \cr
#'    \tab \cr
#'    \code{se} \tab The specified standard error. \cr
#'    \tab \cr
#'    \code{siLevel} \tab The specified support level. \cr
#'    \tab \cr
#'    \code{ciLevel} \tab The specified confidence level. \cr
#'    \tab \cr
#'    \code{priorParams} \tab The specified prior parameters. \cr
#' }
#'
#' @references
#'
#' Pawel, S., Ly, A., and Wagenmakers, E.-J. (2023). Evidential Calibration of
#' Confidence Intervals. The American Statistician.
#' \doi{10.1080/00031305.2023.2216239}
#'
#' Wagenmakers, E.-J., Gronau, Q. F., Dablander, F., and Etz, A. (2020). The
#' support interval. Erkenntnis. \doi{10.1007/s10670-019-00209-z}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## confidence interval of hazard ratio needs to be transformed to log-scale
#' ciHR <- c(0.75, 0.93)
#' ci <- log(ciHR)
#'
#' ## normal prior under the alternative hypothesis H1
#' m <- log(0.8) # prior mean
#' s <- 2 # prior sd
#'
#' ## compute 10 support interval
#' si <- ciCalibrate(ci = ci, method = "SI-normal", priorMean = m,
#'                   priorSD = s, siLevel = 10)
#' si # on logHR scale
#' exp(si$si) # on HR scale
#'
#' ## plot Bayes factor function and support interval
#' plot(si)
#'
#' ## minimum support interval based on local normal priors
#' msi <- ciCalibrate(ci = ci, method = "mSI-normal-local")
#' plot(msi)
#'
#' @export
ciCalibrate <- function(ci = NULL,
                        ciLevel = 0.95,
                        estimate = mean(ci),
                        se = diff(ci) * 0.5 / stats::qnorm(p = 0.5*(1 + ciLevel)),
                        siLevel = 1,
                        method = c("SI-normal", "SI-normal-local", "SI-normal-nonlocal",
                                   "mSI-all", "mSI-normal-local", "mSI-eplogp"),
                        priorMean,
                        priorSD) {
    ## input checks
    stopifnot(
        is.null(ci) |
        ((length(ci) == 2) &
         all(is.numeric(ci)) &
         all(is.finite(ci))),

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

    ## global normal prior under the alternative
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
        ## standard error multiplier
        mSE <- sqrt(log(1 + priorSD^2/se^2) +
                    (estimate - priorMean)^2/(se^2 + priorSD^2) -
                    2*log(siLevel))
        ## Bayes factor function
        bfFun <- function(x) {
            sqrt(1 + priorSD^2/se^2)*exp(-0.5*((estimate - x)^2/se^2 -(estimate - priorMean)^2/
                                          (se^2 + priorSD^2)))
        }
        ## prior parameters
        priorParams <- list(type = method, priorMean = priorMean, priorSD = priorSD)
    }

    ## local normal prior under the alternative
    if (method == "SI-normal-local") {
        ## input checks
        stopifnot(
            length(priorSD) == 1,
            is.numeric(priorSD),
            is.finite(priorSD),
            0 < priorSD
        )
        ## standard error multiplier
        mSE <- sqrt((log(1 + priorSD^2/se^2) - 2*log(siLevel))*(1 + se^2/priorSD^2))
        ## Bayes factor function
        bfFun <- function(x) {
            sqrt(1 + priorSD^2/se^2)*exp(-0.5*(estimate - x)^2/se^2/(1 + se^2/priorSD^2))
        }
        ## prior parameters
        priorParams <- list(type = method, priorSD = priorSD)
    }

    ## non-local normal moment prior under the alternative
    if (method == "SI-normal-nonlocal") {
        ## input checks
        stopifnot(
            length(priorSD) == 1,
            is.numeric(priorSD),
            is.finite(priorSD),
            0 < priorSD
        )
        ## standard error multiplier
        mSE <- sqrt((2*lamW::lambertW0(x = 0.5*exp(0.5)/siLevel*(1 + priorSD^2/se^2)^1.5) - 1)*
            (1 + se^2/priorSD^2))
        ## Bayes factor function
        bfFun <- function(x) {
            (1 + priorSD^2/se^2)^1.5 * exp(-0.5 * (estimate - x)^2/se^2 / (1 + se^2/priorSD^2)) /
                (1 + (estimate - x)^2/se^2 / (1 + se^2/priorSD^2))
        }
        ## prior parameters
        priorParams <- list(type = method, priorSD = priorSD)
    }

    ## class of all priors under the alternative
    if (method == "mSI-all") {
        ## standard error multiplier mSE
        mSE <- sqrt(-2*log(siLevel))
        ## Bayes factor function bfFun
        bfFun <- function(x) {
            exp(-0.5*(estimate - x)^2/se^2)
        }
        ## prior parameters
        priorParams <- list(type = method)
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
        ## prior parameters
        priorParams <- list(type = method)
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
        ## prior parameters
        priorParams <- list(type = method)
    }

    ## compute support interval
    si <- estimate + c(-1, 1)*se*mSE
    if (is.nan(mSE)) {
        warning("Support interval does not exist for specified support level")
    }
    res <- list(si = si, bfFun = bfFun, estimate = estimate, se = se,
                siLevel = siLevel, ciLevel = ciLevel, priorParams = priorParams)
    class(res) <- "supInt"
    return(res)
}


#' Print method for class \code{"supInt"}
#' @method print supInt
#'
#' @description Prints parameter estimate, confidence interval, and support
#'     interval.
#'
#' @param x Object of class \code{"supInt"}
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"supInt"} object
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## compute and print 5 support interval
#' si5 <- ciCalibrate(ci = c(-2, -0.5), method = "SI-normal", priorMean = -0.3,
#'                    priorSD = 2, siLevel = 5)
#' print(si5)
#'
#' @export
print.supInt <- function(x, ...) {
    ## Point estimate with confidence interval
    ci <- x$estimate + c(-1, 1)*x$se*stats::qnorm(p = 0.5*(1 + x$ciLevel))
    cat(paste0("\nPoint Estimate [",
               x$ciLevel*100, "% Confidence Interval] \n"))
    cat(signif(x$estimate, 2),
        " [", paste0(signif(ci, 2), collapse = ","), "]\n", sep = "")

    ## Calibration method
    cat("\nCalibration Method\n")
    if (x$priorParams$type == "SI-normal") {
        cat("Normal prior for parameter under alternative\n",
            "with mean m = ", signif(x$priorParams$priorMean, 2),
            " and standard deviation sd = ", signif(x$priorParams$priorSD, 2), sep = "")
    }
    ## local normal prior under the alternative
    if (x$priorParams$type == "SI-normal-local") {
        cat("Local normal prior for parameter under alternative\n",
            "with standard deviation sd = ", signif(x$priorParams$priorSD, 2),
            sep = "")
    }
    ## nonlical normal moment prior under the alternative
    if (x$priorParams$type == "SI-normal-nonlocal") {
        cat("Nonlocal normal moment prior for parameter under alternative\n",
            "with spread s = ", signif(x$priorParams$priorSD, 2),
            sep = "")
    }
    ## class of all priors under the alternative
    if (x$priorParams$type == "mSI-all") {
        cat("Minimizing support for class of all priors for parameter\n",
            "under alternative", sep = "")
    }
    ## class of local normal priors under the alternative
    if (x$priorParams$type == "mSI-normal-local") {
        cat("Minimizing support for class of local normal priors for parameter\n",
            "under alternative", sep = "")
    }
    ## class of Beta(a, 1), a >= 1 priors for the p-value under the alternative
    if (x$priorParams$type == "mSI-eplogp") {
        cat("Minimizing support for the class of Beta(a, 1), a >= 1 priors\n",
            "for p-value of data under alternative", sep = "")
    }
    cat("\n")

    ## Support interval
    if (x$siLevel < 1) siLevel <- paste0("1/", signif(1/x$siLevel, 2))
    else siLevel <- as.character(signif(x$siLevel, 2))
    if (x$priorParams$type %in% c("mSI-eplogp", "mSI-normal-local", "mSI-all")) {
        siString <- " Minimum Support Interval\n"
    } else {
        siString <- " Support Interval\n"
    }
    cat(paste0("\nk = ", siLevel, siString))
    cat("[", paste0(signif(x$si, 2), collapse = ","), "]\n", sep = "")

    invisible(x)
}

#' Plot method for class \code{"supInt"}
#' @method plot supInt
#'
#' @description Plots Bayes factor function and support interval at the specified support level.
#'
#' @param x Object of class \code{"supInt"}
#' @param xlim Limits of x-axis
#' @param ... Other arguments passed to \code{plot}
#'
#' @return Generates a plot of the Bayes factor function with support interval
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## compute 3 support interval
#' si3 <- ciCalibrate(ci = c(0.3, 1.5), method = "SI-normal", priorMean = 0,
#'                    priorSD = 2, siLevel = 3)
#'
#' ## plot Bayes factor function and support interval
#' plot(si3, xlim = c(0, 1.75))
#'
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
    graphics::points(x = x$estimate, y = x$siLevel,
                     # y = x$bfFun(x = x$estimate),
                     pch = 20)
    graphics::abline(h = 1, col = "#00000033")
    graphics::arrows(x0 = x$si[1], x1 = x$si[2], y0 = x$siLevel, angle = 90,
                     code = 3, length = 0.1)

}
