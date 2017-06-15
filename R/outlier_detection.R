#' @title newfunction1
#' @description newfunction1
#' @export
#' @import tsoutliers
outlier_detection<- function(y, xreg = NULL, cval = NULL, delta = 0.7, n.start = 50,
                             types = c("AO", "LS", "TC"), # c("IO", "AO", "LS", "TC", "SLS")
                             maxit = 1, maxit.iloop = 4, cval.reduce = 0.14286,
                             remove.method = c("en-masse", "bottom-up"),
                             remove.cval = NULL,
                             tsmethod = c("auto.arima", "arima", "stsm"),
                             args.tsmethod = NULL, args.tsmodel = NULL, logfile = NULL)
{
  tsmethod <- match.arg(tsmethod)
  remove.method <- match.arg(remove.method)
  attr.y <- attributes(y)
  n <- length(y)
  yname <- deparse(substitute(y))
  #stopifnot(is.ts(y))

  if (!is.null(args.tsmethod$xreg))
  {
    xreg <- args.tsmethod$xreg
    args.tsmethod$xreg <- NULL # this removes element "xreg" from the list

  }

  if (tsmethod == "stsm")
  {
    if (is.null(args.tsmodel$model))
      args.tsmodel$model <- ifelse(frequency(y) == 1, "local-level", "BSM")

    ##FIXME these defaults only if stsm.method = "maxlik.fd.scoring"

    if (is.null(args.tsmodel$ssd))
      args.tsmodel$ssd <- TRUE
    if (is.null(args.tsmodel$sgfc))
      args.tsmodel$sgfc <- TRUE
    # let "stsm::stsmFit" handle "xreg", not here
    y <- do.call("stsm.model", args = c(list(y = y), args.tsmodel))
    #ylist <- list(m = m)
  }
  if (is.null(args.tsmethod))
  {
    args.tsmethod <- switch(tsmethod,
                            "auto.arima" = list(allowdrift = FALSE, ic = "bic"),
                            "arima" = list(order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1))),
                            "stsm" = list(stsm.method = "maxlik.td.optim", method = "L-BFGS-B",
                                          KF.version = "KFKSDS", KF.args = list(P0cov = TRUE), gr = "numerical")) #hessian = TRUE
    #list(stsm.method = "maxlik.fd.scoring", step = NULL, information = "expected"))
  }

  # default critical value
  # the same is done in functions "locate.outliers.oloop" and "remove.outliers"
  # "cval" is passed as a non-null value from tso() to those functions
  # but keep there this block so that default value is used when those functions
  # are called outside tso()

  if (is.null(cval))
  {
    #n <- length(y)
    if (n <= 50) {
      cval <- 3
    } else
      if (n >= 450) {
        cval <- 4
      } else
        cval <- round(3 + 0.0025 * (n - 50), 2)
  }

  cval0 <- cval
  if (is.null(remove.cval))
    remove.cval <- cval

  # "res0" is used below to generate the output,
  # "res" is overwritten until no more outliers are found
  # "res0" is also used if maxit = 1

  res0 <- res <- tso0(x = y, xreg = xreg, cval = cval,
                      delta = delta, n.start = n.start,
                      types = types, maxit.iloop = maxit.iloop,
                      remove.method = remove.method, remove.cval = remove.cval,
                      tsmethod = tsmethod, args.tsmethod = args.tsmethod,
                      logfile = logfile)

  fit.wo.outliers <- res$fit0 # model without outliers (if maxit>1 res0 may change)
  moall <- res$outliers
  outtimes <- res$times

  iter <- 1
  cval <- round(cval * (1 - cval.reduce), 2)

  if (nrow(moall) > 1)
    while (iter < maxit)
    {
      ##FIXME see move res0 <- res after if(...) break

      if (tsmethod == "stsm")
      {
        ##FIXME TODO create stsm object based on res$yadj as done above
        warning("currently ", sQuote("maxit"), " > 1 is not allowed for ", sQuote("tsmethod=\"stsm\""))
        break
      }
      # save "res" to have a copy of the last fitted model, res$fit;
      # if in the current run no outliers are found then
      # tso0() does not return the fitted model

      res0 <- res

      res <- tso0(x = res$yadj, xreg = xreg, cval = cval,
                  delta = delta, n.start = n.start,
                  types = types, maxit.iloop = maxit.iloop,
                  remove.method = remove.method, remove.cval = remove.cval,
                  tsmethod = tsmethod, args.tsmethod = args.tsmethod,
                  logfile = logfile)

      ##FIXME check
      #discard (remove) duplicates and outliers at consecutive type points (if any)
      #
      #do not discard according to abs(t-stat) because the detection of outliers
      #are based on res$yadj (not the original series); discarding an outlier
      #from a previous iteration would require changing the current res$yadj
      #
      #discard outliers at an observation where an outlier (of the same or other type)
      #was detected in a previous iteration
      id <- which(res$outliers[,"ind"] %in% res0$outliers[,"ind"])
      if (length(id) > 0)
        res$outliers <- res$outliers[id,]
      #discard consecutive outliers of any type, keep the outlier from previous iterations
      id <- which(apply(outer(res$outliers[,"ind"], res0$outliers[,"ind"], "-"), MARGIN=1,
                        FUN = function(x) any(x == 1)))
      if (length(id) > 0)
        res$outliers <- res$outliers[id,]

      if (nrow(res$outliers) == 0)
        break

      moall <- rbind(moall, res$outliers)
      outtimes <- c(outtimes, res$times)

      iter <- iter + 1
    }

  # final model given the detected outliers
  ######################################################################
  if (nrow(moall) > 0)
  {
    #NOTE 'pars' is relevant only for innovational outliers,
    #when 'maxit'>1, see if it would be better to use 'res' instead of 'res0',
    #preferably it should be based on 'pars' from a model for the original data
    #rather than the series adjusted for outliers

    pars <- switch(tsmethod,
                   "auto.arima" = , "arima" = coefs2poly(res0$fit),
                   "stsm" = stsm::char2numeric(res0$fit$model))

    # 'xreg': input regressor variables such as calendar effects (if any)
    # 'xreg.outl': outliers regressor variables detected above (if any)
    # 'xregall': all regressors ('xreg' and 'xreg.outl')

    xreg.outl <- outliers.effects(mo = moall, n = n, weights = FALSE, delta = delta,
                                  pars = pars, n.start = n.start, freq = frequency(y))
    xregall <- cbind(xreg, xreg.outl)
    nms.outl <- colnames(xreg.outl)
    colnames(xregall) <- c(colnames(xreg), nms.outl)

    ##NOTE
    # rerunning "auto.arima" (model selection) may not be necessary at this point

    if (tsmethod == "stsm") {
      fit <- do.call("stsmFit", args = c(list(x = y, xreg = xregall), args.tsmethod))
    } else {
      fit <- do.call(tsmethod, args = c(list(x = y, xreg = xregall), args.tsmethod))
      # this is for proper printing of results from "auto.arima" and "arima"
      fit$series <- yname
    }

    id <- colnames(xreg.outl)
    if (tsmethod == "stsm")
    {
      xregcoefs <- fit$pars[id]
      tstats <- xregcoefs / fit$std.errors[id]
    } else { # method "auto.arima", "arima"
      xregcoefs <- coef(fit)[id]
      tstats <- xregcoefs / sqrt(diag(fit$var.coef)[id])
    }

    moall[,"coefhat"] <- xregcoefs
    moall[,"tstat"] <- tstats




  }
  moall
}

tso0 <- function(x, xreg = NULL, cval = 3.5, delta = 0.7, n.start = 50,
                 types = c("AO", "LS", "TC"), maxit.iloop = 4,
                 remove.method = c("en-masse", "bottom-up"),
                 remove.cval = NULL,
                 tsmethod = c("auto.arima", "arima", "stsm"), args.tsmethod = NULL,
                 args.tsmodel = NULL, logfile = NULL)
{
  # "x" can be either a "ts" object or a "stsm" object;
  # if !inherits(x, "stsm") then two identical objects are stored ("x" and "y")

  y <- if(is.ts(x)) { x } else x@y

  #remove.method <- match.arg(remove.method)
  #tsmethod <- match.arg(tsmethod)
  #remove.method <- match.arg(remove.method)
  fitmethod <- gsub("stsm", "stsmFit", tsmethod)

  if (is.null(remove.cval))
    remove.cval <- cval

  # fit time series model

  fit.wo.outliers <-
    fit <- do.call(fitmethod, args = c(list(x = x, xreg = xreg), args.tsmethod))
  #fit$series <- deparse(substitute(y))

  if (!is.null(logfile))
  {
    cat(paste("model selection:\n"), file = logfile, append = FALSE)
    capture.output(fit, file = logfile, append = TRUE)
  }

  # identify and locate prospective outliers by type
  # given a fitted time series model

  stage1 <- locate.outliers.oloop(y = y, fit = fit, types = types, cval = cval,
                                  maxit.iloop = maxit.iloop, delta = delta, n.start = n.start, logfile = logfile)

  # choose and fit the model including the outlier regressors detected so far
  # (the weights of the outliers is fine tuned, to see it
  # compare 'moall[,"coefhat"]' with 'coef(fit)["oeffi"]') then
  # remove the outliers detected so far if they are not significant in the new model/fit

  if (nrow(stage1$outliers) > 0)
  {
    stage2 <- remove.outliers(x = stage1, y = y, cval = remove.cval,
                              method = remove.method, delta = delta, n.start = n.start,
                              tsmethod.call = fit$call, fdiff = NULL, logfile = logfile)

    #moall <- stage2$outliers
    stopifnot(ncol(stage2$xreg) == length(stage2$xregcoefs))
  } else
    stage2 <- list(xreg = NULL, fit = stage1$fit)

  # final outliers and
  # original series adjusted for the outlier effects

  if (!is.null(stage2$xreg))
  {
    # stage2$fit$xreg is not returned by arima()
    moall <- stage2$outliers
    ##NOTE changed 2016Nov12 after changes in remove.outliers(), "moall" is updated there
    #moall[,"coefhat"] <- stage2$xregcoefs
    #moall[,"tstat"] <- stage2$xregtstats

    oeff <- stage2$xreg %*% cbind(stage2$xregcoefs)
    attributes(oeff) <- attributes(y)
    yadj <- y - oeff

    moall <- moall[,c("type", "ind", "coefhat", "tstat")]
    outtimes <- time(y)[moall[,"ind"]]
    if (frequency(y) > 1)
      outseason <- formatC(as.vector(cycle(y)[moall[,"ind"]]),
                           width = 2, flag="0")

    moall <- cbind(moall[,c("type", "ind")],
                   "time" = if (frequency(y) > 1) paste(floor(outtimes),
                                                        outseason, sep = ":") else outtimes,
                   moall[,c("coefhat","tstat")])

    oind <- order(moall[,"ind"])
    moall <- moall[oind,]
    outtimes <- outtimes[oind]
    rownames(moall) <- NULL

  } else { # no outliers detected
    oeff <- NULL
    yadj <- y
    moall <- data.frame(array(dim = c(0, 4)))
    colnames(moall) <- c("type", "ind", "coefhat", "tstat")
    outtimes <- NULL
  }

  if (!is.null(logfile))
  {
    msg <- paste("\nfinal outliers\n")
    cat(msg, file = logfile, append = TRUE)
    capture.output(moall, file = logfile, append = TRUE)
  }

  structure(list(outliers = moall, y = y, yadj = yadj, cval = cval,
                 fit0 = fit.wo.outliers, # initial model fitted without outliers
                 fit = stage2$fit, effects = oeff, times = outtimes),
            class = "tsoutliers")
}

locate.outliers.oloop <- function(y, fit, types = c("AO", "LS", "TC"),
                                  cval = NULL, maxit.iloop = 4, delta = 0.7, n.start = 50, logfile = NULL)
{
  maxit <- 4 # maxit.oloop
  n <- length(y)
  s <- frequency(y)

  if (is.null(cval))
  {
    if (n <= 50) {
      cval <- 3
    } else
      if (n >= 450) {
        cval <- 4
      } else
        cval <- round(3 + 0.0025 * (n - 50), 2)
  }

  # tail(): take the last element just in case fit$call[[1]] is
  # for example "forecast::auto.arima"
  tsmethod <- ifelse(inherits(fit, "stsmFit"),
                     "stsm", tail(as.character(fit$call[[1]]), 1))
  #s <- frequency(y)
  moall <- data.frame(matrix(nrow = 0, ncol=4,
                             dimnames = list(NULL, c("type", "ind", "coefhat", "tstat"))))
  iter <- 0

  # index of initial residuals

  if (inherits(fit, "Arima")) {
    tmp <- fit$arma[6] + fit$arma[5] * fit$arma[7]
    id0resid <- if (tmp > 1) seq.int(tmp) else c(1, 2)
  } else
    if (inherits(fit, "stsmFit")) {
      id0resid <- seq_len(n - length(fit$model@diffy))
    } else
      stop("unexpected type of fitted model")

  # begin outer loop

  #while (TRUE)
  while (iter < maxit)
  {
    # extract the necessary information from the fitted model
    # parameter estimates and residuals

    pars <- switch(tsmethod,
                   "auto.arima" = , "arima" = coefs2poly(fit),
                   "stsm" = stsm::char2numeric(fit$model))

    ##NOTE by default residuals(fit, standardised = FALSE)
    # only relevant for "stsm" but the argument could set here
    # explicitly, it would be ignored if "fit" is an "Arima" object
    resid <- residuals(fit)

    if (any(abs(na.omit(resid[id0resid])) > 3.5 * sd(resid[-id0resid], na.rm = TRUE)))
    {
      ##FIXME
      # see add factor 3.5 as argument

      # this was necessary since in one series (I think it was hicp[["000000"]])
      # the first residuals were too erratic and caused the procedure to detect
      # too many outliers in the whole series

      resid[id0resid] <- 0
      # warning(paste("the first", tail(id0resid, 1), "residuals were set to zero"))

      ##NOTE
      # if this warning is returned
      # the first observations of the series may need to be inspected
      # for possible outliers, since those observations were ignored by the procedure
    }

    # locate possible outliers

    mo <- locate.outliers.iloop(resid = resid, pars = pars, cval = cval,
                                types = types, maxit = maxit.iloop, delta = delta, n.start = n.start,
                                logfile = logfile)

    # discard outliers identified at consecutive time points (if any)
    #
    # this is done by type of outlier, e.g., two or more consecutive LS are
    # replaced by one LS (the one with the highest abs(t-statistic);
    # it is still possible to get two or more consecutive outliers of different
    # type at consecutive time points;
    # alternatively consecutive outliers of any type could be replaced by a
    # single outlier. For the time being, this is kept in this way so that
    # for example, the following sequence AO1,LS2,LS3,LS4 can collapse to
    # AO1,LS4, i.e., the AO is kept
    #
    # it is not enough to do this in locate.outliers.iloop(), because
    # only the outliers detected at a given interation are checked there;
    # here, the whole set of outliers returned by locate.outliers.iloop()
    # is checked

    if (nrow(mo) > 0)
    {
      rmid <- c(
        find.consecutive.outliers(mo, "IO"),
        find.consecutive.outliers(mo, "AO"),
        find.consecutive.outliers(mo, "LS"),
        find.consecutive.outliers(mo, "TC"),
        find.consecutive.outliers(mo, "SLS"))
      #rmid <- NULL
      #for (type in types)
      #  rmid <- c(rmid, find.consecutive.outliers(mo, type))

      if (length(rmid) > 0)
        # do not use is.null(rmid), 'rmid' may be NULL or 'character(0)'
      {
        #changed in version 0.6-4
        #mo <- mo[-as.numeric(rownames(mo[rmid,])),]
        #mo <- mo[-match(as.numeric(rownames(mo[rmid,])), rownames(mo)),]
        #changed in version 0.6-5
        mo <- mo[-rmid,]
      }
    }

    # remove duplicates (if any)
    # similar to locate.outliers.iloop(), if an outlier is detected at
    # an observation where some type of outliers was already detected in
    # a previous run, the outlier that was detected first is kept

    if (nrow(mo) > 0 && iter > 0)
    {
      id.dups <- na.omit(match(moall[,"ind"], mo[,"ind"]))
      if (length(id.dups) > 0)
        mo <- mo[-id.dups,]
      # no problems with mo[-id.dups,]
      # the two dimensions of a matrix are kept in data.frame "mo" even if
      # "mo" contains one element after removing duplicates
    }

    ##FIXME
    #at this point it could be checked for consecutive outliers of a given type
    #and keep the one detected in a previous iteration
    #if (nrow(mo) > 0 && iter > 0)
    #tmp <- rbind(mo, moall)
    #rmid <- c(
    #  find.consecutive.outliers(tmp, "IO"),
    #  find.consecutive.outliers(tmp, "AO"),
    #  find.consecutive.outliers(tmp, "LS"),
    #  find.consecutive.outliers(tmp, "TC"),
    #  find.consecutive.outliers(tmp, "SLS"))
    #print(rmid)

    if (!is.null(logfile))
    {
      msg <- paste("\noloop, iteration:", iter, "\n")
      cat(msg, file = logfile, append = TRUE)
      capture.output(mo, file = logfile, append = TRUE)
    }

    # this must be here, not before the previous "if" statement
    # since it may modify "mo"

    if (nrow(mo) == 0)
      break

    moall <- rbind(moall, mo)

    # remove the effect of outliers on the data and
    # fit the model for the adjusted series

    oeff <- outliers.effects(mo = mo, n = n, weights = TRUE,
                             delta = delta, pars = pars, n.start = n.start, freq = s)

    # 'y' is overwritten; 'oeff' is based on 'mo' not 'moall'

    y <- y - rowSums(oeff)

    switch(tsmethod,
           ##FIXME
           #if 'fit' includes intercept or drift, pass here (it could be done based on names of coef(fit))

           # do not modify and evaluate the call, i.e. do not run eval(fit$call)
           # since it will run the model selection procedure (if tsmethod = "auto.arima")
           # here we only want to refit the model (not choose or select a model)

           "auto.arima" = fit <- arima(y, order = fit$arma[c(1,6,2)],
                                       seasonal = list(order = fit$arma[c(3,7,4)])),

           # this reuses arguments passed to the optimization method, e.g. method = "CSS",
           # if they were specified when the input object "fit" was created,
           # (for example through argument "args.tsmethod" of function "tso")

           "arima" = {
             fitcall <- fit$call
             fitcall$x <- y
             # rename fitcall$series since it can be a very long character string to store
             #fitcall$series <- "x"
             fit <- eval(fitcall)
           },

           "stsm" = {
             fitcall <- fit$call
             ##NOTE
             # fitcall$x contains the model, not fitcall$m, since now "stsmFit" is called
             # instead of "maxlik.td.optim" and the other functions
             fitcall$x@y <- y
             dy <- fitcall$x@fdiff(y, frequency(y))
             fitcall$x@diffy <- dy
             if (!is.null(fitcall$x@ssd))
               fitcall$x@ssd <- Mod(fft(as.numeric(dy)))^2 / (2*pi*length(dy))
             ##NOTE
             #last parameter estimates, fit$pars, could be used as starting values
             #fitcall$x@pars[] <- fit$pars
             fit <- eval(fitcall)
           }
    )

    if (!is.null(logfile))
    {
      msg <- paste("\nmodel chosen and fitted for the adjusted series:\n")
      cat(msg, file = logfile, append = TRUE)
      capture.output(fit, file = logfile, append = TRUE)
    }

    iter <- iter + 1
  } # end while

  if (iter == maxit)
    warning(paste("stopped when", sQuote("maxit.oloop = 4"), "was reached"))

  # time points with multiple types of potential outliers are not expected
  # since they are removed at each iteration
  # stopifnot(!any(duplicated(moall[,"ind"])))

  if (any(duplicated(moall[,"ind"])))
  {
    # stop for debugging
    # so far this event has not occurred
    stop("unexpected duplicates since they are handled within the loop above")
  }

  # "coefs" is not actually used but keep so far, it may used to see
  # the estimates for external regressors "xreg" which are not included in "pars"

  list(fit = list(coefs = coef(fit), pars = pars,
                  resid = resid, n = n), outliers = moall, iter = iter)
}
