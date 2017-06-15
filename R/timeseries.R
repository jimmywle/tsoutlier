#' @import zoo , forecast

timeseries<-function (data = NA, start = 1, end = numeric(), frequency = 1,
                      deltat = 1, class = if (nseries >
                                              1) c("mts", "ts", "matrix") else "ts", names = if (!is.null(dimnames(data))) colnames(data) else paste("Series",
                                                                                                                                                     seq(nseries)))
{
  if (is.data.frame(data))
    data <- data.matrix(data)
  if (is.matrix(data)) {
    nseries <- ncol(data)
    ndata <- nrow(data)
    dimnames(data) <- list(NULL, names)
  }
  else {
    nseries <- 1
    ndata <- length(data)
  }
  if (ndata == 0)
    stop("'ts' object must have one or more observations")
  if (missing(frequency))
    frequency <- 1/deltat
  else if (missing(deltat))
    deltat <- 1/frequency
  if (frequency > 1 && abs(frequency - round(frequency)) <
      ts.eps)
    frequency <- round(frequency)
  if (length(start) > 1L) {
    start <- start[1L] + (start[2L] - 1)/frequency
  }
  if (length(end) > 1L) {
    end <- end[1L] + (end[2L] - 1)/frequency
  }
  if (missing(end))
    end <- start + (ndata - 1)/frequency
  else if (missing(start))
    start <- end - (ndata - 1)/frequency
  if (start > end)
    stop("'start' cannot be after 'end'")
  nobs <- floor((end - start) * frequency + 1.01)
  if (nobs != ndata)
    data <- if (NCOL(data) == 1) {
      if (ndata < nobs)
        rep_len(data, nobs)
      else if (ndata > nobs)
        data[1L:nobs]
    }
  else {
    if (ndata < nobs)
      data[rep_len(1L:ndata, nobs), ]
    else if (ndata > nobs)
      data[1L:nobs, ]
  }
  attr(data, "tsp") <- c(start, end, frequency)
  if (!is.null(class) && class != "none")
    attr(data, "class") <- class
  data
}


#timeseries conversion
timeseries_conversion<- function (x = NULL, order.by = index(x),
                                  calendar = getOption("zoo.calendar", TRUE))
{




  ## process index "order.by"
  if(length(unique(MATCH(order.by, order.by))) < length(order.by))
    warning(paste("some methods for", dQuote("zoo"),
                  "objects do not work if the index entries in", sQuote("order.by"), "are not unique"))
  index <- ORDER(order.by)
  order.by <- order.by[index]

  if(is.matrix(x) || is.data.frame(x)) x <- as.matrix(x)
  if(is.matrix(x) && sum(dim(x)) < 1L) x <- NULL
  if(missing(x) || is.null(x))
    x <- numeric()
  else if(is.factor(x))
    x <- factor(rep(as.character(x), length.out = length(index))[index],
                levels = levels(x), ordered = is.ordered(x))
  else if(is.matrix(x) || is.data.frame(x))
    x <- (x[rep(1:NROW(x), length.out = length(index)), ,
            drop = FALSE])[index, , drop = FALSE]
  else if(is.atomic(x))
    x <- rep(x, length.out = length(index))[index]
  else stop(paste(dQuote("x"), ": attempt to define invalid zoo object"))


  attr(x, "oclass") <- attr(x, "class")
  attr(x, "index") <- order.by
  attr(x, "frequency") <- frequency
  class(x) <- if(is.null(frequency)) "zoo" else c("zooreg", "zoo")
  return(x)
}
