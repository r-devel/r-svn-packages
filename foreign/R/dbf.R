# Copyright 2000-2001 (c) Nicholas Lewin-Koh 
dbf.read <- function(filen) {
  filen <- path.expand(filen)
  df <- .Call("Rdbfread", as.character(filen), PACKAGE="foreign")
  onames <- names(df)
  inames <- make.names(onames, unique=TRUE)
  names(df) <- inames
  if (!(identical(onames, inames))) {
    for (i in 1:length(onames))
      if (!(identical(onames[i], inames[i]))) 
        cat("Field name: ", onames[i], " changed to: ", inames[i], "\n")
  }
  df <- data.frame(lapply(df,
    function(x) {if(is.character(x)) {factor(x)} else x }))
  df
}


dbf.write <- function(dataframe, filename, factor2char=TRUE) {
# need to check precision, and that factors are converted to character
# how to handle NAs?
  if (!is.data.frame(dataframe)) stop("not a data frame")
  if (!all(complete.cases(dataframe))) stop("NAs not permitted")
  if (any(sapply(dataframe, function(x) !is.null(dim(x)))))
    stop("Can't handle multicolumn columns")
  if (factor2char) {
    dataframe <- as.data.frame(lapply(dataframe, function(x) {
      if (is.factor(x)) {
        x <- as.character(x)
        class(x) <- "AsIs"
        x
     } else x } ))
  } else {
    dataframe <- as.data.frame(lapply(dataframe, function(x) {
      if (is.factor(x)) {
        x <- as.integer(x)
        x
     } else x } ))
  }
  dataframe <- as.data.frame(lapply(dataframe, function(x) {
    if (is.logical(x)) {
      x <- as.integer(x)
      x
    } else x } ))
  m <- ncol(dataframe)
  precision <- integer(m)
  scale <- integer(m)
  dfnames <- names(dataframe)
  for (i in 1:m) {
    nlen <- nchar(dfnames[i])
    if (is.integer(dataframe[,i])) {
      rx <- range(dataframe[,i])
      mrx <- as.integer(max(ceiling(log10(abs(rx))))+3)
      precision[i] <- as.integer(max(nlen, mrx))
      if (precision[i] > 19) precision[i] <- as.integer(19)
      scale[i] <- as.integer(0)
    } else if (is.double(dataframe[,i])) {
      precision[i] <- as.integer(19)
      rx <- range(dataframe[,i])
      mrx <- as.integer(max(ceiling(log10(abs(rx)))))
      scale[i] <- as.integer(precision[i] - ifelse(mrx > 0, mrx+3, 3))
      if (scale[i] > 15) scale[i] <- as.integer(15)
    } else if (is.character(dataframe[,i])) {
      mf <- max(nchar(dataframe[,i]))
      precision[i] <- as.integer(max(nlen, mf))
      if (precision[i] > 254) precision[i] <- as.integer(254)
      scale[i] <- as.integer(0)
    } else stop("unknown column type in data frame")
  }
  invisible( .External("DoWritedbf", as.character(filename), 
    dataframe, as.integer(precision), as.integer(scale), PACKAGE="foreign"))
}
