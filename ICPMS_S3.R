new_signal.sd.matrixes <- function(signal = matrix(), signalSd = matrix()) {
  stopifnot(is.matrix(signal))
  stopifnot(is.matrix(signalSd))
  stopifnot(identical(dim(signal), dim(signalSd)))
  
  structure(list(signal = signal, signalSd = signalSd),
            class = "signal.sd.matrixes"
  )
}

divideSignalSdMatrixes <- function(x, y) {
  UseMethod("divideSignalSdMatrixes", x)
}

subtractSignalSdMatrixes <- function(x, y) {
  UseMethod("subtractSignalSdMatrixes", x)
}

divideSignalSdMatrixes.signal.sd.matrixes <- function(x, y) {
  signal <- x$signal / y$signal
  signalSd <- sqrt((x$signalSd/x$signal)^2 + (y$signalSd/y$signal)^2) * signal
  return(new_signal.sd.matrixes(signal, signalSd))
}

subtractSignalSdMatrixes.signal.sd.matrixes <- function(x, y) {
  signal <- x$signal - y$signal
  signalSd <- sqrt(x$signalSd^2 + y$signalSd^2)
  return(new_signal.sd.matrixes(signal, signalSd))
}

'[.signal.sd.matrixes' = function(x, i, j){
  signal <- x$signal[i, j, drop = FALSE]
  signalSd <- x$signalSd[i, j, drop = FALSE]
  return(new_signal.sd.matrixes(signal, signalSd))
} 

dim.signal.sd.matrixes = function(x) {
  return(dim(x$signal))
}
