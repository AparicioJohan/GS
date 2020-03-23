
#' Check phenotypic Information
#'
#' \code{checkPhen} returns phenotypic information cleaned
#'
#' @param data A data.frame
#'
#' @return data depurated
#' @export
#'
#' @examples
#' # library(sommer)
#' # data(DT_cpdata)
#' # phen <- DT_cpdata
#'
#' # checkPhen(phen)
#'
#' @author Johan Aparicio, \email{j.aparicio@cgiar.org}
checkPhen <- function(data){

  datafull <- type.convert(data)
  name <- names(datafull)[1]

  dataremo <- datafull[,-1]
  dTr <- sapply(dataremo, class)
  dSd <- na.omit(suppressWarnings(apply(dataremo, 2, sd, na.rm=T)))

  dTrf <- paste(names(dTr[dTr=="factor"]), collapse = ", " )
  dTrc <- paste(names(dTr[dTr=="character"]), collapse = ", " )
  dSdr <- paste(names(dSd[dSd==0]), collapse = ", " )

  if(any(dTr=="factor") ) warning(paste("factors were removed from phen: ", dTrf),call. = FALSE)
  if(any(dTr=="character") ) warning(paste("characters were removed from phen: ", dTrc ),call. = FALSE)
  if(any(dSd==0)) warning(paste("traits with zero variance were removed from phen: ",dSdr ),call. = FALSE)

  dataremo <- dplyr::select_if(dataremo, is.numeric)

  sdnull <- function(x) sd(x,na.rm = T)!=0
  dataremo <- dplyr::select_if(dataremo, sdnull )

  datafull <- cbind(datafull[,1],dataremo)
  names(datafull)[1] <- name

  return(datafull)

}
