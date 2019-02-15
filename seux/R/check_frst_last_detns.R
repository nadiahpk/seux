#' @title Check first and last detections
#'
#' @description Performs various checks on the first and last detection dates
#'  -- e.g. first detection is always less than or equal to last detection --
#'  to ensure that they make sense.
#'
#' @param frstDetn Vector of integers, species' year of first detection.
#'
#' @param lastDetn Vector of integers, species' year of last detection.
#'
#' @return A list with elements: flag, message.
#'  \describe{
#'  \item{flag:}{Boolean, is TRUE if the data passed the checks.}
#'  \item{message:}{String, a message explaining why the data didn't pass the checks (can be
#'   used e.g. in the stop function).}
#'  }
#'
#' @export

check_frst_last_detns <- function( frstDetn, lastDetn ) {

    # check the dates are integers

    flag <- TRUE
    message <- ""

    if ( any( sapply(frstDetn, function(v) v %% 1 != 0 ) ) ) {
        flag <- FALSE
        message <- "First-detection dates must be whole numbers."
    }

    if ( flag ) {
        if ( any( sapply(lastDetn, function(v) v %% 1 != 0 ) ) ) {
            flag <- FALSE
            message <- "Last-detection dates must be whole numbers."
        }
    }

    # check no first detection dates are after last detection dates

    if ( flag ) {
        if ( any(frstDetn > lastDetn) ) {
            flag <- FALSE
            message <- "First-detection dates must not be after last-detection dates."
        }
    }

    return( list( flag=flag, message=message ) )
}


