#' @title Get SEUX model inputs
#'
#' @description Turns first and last detection dates into input for the SEUX model.
#'
#' @param frstDetn Vector of integers, species' year of first detection.
#'
#' @param lastDetn Vector of integers, species' year of last detection.
#'
#' @param collapse_timesteps Boolean. If set to FALSE, each year is treated as a timestep
#'  in the model. If set to TRUE, timesteps are collapsed together so that each timestep
#'  has at least one detected extinction. Defaults to TRUE.
#'
#' @param y0 Integer, the start year for the timeseries.
#'  Default is the year of the first detected extinction. 
#'
#' @param yf Integer, the final year for the timeseries.
#'  Default is the year of the last new discovery or the year of the last extinction,
#'  whichever is later.
#'
#' @details This function is mainly to provide a convenient way to impose constraints on the
#'  timeseries extent and collapse timesteps.
#'
#' @return A dataframe with four columns: year, S, E, d
#'  \describe{
#'   \item{year:}{Integers, the year to which each timestep corresponds.}
#'   \item{S:}{Integers, number of detected extant species at each timestep.}
#'   \item{E:}{Integers, number of detected extinct species at each timestep.}
#'   \item{d:}{Integers, number of species discovered at each timestep.}
#'  }
#'
#' @export

get_model_inputs <- function( frstDetn, lastDetn, collapse_timesteps=TRUE, y0=NULL, yf=NULL) {

    # check user input and assign values to y0 and yf if needed
    # ---

    input_okay <- check_frst_last_detns(frstDetn, lastDetn)
    if (! input_okay$flag ){
        stop( input_okay$message )
    }

    if ( is.null(y0) ) { 

        # start the timeseries with the first extinction
        y0 <- min(lastDetn) 

    }

    if ( is.null(yf) ) { 

        # end the timeseries with the last new discovery of a species
        year_last_new_discovery <- max(frstDetn) 
        year_last_extn <- max(Filter(function(x) x < max(lastDetn), lastDetn)) + 1
        yf <- max(c( year_last_new_discovery, year_last_extn ))

    }

    if ( y0 > yf ) { stop("Start year in timeseries must be before final year: ", as.character(y0), " ", as.character(yf) ) }

    # calculations
    # ---

    # get years vector and S and E vectors
    if ( collapse_timesteps ) {

        # years between y0 and yf, with at least one last-detection, plus y0 and yf
        years <- sort(unique(Filter( function(y) y <= yf & y >= y0, append( lastDetn, c(y0, yf) ) )))

        S <- sapply( years, function(y) sum(frstDetn <= y & lastDetn >= y) )  # detected extant species
        E <- sapply( years, function(y) sum(lastDetn < y) )                  # detected extinct species

    } else {

        years <- y0:yf

        S <- sapply( years, function(y) sum(frstDetn <= y & lastDetn >= y) )  # detected extant species
        E <- sapply( years, function(y) sum(lastDetn < y) )                  # detected extinct species

    }

    # discoveries at each timestep
    T_idx <- length(S)
    d <- S[-1] - S[-T_idx] + E[-1] - E[-T_idx] # in Python: d = S[1:] - S[:-1] + E[1:] - E[:-1]
    d <- c(d,NA) # there is no d_T

    return( data.frame( year=years, S=S, E=E, d=d ) )

}
