#' @title Get old SEUX estimates
#'
#' @description Obtains the estimates of U and X from the method in Chisholm et al. (2016).
#'
#' @param S Vector of integers, number of detected extant species
#'  at each timestep.
#'
#' @param E Vector of integers, number of detected extinct species
#'  at each timestep.
#'
#' @param U_T Integer, the number of undetected extant species remaining at the end of the
#'  timeseries. Typically an assumption. Default value is 0.
#'
#' @details This function uses a simplified method for obtaining the estimates.
#'  The number of species discovered at the previous timestep
#'
#'   \eqn{ d_{t-1} = \nu_t U_{t-1} }
#'
#'  can be obtained from the data. Then the 3rd line from Eq. 1 in 
#'  Chisholm et al. (2016) can be rerranged
#'
#'   \eqn{ U_{t-1} = ( U_t + d_{t-1} ) / ( 1 - \mu_{t-1} ) }
#'
#'  Given an assumed number of undetected species remaining at the
#'  end of the timeseries, \eqn{U_T} (typically assumed \eqn{U_T = 0}),
#'  then U_t can be solved by working backwards through the timeseries.
#'  
#' @return A dataframe with three columns: U_old, X_old, mu_estimate
#'  \describe{
#'   \item{U_old:}{Floats, an estimate of the number of undetected extant species at each timestep.}
#'   \item{X_old:}{Floats, an estimate of the number of undetected extinct species at each timestep.}
#'   \item{mu_estimate:}{Floats, an estimate of the extinction rate at each timestep (calculated from detected species' extinctions. Note that, at the final timestep, an estimate cannot be made (NA value).}
#'  }
#'
#' @export

get_old_estimate <- function( S, E, U_T=0) {

    # check user inputs
    # ---

    if ( ! is.vector(S) ) { stop("Input S needs to be a vector.") }
    if ( ! is.vector(E) ) { stop("Input E needs to be a vector.") }
    if ( length(E) != length(S) ) { stop("Input vectors S and E must be have length.") }
    if ( any(S %% 1 != 0) ) { stop("Numbers of detected extant (S) needs to be whole no.") }
    if ( any(E %% 1 != 0) ) { stop("Numbers of detected extinct (E) needs to be whole no.") }


    # calculations
    # ---

    T_idx <- length(S)

    mu <- (E[-1] - E[-T_idx]) / S[-T_idx]   # estimate of extinction rate at each timestep
    U <- replicate( T_idx, 0 )              # place to store U_t estimates calculated the old way
    U[T_idx] <- U_T

    # if d == NULL
    d <- S[-1] - S[-T_idx] + E[-1] - E[-T_idx]

    for ( t in T_idx:2 ) {

        # This equation comes from rearranging the 3rd line of Eq. 1
        # in Chisholm et al. (2016), where we have used
        #  \nu_t U_{t-1} = \d_{t-1}

        U[t-1] <- ( U[t] + d[t-1] ) / ( 1 - mu[t-1] )

    }

    # calculate X_t from U_t
    N <- S[1] + E[1] + U[1]
    X <- N - E - S - U

    mu <- c(mu,NA) # there is no mu_T

    return( data.frame( U_old=U, X_old=X, mu_estimate=mu ) )

}
