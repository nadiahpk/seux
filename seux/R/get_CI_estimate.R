#' @title Get confidence intervals and estimates
#'
#' @description Obtains confidence intervals and for U and X using sampling
#'  of the mid-P function for the model (e.g. the central hypergeometric model)
#'
#' @param S Vector of integers, the number of detected extant species at each timestep.
#'
#' @param E Vector of integers, the number of detected extinct species at each timestep.
#'
#' @param nreps Integer, the number of random samples of U bounds to take to calculate the CIs.
#'  Default is 10000.
#'
#' @param percentile Float [1,100], the percentile of the random sample of U bounds that
#'  is used to calculate the CI.
#'  Default is 95 (for a 95\% CI).
#'
#' @param U_T Integer, the number of undetected extant species remaining at the end of the
#'  timeseries. Typically an assumption. Default is 0.
#'
#' @param midP_fnc Closure, a function to be used as the mid-P function.
#'  Default will use central_hyper_midP, which is the mid-P function for the
#'  central hypergeometric SEUX model.
#'
#' @details This function uses a random sampling scheme to obtain a sample of U bound values
#'  at each timestep, through repeated calls to the function find_U0_bnd.
#'  The confidence intervals then approximated by taking percentiles of the sample.
#'  
#'  When this function is used in combination with a mid-P function, the random sampling
#'  scheme is equivalent to random sampling of fuzzy confidence bounds. Because the data
#'  is discrete, the coverage of the confidence intervals will not exactly match the nominal
#'  level, particularly when U is low (at the end of the timeseries). 
#'
#' return( data.frame( X_lo=X_lo, X_hi=X_hi, X_mean=X_mean, U_lo=U_lo, U_hi=U_hi, U_mean=U_mean ) )
#' @return A data frame with 6 columns.
#'  \describe{
#'   \item{X_lo:}{Integers, the lower confidence bounds at each timestep on the number of undetected extinct species.}
#'   \item{X_hi:}{Integers, the upper confidence bounds at each timestep on the number of undetected extinct species.}
#'   \item{X_mean:}{Float, an estimate at each timestep of the number of undetected extinct species.}
#'   \item{U_lo:}{Integers, the lower confidence bounds at each timestep on the number of undetected extant species.}
#'   \item{U_hi:}{Integers, the upper confidence bounds at each timestep on the number of undetected extant species.}
#'   \item{U_mean:}{Float, an estimate at each timestep of the number of undetected extant species.}
#'  }
#'
#' @importFrom stats runif
#' @importFrom stats quantile
#'
#' @export


get_CI_estimate <- function( S, E, nreps=10000, percentile=95, U_T=0, midP_fnc=NULL ) {

    # check user inputs
    # ---

    if ( ! is.vector(S) ) { stop("Input S needs to be a vector.") }
    if ( ! is.vector(E) ) { stop("Input E needs to be a vector.") }

    if ( ! is.numeric(U_T) ) {
        stop("Input U_T needs to be a whole number.")
    } else {
        if ( U_T %% 1 != 0 ) { stop("Input U_T needs to be a whole number: ", as.character(U_T) ) }
    }

    if ( is.null(midP_fnc) ) {
        midP_fnc <- central_hyper_midP
    } else {
        if ( typeof(midP_fnc) != "closure" ) { stop("midP_fnc needs to a function.") }
    }

    if ( nreps %% 1 != 0 ) { stop("Input nreps needs to be whole no: ", as.character(nreps)) }
    if ( ! is.numeric(percentile) ) { stop("Percentile needs to be a percentage number: ", as.character(percentile)) }
    if ( percentile < 1 ) { stop("Percentile needs to be a percentage (e.g. 95): ", as.character(percentile)) }

    if ( length(E) != length(S) ) { stop("Input vectors S and E must be have length.") }

    # calculations
    # ---

    T_idx <- length(S)
    d <- S[-1] - S[-T_idx] + E[-1] - E[-T_idx] # no. spp discovered each timestep

    # there is an absolute lower bound on U, which is sum_{i=t}^{T} d
    U_min_possible <- c( rev(cumsum(rev( d[1:length(d)] ))), 0 )

    # a place to store our replicates
    UM <- matrix( NA, nreps, T_idx )
    XM <- matrix( NA, nreps, T_idx )

    for ( nrep in 1:nreps ) {

        # Initialise an empty vector of undetected extant spp at each timestep,
        # with our assumed value of U_T at the end
        U <- replicate( T_idx, 0 ) # [0]*T_idx
        U[T_idx] <- U_T
        impossibleFlag <- FALSE

        # work backwards in time
        for ( t in T_idx:2 ) {

            alpha <- runif(1)
            S0 <- S[t-1]
            S1 <- S[t]
            d0 <- d[t-1]
            U1 <- U[t]

            res <- find_U0_bnd(midP_fnc, alpha, S0, S1, U1, d0, impossibleFlag)
            U[t-1] <- res$U0
            impossibleFlag <- res$impossibleFlag

        }

        # In order to give intervals that don't cross into negative values
        # swap U values for the minimum possible U where necessary.
        # This will still give the correct coverage because the true U values
        # must be greater than or equal to U_min_possible
        U <- pmax( U, U_min_possible )


        # calculate X_t from U_t
        N <- S[1] + E[1] + U[1] # assumes X(time=start) = 0
        X <- N - E - S - U

        # an absolute lower bound on X is 0
        X <- pmax(X,0)

        # store results
        UM[nrep,] <- U
        XM[nrep,] <- X

    }

    # calculate statistics on sample

    q_lo <- ( 1 - percentile/100 ) / 2
    q_hi <- 1 - ( 1 - percentile/100 ) / 2

    X_lo   <- apply( XM, 2, function(x) quantile(x,q_lo) )
    X_hi   <- apply( XM, 2, function(x) quantile(x,q_hi) )
    X_mean <- apply( XM, 2, function(x) mean(x) )

    U_lo   <- apply( UM, 2, function(x) quantile(x,q_lo) )
    U_hi   <- apply( UM, 2, function(x) quantile(x,q_hi) )
    U_mean <- apply( UM, 2, function(x) mean(x) )

    return( data.frame( X_lo=X_lo, X_hi=X_hi, X_mean=X_mean, U_lo=U_lo, U_hi=U_hi, U_mean=U_mean ) )
}
