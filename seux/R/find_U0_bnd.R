#' @title Find U0 bound
#'
#' @description Finds the U0 bound, which is the greatest value of U0 such that
#'  the mid-P function returns midP_fnc(U0_bnd, ...) > alpha
#'
#' @param midP_fnc Closure, a function to be used as the mid-P function.
#'
#' @param alpha Float [0,1], analogous to a confidence level.
#'
#' @param S0 Integer, number of detected extant species
#'  at the previous timestep.
#'
#' @param S1 Integer, number of detected extant species
#'  at the current timestep.
#'
#' @param U1 Integer, number of undetected extant species
#'  at the current timestep.
#'
#' @param d0 Integer, number of species detected during
#'  the timestep.
#'
#' @param impossibleFlag Boolean, indicates whether or not the previous timestep
#'  of the sampling algorithm had sampled an impossible value (i.e. U1) for the bound.
#'  Defaults to FALSE.
#'
#' @details This function is used in a random sampling scheme (i.e. with random alpha)
#'  to obtain a sample of U bounds at a particular timestep. The confidence intervals
#'  can then be approximated by taking percentiles of the sample.
#'
#'  When this function is used in combination with a mid-P function, the random sampling
#'  scheme is equivalent to random sampling of fuzzy confidence bounds. Because the data
#'  is discrete, it can return a bound that is physically impossible i.e. U0_bnd < U1 + d0.
#'  This is necessary in order for the actual coverage to match the nominal.
#'
#'  If the bound is impossible, this will be indicated with the impossibleFlag returned.
#'  If an impossible bound should not be sampled (e.g. because the U1 is already in the
#'  impossible region), then the impossible bound can be prevented by setting impossibleFlag = TRUE
#'  in the input.
#'
#' @return A list with elements: U0_bnd, impossibleFlag.
#'  \describe{
#'  \item{U0_bnd:}{Integer, the U bound at alpha.}
#'  \item{impossibleFlag:}{Boolean, is TRUE if the U bound is impossible.}
#'  }
#'
#' @examples
#'  S0 <- 100; S1 <- 96; d0 <- 1    # from the data
#'  U1 <- 200;                      # from an earlier timestep in the sampling algorithm
#'  alpha <- 0.5;                   # a randomly-sampled confidence level
#'  midP_fnc <- central_hyper_midP
#'  res <- find_U0_bnd(midP_fnc, alpha, S0, S1, U1, d0)
#'  U0_bnd <- res$U0_bnd            # the U bound at alpha=0.5
#'
#' @export

find_U0_bnd <- function(midP_fnc, alpha, S0, S1, U1, d0, impossibleFlag = FALSE){ 


    # check what the user gave us
    # ---

    if ( typeof(midP_fnc) != "closure" ) { stop("midP_fnc needs to be a function (try midP_fnc <- central_hyper_midP).") }

    if ( ! is.numeric(alpha) ){
        stop("alpha needs to be a number between zero and one.")
    } else {
        if ( alpha < 0 | alpha > 1 ) { stop("alpha needs to be between zero and one: ", as.character(alpha)) }
    }

    if ( ! is.numeric(S0) ){
        stop("Number of detected extant in previous timestep (S0) needs to be a whole number.")
    } else {
        if ( S0 %% 1 != 0 ) { 
            stop("Number of detected extant in previous timestep (S0) needs to be whole no: ", as.character(S0))
        }
    }

    if ( ! is.numeric(S1) ){
        stop("Number of detected extant in current timestep (S1) needs to be a whole number.")
    } else {
        if ( S1 %% 1 != 0 ) { 
            stop("Number of detected extant in current timestep (S1) needs to be whole no: ", as.character(S1))
        }
    }

    if ( ! is.numeric(U1) ){
        stop("Number of undetected extant in previous timestep (U1) needs to be a whole number.")
    } else {
        if ( U1 %% 1 != 0 ) { 
            stop("Number of undetected extant in previous timestep (U1) needs to be whole no: ", as.character(U1))
        }
    }

    if ( ! is.numeric(d0) ){
        stop("Number of species discovered in previous timestep (d) needs to be a whole number.")
    } else {
        if ( d0 %% 1 != 0 ) { 
            stop("Number of species discovered in previous timestep (d) needs to be whole no: ", as.character(d0))
        }
    }


    # calculations
    # ---

    # first, check if we're in the situation where we wouldn't accept the minimum possible value of U0

    min_poss_U0 <- U1 + d0 # the minimum possible value of U0 in reality

    if ( midP_fnc(min_poss_U0, S0, S1, U1, d0) < alpha ) {

        if ( impossibleFlag ){

            # don't take two steps in to the impossible region
            U0_bnd <- min_poss_U0

        } else {

            # U0 is set to the 'impossible' value, one less than minimum value, end search
            U0_bnd <- min_poss_U0-1
            impossibleFlag <- TRUE

        }

    } else {

        # The search for U0_bnd has two phases:
        #   1. Find an interval [U0_lo, U0_hi] within which U0_bnd lies.
        #   2. A binary search between our upper and lower search bounds, U0_lo and U0_hi, to find U0_bnd.

        # 1. find an interval [U0_lo, U0_hi] within which U0_bnd lies
        # ----

        # To find the interval, we start at the minimum possible value for U0_bnd,
        # and work our way upwards until we find a U0_hi such that midP_fnc(U0_hi, ...) < alpha.
        # The step-size from U0_lo to the next U0 value is doubled each time, to account
        # for the fact that there is no upper limit to the U0 value, and that a randomly-chosen
        # alpha may be quite small.

        # initialise for the loop
        U0_lo <- min_poss_U0 # our first lower bound is the minimum possible value of U0
        step_size <- 1       # this will double at each step
        U0_hi_found <- FALSE

        while ( ! U0_hi_found ) {

            U0_hi <- U0_lo + step_size  # prospective upper search bound

            if ( midP_fnc(U0_hi, S0, S1, U1, d0) < alpha ) {

                U0_hi_found <- TRUE     # we have found an upper search bound

            } else {

                U0_lo <- U0_hi               # we didn't find an upper search bound this time, so this is new lower search bound 
                step_size <- step_size * 2   # double the step size for next time

            }
        }


        # 2. binary search between our upper and lower search bounds, U0_lo and U0_hi
        # ----

        # Once an interval [U0_lo, U0_hi] has been identified, we use a binary search
        # to find the value of U0_lo such that
        #   midP_fnc(U0_lo, ...) > alpha  and  midP_fnc(U0_lo+1, ...) < alpha.
        # This value of U0_lo is U0_bnd.
        # The binary search works by iteratively taking the mid-point U0_mid between U0_lo and U0_hi,
        # checking if it is a new low or high bound on the interval, and checking the condition above
        # to terminate the search.

        # have we already found bound?
        U0_bnd_found <- U0_hi - U0_lo == 1 # if they're one apart, then we've found the bound

        while ( ! U0_bnd_found ) {

            U0_mid <- floor( (U0_lo+U0_hi) / 2 )             # the midpoint between U0_lo and U0_hi
            alpha_mid <- midP_fnc(U0_mid, S0, S1, U1, d0)    # the alpha value at this mid-point

            # Now we must figure out if U0_mid is: (1) the actual bound, (2) a new lower search bound 
            # for the search interval, or (3) a new upper search bound.

            if ( alpha_mid == alpha ) { 
                
                # it's the actual bound (unlikely to happen)
                U0_bnd_found <- TRUE
                U0_lo <- U0_mid # the bound is stored in U0_lo

            } else {

                if ( alpha_mid < alpha ) {

                    # the mid-point is a new upper bound
                    U0_hi <- U0_mid

                } else {
                    
                    # the mid-point is a new lower bound
                    U0_lo <- U0_mid

                }

                # if the interval bounds are one apart, then we've found U0_bnd
                U0_bnd_found <- U0_hi - U0_lo == 1

            }

        }

        # When we exit the while loop above, then the bound has been found.
        # It is stored as the lower bound of the search interval, U0_lo, i.e. 
        # the greatest value of U0 for which midP_fnc >= alpha.

        U0_bnd <- U0_lo

        # Check if we've moved out of the impossible region
        if (U0_bnd > min_poss_U0 ) { impossibleFlag <- FALSE }
        # note that, if we went in with impossibleFlag == True, and U0_bnd == min_poss_U0, then we remain impossibleFlag == True

    }

    # It is sometimes possible to choose a negative number of undetected species
    # (as an impossible bound) near the end of the timeseries. Disallow this,
    # as we cannot evaluate the hypergeometric model at the next timestep with a 
    # negative value.

    if ( U0_bnd < 0 ) { U0_bnd <- 0 }

    return( list( U0_bnd=U0_bnd, impossibleFlag=impossibleFlag ) )

}


