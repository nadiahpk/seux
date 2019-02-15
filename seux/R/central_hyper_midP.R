#' @title Central hypergeometric distribution mid-P function
#'
#' @description Calculates the mid-P function's value for the 
#'  central hypergeometric SEUX model. 
#'
#' @param U0 Integer, number of undetected extant species
#'  at the previous timestep.
#'
#' @param S0 Integer, number of detected extant species
#'  at the previous timestep.
#'
#' @param S1 Integer, number of detected extant species
#'  at the current timestep.
#'
#' @param d0 Integer, number of species detected during
#'  the timestep.
#'
#' @param U1 Integer, number of undetected extant species
#'  at the current timestep.
#'
#' @details The mid-P function is 
#'
#'  \eqn{ P_M(U_b \geq U^*|\phi) = 0.5 [ P(\Phi \leq \phi | U^*) + 1 - P(\Phi \geq \phi | U^*) ] }.
#'
#'  Sampling according the mid-P function is equivalent to sampling fuzzy confidence bounds.
#'
#' @return alpha: Float [0,1], analogous to a confidence level.
#'
#' @examples
#'  S0 <- 100; S1 <- 96; d0 <- 1    # from the data
#'  U1 <- 200;                      # from an earlier timestep in the algorithm
#'  central_hyper_midP(211, S0, S1, U1, d0) # 0.5491703
#'  central_hyper_midP(212, S0, S1, U1, d0) # 0.4838616
#'
#' @importFrom stats phyper
#'
#' @export


central_hyper_midP <- function(U0, S0, S1, U1, d0){

    # Calculate the mid-P function's value for the central hypergeometric model. 
    # The mid-P function is:
    #   P_M(U_b \geq U^*|\phi) = 0.5 [ P(\Phi \leq \phi | U^*) + 1 - P(\Phi \geq \phi | U^*) ].

    alpha <- 0.5 * ( phyper( U1+d0, U0, S0, U1+S1 ) + phyper( U1+d0-1, U0, S0, U1+S1 ) )

    return(alpha)

}

