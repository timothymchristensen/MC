#
# Transform parameters to have full support via logit transform
#
par_trans <- function(pars,lower_bounds,upper_bounds){
    #
    trans_pars <- log((pars-lower_bounds)/(upper_bounds - pars))
    #
    return(trans_pars)
    #
}
#
par_inv_trans <- function(pars_trans,lower_bounds,upper_bounds){
    #
    pars <- (lower_bounds + upper_bounds*exp(pars_trans))/(1+exp(pars_trans))
    #
    return(pars)
    #
}