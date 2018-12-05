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
    if(any(is.nan(pars))){pars[is.nan(pars)]<-upper_bounds[is.nan(pars)]}
    if(any(is.infinite(pars))){pars[is.infinite(pars)]<-upper_bounds[is.infinite(pars)]}
    #
    return(pars)
    #
}