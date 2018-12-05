#
# KL distance objective function to compute identified set corresponding to true probabilities p0
#
klobj <- function(p0,pars,lower_bounds,upper_bounds,ssame){
  #
  pars <- (lower_bounds + upper_bounds*exp(pars))/(1+exp(pars))
  if(any(is.nan(pars))){pars[is.nan(pars)]<-upper_bounds[is.nan(pars)]}
  if(any(is.infinite(pars))){pars[is.infinite(pars)]<-upper_bounds[is.infinite(pars)]}
  #
  ppars <- par_to_prob(pars,ssame)
  #
  ret <- sum(p0[ppars>0]*log(p0[ppars>0]/ppars[ppars>0]))
  #
  return(ret)
}
#
# KL distance bjective function with some fixed parameters
#
klobj_fixed <- function(free_pars,fixed_pars,fixed_ix,p0,lower_bounds,upper_bounds,ssame){
  #
  pars <- numeric(length(lower_bounds))
  #
  pars[-fixed_ix] <- free_pars
  pars[fixed_ix] <- fixed_pars
  #
  ret <- klobj(p0,pars,lower_bounds,upper_bounds,ssame)
  #
  return(ret)
}
#
# Minimized KL distance objective function with a subset of parameters fixed (note: all should be in transformed form)
#
kldist_fixed <- function(fixed_pars,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame){
  #
  opt <- nlm(klobj_fixed,p=initial_vals,fixed_pars=fixed_pars,fixed_ix=fixed_ix,p0=p0,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame=ssame)
  #
  return(opt$minimum)
}
#
# Minimized KL distance objective function with a subset of parameters fixed (note: fixed should be in raw form; initial should be in transformed form)
#
kldist_fixed_raw <- function(fixed_pars,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame){
  #
  fixed_pars <- log((fixed_pars-lower_bounds[fixed_ix])/(upper_bounds[fixed_ix] - fixed_pars))
  opt <- nlm(klobj_fixed,p=initial_vals,fixed_pars=fixed_pars,fixed_ix=fixed_ix,p0=p0,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame)
  opt <- nlm(klobj_fixed,p=opt$estimate,fixed_pars=fixed_pars,fixed_ix=fixed_ix,p0=p0,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame)
  #
  return(opt$minimum)
}
#
# Find minimum point for Delta_1 or beta_1 within tolerance
#
find_int_lower <- function(p0,fixed_par,fixed_ix,initial_vals,lower_bounds,upper_bounds,ssame,tol=1e-5){
  #
  grid_lower <- lower <- lower_bounds[fixed_ix]
  grid_upper <- upper <- fixed_par
  stepsize <- (grid_upper-grid_lower)/2
  #
  kl_lower <- kldist_fixed_raw(grid_lower,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame)
  kl_upper <- kldist_fixed_raw(grid_upper,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame)
  #
  while(stepsize>0.0005){
    #
    grid <- seq(grid_lower,grid_upper,by=stepsize)
    grid_short <- grid[-c(1,length(grid))]
    if(length(grid_short)==1){
      kl_grid <- c(kl_lower,kldist_fixed_raw(grid_short,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame),kl_upper)
    }else{
      kl_grid <- c(kl_lower,apply(matrix(grid_short),1,function(x) kldist_fixed_raw(x,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame)),kl_upper)
    }#
    if(any(kl_grid<tol)){
      #
      upper <- min(grid[kl_grid<tol])
      #
      if(upper==grid_lower){
        lower <- upper
      }else{
        lower <- max(grid[grid<upper])
      }
      #
      grid_lower <- lower
      grid_upper <- upper
      stepsize <- (grid_upper-grid_lower)/2
      #
      kl_lower <- kl_grid[grid==lower]
      kl_upper <- kl_grid[grid==upper]
      #
    }else{
      # what happens if all above tolerance
      stepsize <- stepsize/2
      lower <- upper <- grid_lower
      #
      kl_lower <- kl_grid[grid==grid_lower]
      kl_upper <- kl_grid[grid==grid_upper]
    }
  }
  #
  ret <- (upper+lower)/2
  #
  return(ret)
}
#
# Find maximum point for Delta_1 or beta_1 within tolerance
#
find_int_upper <- function(p0,fixed_par,fixed_ix,initial_vals,lower_bounds,upper_bounds,ssame,tol=1e-5){
  #
  grid_lower <- lower <- fixed_par
  grid_upper <- upper <- upper_bounds[fixed_ix]
  stepsize <- (grid_upper-grid_lower)/2
  #
  kl_lower <- kldist_fixed_raw(grid_lower,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame)
  kl_upper <- kldist_fixed_raw(grid_upper,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame)
  #
  while(stepsize>0.0005){
    #
    grid <- seq(grid_lower,grid_upper,by=stepsize)
    grid_short <- grid[-c(1,length(grid))]
    if(length(grid_short)==1){
      kl_grid <- c(kl_lower,kldist_fixed_raw(grid_short,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame),kl_upper)
    }else{
      kl_grid <- c(kl_lower,apply(matrix(grid_short),1,function(x) kldist_fixed_raw(x,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds,ssame)),kl_upper)
    }
    #
    if(any(kl_grid<tol)){
      #
      lower <- max(grid[kl_grid<tol])
      #
      if(lower==grid_upper){
        upper <- lower
      }else{
        upper <- min(grid[grid>lower])
      }
      #
      grid_lower <- lower
      grid_upper <- upper
      stepsize <- (grid_upper-grid_lower)/2
      #
      kl_lower <- kl_grid[grid==lower]
      kl_upper <- kl_grid[grid==upper]
      #
    }else{
      # what happens if all above tolerance
      stepsize <- stepsize/2
      lower <- upper <- grid_upper
      #
      kl_lower <- kl_grid[grid==grid_lower]
      kl_upper <- kl_grid[grid==grid_upper]
    }
  }
  #
  ret <- (upper+lower)/2
  #
  return(ret)
}
#
# END



