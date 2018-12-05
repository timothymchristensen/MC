#
# KL distance objective function to compute identified set corresponding to true probabilities p0
#
klobj <- function(p0,pars,lower_bounds,upper_bounds){
  #
  pars <- (lower_bounds + upper_bounds*exp(pars))/(1+exp(pars))
  if(any(is.nan(pars))){pars[is.nan(pars)]<-upper_bounds[is.nan(pars)]}
  #
  ppars <- par_to_prob(pars)
  #
  ret <- sum(p0*log(p0/ppars))
  #
  return(ret)
}
#
# KL distance bjective function with some fixed parameters
#
klobj_fixed <- function(free_pars,fixed_pars,fixed_ix,p0,lower_bounds,upper_bounds){
  #
  pars <- numeric(6)
  #
  pars[-fixed_ix] <- free_pars
  pars[fixed_ix] <- fixed_pars
  #
  ret <- klobj(p0,pars,lower_bounds,upper_bounds)
  #
  return(ret)
}
#
# Minimized KL distance objective function with a subset of parameters fixed (note: all should be in transformed form)
#
kldist_fixed <- function(fixed_pars,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds){
  #
  opt <- nlm(klobj_fixed,p=initial_vals,fixed_pars=fixed_pars,fixed_ix=fixed_ix,p0=p0,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  #
  return(opt$minimum)
}
#
# Minimized KL distance objective function with a subset of parameters fixed (note: fixed should be in raw form; initial should be in transformed form)
#
kldist_fixed_raw <- function(fixed_pars,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds){
  #
  fixed_pars <- log((fixed_pars-lower_bounds[fixed_ix])/(upper_bounds[fixed_ix] - fixed_pars))
  opt <- nlm(klobj_fixed,p=initial_vals,fixed_pars=fixed_pars,fixed_ix=fixed_ix,p0=p0,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  #
  return(opt$minimum)
}
#
# Find minimum point for Delta_1 or beta_1 within tolerance
#
find_Delta_1_lower <- function(p0,fixed_ix,initial_vals,lower_bounds,upper_bounds,tol){
  #
  grid_lower <- lower_bounds[fixed_ix]
  grid_upper <- upper_bounds[fixed_ix]
  stepsize <- (grid_upper-grid_lower)/2
  #
  kl_lower <- kldist_fixed_raw(grid_lower,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
  kl_upper <- kldist_fixed_raw(grid_upper,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
  #
  while(stepsize>0.0001){
    #
    grid <- seq(grid_lower,grid_upper,by=stepsize)
    grid_short <- seq(grid_lower+stepsize,grid_upper-stepsize,by=stepsize)
    if(length(grid_short)==1){
      kl_grid <- c(kl_lower,kldist_fixed_raw(grid_short,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds),kl_upper)
    }else{
      kl_grid <- c(kl_lower,apply(matrix(grid_short),1,function(x) kldist_fixed_raw(x,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)),kl_upper)
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
# Find maximum point for Delta_1 within tolerance
#
find_Delta_1_upper <- function(p0,initial_vals,lower_bounds,upper_bounds,tol){
  #
  fixed_ix <- 1
  #
  # try upper limit of interval
  kl_upper <- kldist_fixed_raw(upper_bounds[1],fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
  #
  if(kl_upper<tol){
    ret <- upper_bounds[1]
  }else{
    #
    stepsize <- 0.0001
    upper_temp <- upper_bounds[1]
    #
    while(kl_upper>tol&&upper_temp>lower_bounds[1]+stepsize){
      #
      upper_temp <- upper_temp-stepsize
      kl_upper <- kldist_fixed_raw(upper_temp,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
    }
    #
    ret <- upper_temp
  }
  #
  return(ret)
}
#
# Find minimum point for Delta_1 or beta_1 within tolerance
#
find_beta_1_lower <- function(p0,fixed_ix,initial_vals,lower_bounds,upper_bounds,tol){
  #
  grid_lower <- lower_bounds[fixed_ix]
  grid_upper <- upper_bounds[fixed_ix]
  stepsize <- (grid_upper-grid_lower)/2
  #
  kl_lower <- kldist_fixed_raw(grid_lower,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
  kl_upper <- kldist_fixed_raw(grid_upper,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
  #
  while(stepsize>0.0001){
    #
    grid <- seq(grid_lower,grid_upper,by=stepsize)
    grid_short <- seq(grid_lower+stepsize,grid_upper-stepsize,by=stepsize)
    if(length(grid_short)==1){
      kl_grid <- c(kl_lower,kldist_fixed_raw(grid_short,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds),kl_upper)
    }else{
      kl_grid <- c(kl_lower,apply(matrix(grid_short),1,function(x) kldist_fixed_raw(x,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)),kl_upper)
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
  ret <- upper
  #
  return(ret)
}
#
# Find maximum point for Delta_1 or beta_1 within tolerance
#
find_beta_1_upper <- function(p0,fixed_ix,initial_vals,lower_bounds,upper_bounds,tol){
  #
  grid_lower <- lower_bounds[fixed_ix]
  grid_upper <- upper_bounds[fixed_ix]
  stepsize <- (grid_upper-grid_lower)/2
  #
  kl_lower <- kldist_fixed_raw(grid_lower,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
  kl_upper <- kldist_fixed_raw(grid_upper,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)
  #
  while(stepsize>0.0001){
    #
    grid <- seq(grid_lower,grid_upper,by=stepsize)
    grid_short <- seq(grid_lower+stepsize,grid_upper-stepsize,by=stepsize)
    if(length(grid_short)==1){
      kl_grid <- c(kl_lower,kldist_fixed_raw(grid_short,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds),kl_upper)
    }else{
      kl_grid <- c(kl_lower,apply(matrix(grid_short),1,function(x) kldist_fixed_raw(x,fixed_ix,p0,initial_vals,lower_bounds,upper_bounds)),kl_upper)
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
  ret <- lower
  #
  return(ret)
}
#
# END



