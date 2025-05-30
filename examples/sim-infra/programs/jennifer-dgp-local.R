# Generate data from the DGP.
genData <- function(
    world='medicare',
    nT = 1000,
    nC = 2000,
    sigy_multiplier=1,
    sigu_multiplier=1,
    sigv_multiplier=1,
    sigt_multiplier=1,
    nonlin_frac = NULL,
    r2_w_tau=0,
    r2_w_mu=0,
    r2_w_pi=0,
    rho=0,
    r2_error=NULL,
    weights=TRUE,
    trt_eff_scenario='het',
    seed=NULL,
    uv_dist='normal',
    uv_dist_df=NULL, targeted_sel=TRUE, reparam=TRUE, center_ate=FALSE, fast_r2a=FALSE, driver_sheet='sim_driver_rereparam'){
  #-------------------------------------------------------------------------------------------------
  # INPUTS:
  # world   "medicare" or "edu".  Defines baseline configuration values. (ICC, resid r-squared, etc.)
  # nT:     Number of treated units.
  # nC:     Number of control units.
  # r2trt_multiplier: Multiplier for r2trt = var(taux)/(var(taux) + var(mux)) ratio from driver file.
  # r2uv_multiplier:  Multiplier for r2uv = sig2v/(sig2v + sig2u) ratio from driver file.
  # rho:    Correlation coefficient for ui and vi.
  # weights:  Indicator for using weighted data.
  #           If FALSE, weights are set to 1.
  #           If TRUE, weights are sampled from .RDS files with practice or school-size quantiles.
  # trt_eff_scenario: "het" or "homog" to determine treatment effect function.
  # seed:  Optional seed for replicating datasets.
  # uv_dist: "normal" or "t".  Sets bivariate distribution for (ui,vi).
  # uv_dist_df: Degrees of freedom if uv_dist=t.  Ignored if uv_dist=normal.
  #
  # OUTPUTS:
  # One dataframe generated from the data-generating process as specified by driver file,
  # practice or school size quantile file, and user inputs.
  
  #- Error checking --------------------------------------------------------------------------------
  
  # Set seed if specified.
  if(!is.null(seed)) set.seed(seed)
  
  # Error checking for treatment effect scenario input.
  if(!trt_eff_scenario %in% c('homog','het')){
    stop('trt_eff_scenario must be het or homog.')
  }
  
  # Error checking for inputting an r2trt_multiplier with homogeneous treatment effects.
  # (r2trt = 0 by definition for homog treatment effects; no variability in x-related trt effect.)
  if(trt_eff_scenario=='homog' & sigt_multiplier!=1){
    warning('sigt_multiplier not applied for homogeneous treatment effects.')
  }
  
  # Error checking for random effects distribution.
  if(!uv_dist %in% c('normal','t')){
    stop("uv_dist must be set to 'normal' or 't'.")
  }
  
  # Warning for ignoring degrees of freedom if uv ~ MVN.
  if(uv_dist=='normal' & !is.null(uv_dist_df)){
    warning('uv_dist_df will be ignored, since uv_dist is normal.')
  }
  
  # Error checking for valid degrees of freedom input if uv_dist=t.
  if(uv_dist=='t'){
    if(is.null(uv_dist_df)) stop('uv_dist_df cannot be null.  enter value 3 or greater.')
    if(!is.null(uv_dist_df) & uv_dist_df<1) stop('uv_dist_df must be 3 or greater.')
  }
  
  if ((r2_w_tau!=0 | r2_w_mu!=0) & !is.null(nonlin_frac)) {
    stop('Can\'t do both correlation with weight and nonlinearity')
  }
  
  if ((r2_w_tau!=0 | r2_w_mu!=0) & !weights) {
    stop('could do correlation with weights and no weights but that is dumb')
  }
  
  if (!is.null(r2_error)) {
    if(r2_error<0 | r2_error>1) {
      stop('r2_error must be [0,1]')
    }
    if (sigu_multiplier!=1 | sigy_multiplier!=1) {
      stop('r2_error only works with sigu/sigy=1')
    }
  }
  
  #- Setup -----------------------------------------------------------------------------------------
  
  # Read the driver file of baseline real values, depending on world (medicare vs education).
  # Sets real baseline values, which are them modified by multipliers input by the user.
  driver = read_excel(path='data/sim_study_config_medicare_and_edu.xlsx',
                      sheet=driver_sheet) %>%
    filter(sim_world==world)
  
  # Extract parameters from driver file
  sig_y_bene = driver %>% filter(name=='sigma_y')    %>% select(value) %>% pull()
  sig_u      = driver %>% filter(name=='sigma_u')    %>% select(value) %>% pull()
  sig_v      = driver %>% filter(name=='sigma_v')    %>% select(value) %>% pull()
  sigma_taux = driver %>% filter(name=='sigma_taux') %>% select(value) %>% pull()
  r2a        = driver %>% filter(name=='r2a')        %>% select(value) %>% pull()
  
  # Calculate total sample size
  N = (nT + nC)
  
  # Draw practice sizes and calculate weights.
  if(world=='medicare'){
    # Real practice sizes from CPC+.
    real_prac_size_q <- readRDS('data/RealDataPracticeSizeQuantiles.RDS')
    size = sample(quantile(ecdf(real_prac_size_q), seq(0,1,length=N)))
  } else{
    # Real # of students in traditional tested grades (3-8) across all schools with some students in those grades.
    real_edu_size_q <- readRDS('data/sch-size-quantiles.RDS') %>%
      select(q_tested) %>%
      pull()
    size = sample(quantile(ecdf(real_edu_size_q), seq(0,1,length=N)))
  }
  
  if(weights==TRUE){
    w = size # Use un-standardized weights, like we feed to iBCF for real data.
  }else{
    w=rep(mean(size),N)
    size=w
  }
  
  #- Random components -----------------------------------------------------------------------------
  
  # Generate random components.
  # These are bivariate normal (or t) draws with correlation rho, and df uv_dist_df (if t).
  if(uv_dist=='normal'){
    sig_uv = matrix(c(sig_u*sig_u, sig_u*sig_v*rho, sig_u*sig_v*rho, sig_v*sig_v), ncol=2, byrow=T)
    uv = rmvnorm(N, mean=c(0,0), sigma=sig_uv)
  } else{
    # If t, Var = (df/(df-2)*Scale matrix --> sig_uv is scale matrix, so multiply by (df-2)/df.
    sig_uv = matrix(c(sig_u*sig_u, sig_u*sig_v*rho, sig_u*sig_v*rho, sig_v*sig_v),
                    ncol=2, byrow=T) * (uv_dist_df-2)/uv_dist_df
    uv = rmvt(N, sigma=sig_uv, delta=c(0,0), df=uv_dist_df)
  }
  
  u = uv[,1]
  v = uv[,2]
  
  # Generate practice-level residuals.
  e = rnorm(N, 0, sig_y_bene/sqrt(w))
  
  #- Covariates and treatment assignment  ----------------------------------------------------------
  
  # Generate covariates.  All are standard normal.
  p = 5
  x = data.frame(matrix(rnorm(N*p), nrow=N))
  colnames(x) = paste0('x',1:p)
  
  # Propensity score function and treatment assignment.
  #For regular DGP, do this before mux/taux
  if (is.null(nonlin_frac)) {
    if(targeted_sel==TRUE){
      # q creates targeted selection using control covariates x1, x2.
      # Calibrate intercept approximately s.t. mean of propensity score aligns with nT vs nC ratio.
      q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) +  runif(N,-.5,.5)/10
    } else {
      q = sample(c(-1,1), size=N, replace=TRUE) +  runif(N,-.5,.5)/10
    }
    
    if (r2_w_pi!=0) {
      #Just do a basic step function
      q_w <- as.numeric(cut(w,c(0,quantile(w,c(.25,.5,.75)), Inf)))
      q <- sqrt(1-r2_w_pi)*q/sd(q) + sqrt(r2_w_pi) * (q_w-mean(q_w))/sd(q_w)
    }
    
    #Why is this NC+NC?
    # Empirically setting mean q to 1/4 gets mean p close to 1/3 but I don't think that's mechanical
    q = q - mean(q) + qnorm(nT/(nC+nC))
    pi =  pnorm(q)
    #Select treatment with probability proportional to pscore
    #NB: I think this means the scaling of the pscore doesn't really matter?
    z <- rep(0, N)
    z[sample(1:N, nT, prob=pi)] <- 1
  }
  
  # - Mux and taux definition and calibration ------------------------------------------------------
  
  #Base DGP
  if (is.null(nonlin_frac)) {
    # Define mux function. Uses control covariates x1, x2, x5.
    mux = ifelse(x[,2]>0,-1,1)*6 + abs(x[,1]-1) + 3*x[,5]
    
    if(r2_w_mu!=0) {
      mu_w <- as.numeric(cut(w,c(0,quantile(w,c(.25,.5,.75)), Inf)))
      mux <- sqrt(1-r2_w_mu)*mux + sqrt(r2_w_mu)*sd(mux)*mu_w/sd(mu_w)
    }
    
    # Define treatment effects function taux. For het, uses treatment covariates x3, x4.
    # For het, center tau so that it's sensible to use 0 as a threshold for P(favorable impact).
    if(trt_eff_scenario=='homog'){
      taux = rep(3, N)
    } else {
      taux = 1 + x[,3] + .75*x[,4]
      taux <- taux * sigma_taux/sd(taux)
    }
    
    if(r2_w_tau!=0) {
      tau_w <- as.numeric(cut(w,c(0,quantile(w,c(.25,.5,.75)), Inf)))
      taux <- sqrt(1-r2_w_tau)*taux + sqrt(r2_w_tau)*sd(taux)*tau_w/sd(tau_w)
    }
    
  } else {
    #DGP with separate linear/nonlinear portions
    lin_mux <- x[,1] -2*x[,2] + 3*x[,5]
    #Stupid functions designed to each have about the same SD, and to have no linear R2
    #abs(x-1) for example does have a pretty linear relationship, because most of the mass is <1
    nonlin_mux <- abs(x[,1]) + 
      .5*(x[,2] < qnorm(.25)) + 
      -1*(x[,2] > qnorm(.25) & x[,2] < qnorm(.5)) +
      1*(x[,2] > qnorm(.5) & x[,2] < qnorm(.75)) +
      -.5*(x[,2] > qnorm(.75)) +
      (x[,5]+.5)^3/12 + (x[,5]-.5)^2/3
    intxn_mux <- x[,1]*x[,2] + 2*(x[,2]+x[,5]>0 & x[,2]+x[,5]<1)
    
    #Standardize just so we can (roughly) control what proportion is linear vs not
    lin_mux    <- (lin_mux - mean(lin_mux)) / sd(lin_mux)
    nonlin_mux <- (nonlin_mux - mean(nonlin_mux)) / sd(nonlin_mux)
    intxn_mux  <- (intxn_mux - mean(intxn_mux)) / sd(intxn_mux)
    
    #Actual SD of mux doesn't matter since it later gets scaled to match r2_a, but get it close to start
    v_resid <- sig_u^2 + sig_v^2 + sig_y_bene^2/mean(size)
    approx_sd_mux <- sqrt((r2a*v_resid + (r2a-1)*sigma_taux^2) / (1-r2a))
    #Nonlinearity should be 2:1 from main effects vs interactions
    mux <- approx_sd_mux * (sqrt(1-nonlin_frac)*lin_mux + sqrt(2*nonlin_frac/3)*nonlin_mux + sqrt(nonlin_frac/3)*intxn_mux)
    
    if(trt_eff_scenario=='homog'){
      taux = rep(3, N)
    } else {
      lin_taux <- x[,3] + .75*x[,4]
      nonlin_taux <- -.5*(x[,3] < qnorm(.25)) + 
        1*(x[,3] > qnorm(.25) & x[,3] < qnorm(.5)) +
        -1*(x[,3] > qnorm(.5) & x[,3] < qnorm(.75)) +
        .5*(x[,3] > qnorm(.75)) +
        1/(abs(x[,4])+.25)
      intxn_taux <- x[,3]*x[,4]
      
      lin_taux    <- (lin_taux - mean(lin_taux)) / sd(lin_taux)
      nonlin_taux <- (nonlin_taux - mean(nonlin_taux)) / sd(nonlin_taux)
      intxn_taux  <- (intxn_taux - mean(intxn_taux)) / sd(intxn_taux)
      taux <- sqrt(1-nonlin_frac)*lin_taux + sqrt(2*nonlin_frac/3)*nonlin_taux + sqrt(nonlin_frac/3)*intxn_taux
      
      taux <- 3 + (taux-mean(taux)) * sigma_taux/sd(taux)
    }
  }
  
  #For nonlinear DGP, do treatment assignment based on mux
  if (!is.null(nonlin_frac)) {
    if(targeted_sel==TRUE){
      # q creates targeted selection
      # First take std mu, then add random noise
      q = ((mux-mean(mux))/sd(mux))/2 + runif(N,-1,1)
    } else {
      q = sample(c(-1,1), size=N, replace=TRUE) +  runif(N,-.5,.5)/10
    }
    #Why is this NC+NC?
    # Empirically setting mean q to 1/4 gets mean p close to 1/3 but I don't think that's mechanical
    q = q - mean(q) + qnorm(nT/(nT+nC))
    pi =  pnorm(q)
    #Select treatment with probability proportional to pscore
    #NB: I think this means the scaling of the pscore doesn't really matter?
    z <- rep(0, N)
    z[sample(1:N, nT, prob=pi)] <- 1
  }
  
  # Calibrate mux to get desired r2a for trt obs.
  # r2a = 1 - RSS / SST = 1 - var(resids)/var(y),
  # where y = mux*multiplier + u + (taux+v)*z + e.
  y = mux + u + (taux + v) * z + e # Set y with the current mux and taux, to calc resids.
  resids = u + v*z + e # Calculate reids.
  
  multiplier=1
  temp_r2a = 1 - weighted.var(resids[z==1], w[z==1]) / weighted.var(y[z==1], w[z==1])
  
  # If r2a is smaller than desired, increase multiplier until reach desired r2a.
  # If r2a is larger than desired, decrease multiplier until reach desired r2a.
  # Fast mode is useful for testing, but it is decently imprecise (e.g. IQR is 47-52% when targeting 50)
  # I think that's due to incidental correlation between terms
  if (fast_r2a) {
    var_resid = var(u[z==1] + v[z==1] + e[z==1])
    var_explained = var_resid * r2a / (1-r2a)
    sigma_mux = sqrt(var_explained - sigma_taux^2)
    if (sigma_mux <=0) stop ('something funky in fast_r2a')
    mux <- mux * sigma_mux/sd(mux)
  } else {
    #Precalculate so loops are faster
    idx_trt <- z==1
    w_trt <- w[idx_trt]
    res_trt <- resids[idx_trt]
    wv_resid <- weighted.var(res_trt, w_trt)
    mux_trt <- mux[idx_trt]
    y_trt <- y[idx_trt]
    y_less_mux_trt <- (u + (taux+v)*z + e)[idx_trt]
    
    if(temp_r2a < r2a){
      while(temp_r2a < r2a){
        multiplier = multiplier + .001  # Increment multiplier
        y_trt = mux_trt*multiplier + y_less_mux_trt # Calculate new y with multiplier
        temp_r2a = 1 - wv_resid / weighted.var(y_trt, w_trt) # Recalc r2a
      }
    } else if(temp_r2a > r2a){
      while(temp_r2a > r2a){
        multiplier = multiplier - .001  # Increment multiplier
        y_trt = mux_trt*multiplier + y_less_mux_trt # Calculate new y with multiplier
        temp_r2a = 1 - wv_resid / weighted.var(y_trt, w_trt) # Recalc r2a
      }
    } else{
      multiplier=1
    }
  }
  
  mux = multiplier*mux
  
  #Now that everything has been scaled correctly relative to one another, apply our options for differential scaling
  #doing this after calibration, so e.g. we don't have an increase in sigma_u also scaling up sd(mux) to maintain r2
  if (!is.null(r2_error)) {
    avg_w <- mean(w)
    #When using r2_error the idea is to keep the average total variance the same
    tot_var <- sig_u^2 + sig_y_bene^2/avg_w
    sigu_multiplier <- sqrt(r2_error*tot_var)/sig_u
    sigy_multiplier <- sqrt((1-r2_error)*tot_var*avg_w)/sig_y_bene
  }
  
  e    <-    e * sigy_multiplier
  u    <-    u * sigu_multiplier
  v    <-    v * sigv_multiplier
  taux <- taux * sigt_multiplier
  
  #- Assemble and return data frame ----------------------------------------------------------------
  
  # Add means to mux and taux.
  if(world=='medicare'){
    mux = mux - mean(mux) + 70    # Represents annual increasing cost of medical care on avg.
    
    if(center_ate==FALSE){
      taux = taux - mean(taux) + 20 # ATE = 20 (negative for decreasing expenditures)
    }
  } else{
    mux = mux - mean(mux)
    if(center_ate==FALSE){
      taux = taux - mean(taux) + .2 # .04 is value that gives similar effect size to Medicare, for testing.
    }
  }
  
  # Update y to use the adjusted mux and taux.
  y = mux + u + (taux + v) * z + e
  
  # Assemble dataframe.
  d = data.frame('id' = 1:N,
                 x, size, w, pi,
                 'mux' = mux, u, 'mu'=mux+u,
                 'z_str' = ifelse(z==1,'trt','ctrl'), 'z'=as.numeric(z),
                 'taux'=taux, v, 'tau'=taux+v,
                 y, e,
                 'sigma_u' = as.numeric(sig_u)*sigu_multiplier,
                 'sigma_v' = sig_v*sigv_multiplier,
                 'rho' = rho,
                 'sigma_y' = sig_y_bene*sigy_multiplier,
                 'sigma_taux' = sigma_taux*sigt_multiplier,
                 'r2a_input' = r2a)
  
  # Check the actual r2 values in the generated dataset.
  r2_check = dataSummy(d)
  d$r2a_check = r2_check$r2a
  d$r2uv_check = r2_check$r2uv
  d$r2trt_check = r2_check$r2trt
  
  # Use reparam parameterization
  if(reparam){
    
    # Calculate reparam sigma_u and sigma_v, and ui and vi.
    d = d %>%
      mutate(
        a = u,
        b = u+v,
        sigma_a = sigma_u,
        sigma_b = sqrt(sigma_a^2 + sigma_v^2 + 2*rho*sigma_a*sigma_v)
      )
  }
  
  row.names(d)=NULL
  return(d)
}


#- Utility function to summarize calibration of generated dataset ----------------------------------

dataSummy = function(data){
  
  # Residual R-squared for treated.  (What we were calling r2a.)
  r2a = 1 - weighted.var(data$u[data$z==1] + data$v[data$z==1] + data$e[data$z==1], data$w[data$z==1]) /
    weighted.var(data$y[data$z==1], data$w[data$z==1])
  
  #r2uv.  sig2v / (sig2v + sig2u)
  r2uv = weighted.var(data$v, data$w) / (weighted.var(data$v, data$w) + weighted.var(data$u, data$w))
  
  #r2trt.  var(taux) / (var(taux) + sig2v), the proportion of trt var explained by x's.
  r2trt = weighted.var(data$taux[data$z==1], data$w[data$z==1]) /
    (weighted.var(data$taux[data$z==1], data$w[data$z==1]) + weighted.var(data$v[data$z==1], data$w[data$z==1]))
  
  out = cbind.data.frame(
    'r2a' = r2a,
    'r2uv' = r2uv,
    'r2trt' = r2trt) %>%
    mutate(across(where(is.numeric), round,3))
  
  
  return(out)
  
}
