#' @importFrom stats approxfun lm qchisq quantile sd
#' @importFrom RcppParallel RcppParallelLibs
Rcpp::loadModule(module = "TreeSamples", TRUE)

.ident <- function(...){
    # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
    args <- c(...)
    if( length( args ) > 2L ){
        #  recursively call ident()
        out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
    }else{
        out <- identical( args[1] , args[2] )
    }
    return( all( out ) )
}

.cp_quantile = function(x, num=10000, cat_levels=8){
    nobs = length(x)
    nuniq = length(unique(x))

    if(nuniq==1) {
        ret = x[1]
        warning("A supplied covariate contains a single distinct value.")
    } else if(nuniq < cat_levels) {
        xx = sort(unique(x))
        ret = xx[-length(xx)] + diff(xx)/2
    } else {
        q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
        ind = seq(min(x),max(x),length.out=num)
        ret = q(ind)
    }

    return(ret)
}

.get_chain_tree_files = function(tree_path, chain_id){
    if (is.null(tree_path)){
        out <- list(
            "con_trees" = toString(character(0)),
            "mod_trees" = toString(character(0))
        )
    } else{
        out <- list("con_trees" = paste0(tree_path,'/',"con_trees.", chain_id, ".txt"),
                    "mod_trees" = paste0(tree_path,'/',"mod_trees.", chain_id, ".txt"))
    }
    return(out)
}

.get_do_type = function(n_cores, log_file){
    if(n_cores>1){
        cl <- parallel::makeCluster(n_cores, outfile=log_file)

        cat(sprintf("Running in parallel, saving BCF logs to %s \n", log_file))
        doParallel::registerDoParallel(cl)
        `%doType%`  <- foreach::`%dopar%`
    } else {
        cl <- NULL
        `%doType%`  <- foreach::`%do%`
    }

    do_type_config <- list('doType'  = `%doType%`,
                           'n_cores' = n_cores,
                           'cluster' = cl)

    return(do_type_config)
}

.cleanup_after_par = function(do_type_config){
    if(do_type_config$n_cores>1){
        parallel::stopCluster(do_type_config$cluster)
    }
}

.get_info_from_chains = function(chains) {
    list(nchain = length(chains),
         ndraw  = nrow(chains[[1]]$yhat),
         nobs   = ncol(chains[[1]]$yhat),
         abcf   = chains[[1]]$abcf,
         ibcf   = chains[[1]]$ibcf)
}

.extract_matrix_from_chains = function(chains, par) {
    do.call(what = rbind, args = lapply(chains, `[[`, par))
}

.extract_vector_from_chains = function(chains, par) {
    unlist(lapply(chains, `[[`, par))
}

.extract_value_from_chains = function(chains, par) {
    chains[[1]][[par]]
}

.extract_coda_chains = function(chains) {
    info = .get_info_from_chains(chains)
    mcmcs <- lapply(chains, function(x) {

        scalars <- data.frame("tau_bar"   = matrixStats::rowWeightedMeans(x$tau,  x$w),
                              "mu_bar"    = matrixStats::rowWeightedMeans(x$mu,   x$w),
                              "yhat_bar"  = matrixStats::rowWeightedMeans(x$yhat, x$w),
                              "mu_scale"  = x$mu_scale,
                              "tau_scale" = x$tau_scale,
                              "b0"        = x$b0,
                              "b1"        = x$b1,
                              "delta_mu"  = x$delta_mu)
        if (info$ibcf) {
            addl_scalars <- data.frame("sigma_y"   = x$sigma_y,
                                       "sigma_u"   = x$sigma_u,
                                       "sigma_v"   = x$sigma_v,
                                       "rho"       = x$rho)
        } else if (info$abcf) {
            addl_scalars <- data.frame("sigma_y"   = x$sigma_y,
                                       "sigma_u"   = x$sigma_u)
        } else {
            addl_scalars <- data.frame("sigma"     = x$sigma)
        }
        return(coda::as.mcmc(cbind(scalars, addl_scalars)))
    })
    return(coda::as.mcmc.list(mcmcs))
}

.get_components_from_chains <- function(chains) {
    info = .get_info_from_chains(chains)
    ret <- list()

    #Matrices
    ret$yhat <- .extract_matrix_from_chains(chains, 'yhat')
    ret$mu   <- .extract_matrix_from_chains(chains, 'mu')
    ret$tau  <- .extract_matrix_from_chains(chains, 'tau')

    if (info$ibcf) {
        ret$u             <- .extract_matrix_from_chains(chains, 'u')
        ret$v             <- .extract_matrix_from_chains(chains, 'v')
        ret$sigma_y       <- .extract_vector_from_chains(chains, 'sigma_y')
        ret$sigma_u       <- .extract_vector_from_chains(chains, 'sigma_u')
        ret$sigma_v       <- .extract_vector_from_chains(chains, 'sigma_v')
        ret$rho           <- .extract_vector_from_chains(chains, 'rho')
        #This is super esoteric, and also giant since it's a N*draws matrix, so not actually saving
        #Ditto the unfixed uv total residuals
        #ret$sigma_i       <- .extract_matrix_from_chains(chains, 'sigma_i')
        ret$acceptance    <- .extract_matrix_from_chains(chains, 'acceptance')
        rownames(ret$acceptance) <- paste('chain',1:info$nchain,sep='_')
    } else if (info$abcf) {
        ret$u             <- .extract_matrix_from_chains(chains, 'u')
        ret$sigma_y       <- .extract_vector_from_chains(chains, 'sigma_y')
        ret$sigma_u       <- .extract_vector_from_chains(chains, 'sigma_u')
        ret$acceptance    <- .extract_matrix_from_chains(chains, 'acceptance')
        rownames(ret$acceptance) <- paste('chain',1:info$nchain,sep='_')
    } else {
        ret$sigma         <- .extract_vector_from_chains(chains, 'sigma')
    }

    ret$mu_scale  <- .extract_vector_from_chains(chains, 'mu_scale')
    ret$delta_mu  <- .extract_vector_from_chains(chains, 'delta_mu')
    ret$tau_scale <- .extract_vector_from_chains(chains, 'tau_scale')
    ret$b0        <- .extract_vector_from_chains(chains, 'b0')
    ret$b1        <- .extract_vector_from_chains(chains, 'b1')

    for (par in c('sdy','con_sd','mod_sd','muy','y','z','w','perm','include_pi','abcf', 'ibcf', 'random_seed')) {
        ret[[par]] <- .extract_value_from_chains(chains, par)
    }

    if (info$ibcf) {
        ret[['ate_prior_sd']] <- .extract_value_from_chains(chains, 'ate_prior_sd')
    }

    return(ret)
}

#' Fit Bayesian Causal Forests
#'
#' @references Hahn, Murray, and Carvalho(2017). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.
#'  https://arxiv.org/abs/1706.09523. (Call citation("bcf") from the
#' command line for citation information in Bibtex format.)
#'
#' @details Fits the Bayesian Causal Forest model (Hahn et. al. 2018): For a response
#' variable y, binary treatment z, and covariates x, we return estimates of mu, tau, and sigma in
#' the model
#' \deqn{y_i = \mu(x_i, \pi_i) + \tau(x_i, \pi_i)z_i + \epsilon_i}
#' where \eqn{\pi_i} is an (optional) estimate of the propensity score \eqn{\Pr(Z_i=1 | X_i=x_i)} and
#' \eqn{\epsilon_i \sim N(0,\sigma^2)}
#'
#' Some notes:
#' \itemize{
#'    \item x_control and x_moderate must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame
#'    \item sd_control and sd_moderate are the prior SD(mu(x)) and SD(tau(x)) at a given value of x (respectively). If
#'    use_muscale = FALSE, then this is the parameter \eqn{\sigma_\mu} from the original BART paper, where the leaf parameters
#'    have prior distribution \eqn{N(0, \sigma_\mu/m)}, where m is the number of trees.
#'    If use_muscale=TRUE then sd_control is the prior median of a half Cauchy prior for SD(mu(x)). If use_tauscale = TRUE,
#'    then sd_moderate is the prior median of a half Normal prior for SD(tau(x)).
#'    \item By default the prior on \eqn{\sigma^2} is calibrated as in Chipman, George and McCulloch (2008).
#'
#'
#' }
#' @param y Response variable
#' @param z Treatment variable
#' @param x_control Design matrix for the prognostic function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param pihat Length n estimates of propensity score
#' @param w An optional vector of weights. When present, BCF fits a model \eqn{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function.
#' @param random_seed A random seed passed to R's set.seed
#' @param n_chains  An optional integer of the number of MCMC chains to run
#' @param n_cores An optional integer of the number of cores to run your MCMC chains on
#' @param n_threads An optional integer of the number of threads to parallelize within chain bcf operations on. Values greater than 1 tend to reduce performance, unless you have truly massive within-chain, within-iteration calculations, like a ginormous dataset
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in. The chain will run for nsim*nthin iterations after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param ntree_moderate Number of trees in tau(x)
#' @param sd_moderate SD(tau(x)) marginally at any covariate value (or its prior median if use_tauscale=TRUE)
#' @param base_moderate Base for tree prior on tau(x) trees (see details)
#' @param power_moderate Power for the tree prior on tau(x) trees (see details)
#' @param save_tree_directory Specify where trees should be saved. Keep track of this for predict(). Defaults to working directory. Setting to NULL skips writing of trees.
#' @param continuous_tree_save Whether trees are written continuously or only once at the end
#' @param log_file file where BCF should save its logs when running multiple chains in parallel. This file is not written to when only running one chain.
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to
#' include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#' @param simplified_return Whether to return just the raw_chains object (which contains all relevant information), or to pre-calculate other output for ease of use (e.g. fit$tau).
#' All output can be recreated by calling abcf:::.get_components_from_chains(fit$raw_chains)
#' @param verbose Integer, whether to print log of MCMC iterations, defaults to 1 - basic logging of iteration progress.
#' Setting to 0 disables logging, while setting to 2 enables logging of detailed statistics each iteration,
#' and setting to 3 enables logging of individual trees.
#' @param block_b0_b1 Whether to constrain b0 and b1 such that b0 = -b1. Better mixing at the expense of assuming equal variance across treatment and control
#' @param abcf Boolean; whether to estimate the aggregate BCF (aBCF) model, including individual random effects
#' @param sigu_hyperprior Prior median for prior SD of idiosyncratic u terms in aBCF and iBCF
#' @param batch_size aBCF and iBCF only. Size of batches for adaptive Metropolis Hastings sampling of sigma_u (aBCF) as well as sigma_v/rho (iBCF only)
#' @param ibcf Boolean; whether to estimate the individualized BCF (iBCF) model, including individual random effects and random treatment effects.
#' @param ate_prior_sd Prior SD of the treatment effect, used to form a hyperprior for the SD of idiosyncratic treatment effects
#' @return A fitted bcf object that is a list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual-level treatment effect estimates}
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of prognostic function E(Y|Z=0, x=x) estimates}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#'
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=2000, nsim=2000)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'
#'}
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' pihat = pnorm(q)
#'
#' # nburn and nsim should be much larger, at least a few thousand each
#' # The low values below are for CRAN.
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=100, nsim=10)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'}
#'
#' @useDynLib abcf
#' @export
bcf <- function(y, z, x_control, x_moderate=x_control, pihat, w = NULL,
                random_seed = sample.int(.Machine$integer.max, 1),
                n_chains = 4,
                n_cores  = n_chains,
                n_threads = 1, #max number of threads, minus a arbitrary holdback, over the number of cores. Tends to be much slower to use a value over 1, unless there are massive within-iteration calculations
                nburn, nsim, nthin = 1, update_interval = 100,
                ntree_control = 200,
                sd_control = NULL,      #Still defaults to 2*sdy, but we do that latter so we can account for weights
                base_control = 0.95,
                power_control = 2,
                ntree_moderate = 50,
                sd_moderate = NULL,     #Still defaults to sdy, but we do that latter so we can account for weights
                base_moderate = 0.25,
                power_moderate = 3,
                save_tree_directory = '.',
                keep_trees = FALSE,
                continuous_tree_save=FALSE,
                log_file=file.path('.',sprintf('bcf_log_%s.txt',format(Sys.time(), "%Y%m%d_%H%M%S"))),
                nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                include_pi = "control", use_muscale=TRUE, use_tauscale=TRUE,
                simplified_return=FALSE,
                verbose=1,
                block_b0_b1=FALSE,
                abcf=TRUE,
                sigu_hyperprior = NULL,
                batch_size = 100,
                ibcf=FALSE,
                ate_prior_sd = NULL
) {


    if(is.null(w)){
        w <- matrix(1, ncol = 1, nrow = length(y))
    }

    pihat = as.matrix(pihat)
    if(!.ident(length(y),
               length(z),
               length(w),
               nrow(x_control),
               nrow(x_moderate),
               nrow(pihat))
    ) {
        stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
             "length(z): ", length(z), "\n",
             "length(w): ", length(w), "\n",
             "nrow(x_control): ", nrow(x_control), "\n",
             "nrow(x_moderate): ", nrow(x_moderate), "\n",
             "nrow(pihat): ", nrow(pihat),"\n"
        )
    }

    if(any(is.na(y))) stop("Missing values in y")
    if(any(is.na(z))) stop("Missing values in z")
    if(any(is.na(w))) stop("Missing values in w")
    if(any(is.na(x_control))) stop("Missing values in x_control")
    if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
    if(any(is.na(pihat))) stop("Missing values in pihat")
    if(any(!is.finite(y))) stop("Non-numeric values in y")
    if(any(!is.finite(z))) stop("Non-numeric values in z")
    if(any(!is.finite(w))) stop("Non-numeric values in w")
    if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
    if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
    if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")
    if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")
    if(!(keep_trees %in% c(TRUE,FALSE))) stop("keep_trees must be TRUE or FALSE")
    if(!(continuous_tree_save %in% c(TRUE,FALSE))) stop("continuous_tree_save must be TRUE or FALSE")
    if (keep_trees & !is.null(save_tree_directory)) stop("Can\'t both save trees and keep trees; set save_tree_directory to NULL")
    if (keep_trees & continuous_tree_save) stop("Can\'t both keep trees and write them continuously; set continuous_tree_save to FALSE")
    if(!use_tauscale & block_b0_b1) stop('Can\'t block b0 and b1 if tauscale is not used')
    if(!(abcf %in% c(TRUE,FALSE))) stop("abcf must be TRUE or FALSE")
    if(!(ibcf %in% c(TRUE,FALSE))) stop("ibcf must be TRUE or FALSE")
    if(abcf|ibcf) {
        if(round(batch_size)!=batch_size | batch_size<1) stop("batch_size must be an integer larger than 0")
    }
    if (ibcf) {
        if(is.null(ate_prior_sd)) stop("must supply ate_prior_sd when using iBCF")
        if(!is.numeric(ate_prior_sd) | ate_prior_sd<=0) stop("ate_prior_sd must be a positive number")
    }
    if(abcf & ibcf) {
        stop('Can\'t do both aBCF and iBCF')
    }
    if(!(verbose %in% 0:4)) stop("verbose must be an integer from 0 to 4")

    if(length(unique(y))<5) warning("y appears to be discrete")

    if(nburn<0) stop("nburn must be positive")
    if(nsim<0) stop("nsim must be positive")
    if(nthin<0) stop("nthin must be positive")
    if(nthin>nsim+1) stop("nthin must be < nsim")
    if(nburn<1000) warning("A low (<1000) value for nburn was supplied")
    if(nsim*nburn<1000) warning("A low (<1000) value for total iterations after burn-in was supplied")
    if (use_tauscale & !identical(x_control, x_moderate)) {
        warning("Different covariate matrices supplied to x_control and x_moderate, but tau_scale is set to TRUE. When use_tauscale is TRUE, all covariates in x_moderate can still affect mu (but covariates in x_control cannot affect tau)")
    }
    if ((abcf|ibcf) & length(unique(w))==1) {
        warning('aBCF and iBCF models are not identified without weights')
    }


    ### TODO range check on parameters

    ###
    x_c = matrix(x_control, ncol=ncol(x_control))
    x_m = matrix(x_moderate, ncol=ncol(x_moderate))

    if(include_pi=="both" | include_pi=="control") {
        x_c = cbind(x_control, pihat)
    }
    if(include_pi=="both" | include_pi=="moderate") {
        x_m = cbind(x_moderate, pihat)
    }
    cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
    cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))

    sdy = sqrt(Hmisc::wtd.var(y, w))
    muy = stats::weighted.mean(y, w)
    yscale = (y-muy)/sdy

    if(is.null(lambda)) {
        if(is.null(sighat)) {
            lmf = lm(yscale~z+as.matrix(x_c), weights = w)
            sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
        }
        qchi = qchisq(1.0-sigq,nu)
        lambda = (sighat*sighat*qchi)/nu
    }

    #If prior sds are not given, default them to scale off of the weighted sdy
    if (is.null(sd_control)) {
        con_sd <- 2
    } else {
        con_sd = sd_control/sdy
    }

    if (is.null(sd_moderate)) {
        mod_sd <- 1/ifelse(use_tauscale,0.674,1)
    } else {
        mod_sd = sd_moderate/sdy/ifelse(use_tauscale,0.674,1)
    }

    #If hyperprior sd isn't given, scale them off of the prior sds
    if (is.null(sigu_hyperprior)) {
        sigu_hyperprior <- con_sd/3
    } else {
        #user-entered values will should be prior medians, so convert to scale using 0.674
        sigu_hyperprior <- sigu_hyperprior/sdy/0.674
    }

    if (ibcf) {
        ate_prior_sd <- ate_prior_sd/sdy
    } else {
        ate_prior_sd <- 1
    }

    dir = tempdir()

    perm = order(z, decreasing=TRUE)

    RcppParallel::setThreadOptions(numThreads=n_threads)

    do_type_config <- .get_do_type(n_cores, log_file)
    `%doType%` <- do_type_config$doType

    chain_out <- foreach::foreach(iChain=1:n_chains) %doType% {

        this_seed = random_seed + iChain - 1

        cat("Calling bcfoverparRcppClean From R\n")
        set.seed(this_seed)

        tree_files = .get_chain_tree_files(save_tree_directory, iChain)

        fitbcf = bcfoverparRcppClean(y_ = yscale[perm], z_ = z[perm], w_ = w[perm],
                                     x_con_ = t(x_c[perm,,drop=FALSE]), x_mod_ = t(x_m[perm,,drop=FALSE]),
                                     x_con_info_list = cutpoint_list_c,
                                     x_mod_info_list = cutpoint_list_m,
                                     burn = nburn, nd = nsim, thin = nthin,
                                     ntree_mod = ntree_moderate, ntree_con = ntree_control,
                                     lambda = lambda, nu = nu,
                                     con_sd = con_sd,
                                     mod_sd = mod_sd, # if HN make sd_moderate the prior median
                                     mod_alpha = base_moderate,
                                     mod_beta = power_moderate,
                                     con_alpha = base_control,
                                     con_beta = power_control,
                                     treef_con_name_ = tree_files$con_trees,
                                     treef_mod_name_ = tree_files$mod_trees,
                                     keep_trees = keep_trees,
                                     continuous_tree_save=continuous_tree_save,
                                     status_interval = update_interval,
                                     use_mscale = use_muscale, use_bscale = use_tauscale,
                                     b_half_normal = TRUE,
                                     abcf=abcf, ibcf=ibcf,
                                     batch_size=batch_size, acceptance_target=0.44,
                                     verbose=verbose,
                                     block_b0_b1=block_b0_b1,
                                     sigu_hyperprior=sigu_hyperprior,
                                     ate_prior_sd=ate_prior_sd)

        cat("bcfoverparRcppClean returned to R\n")

        ac = fitbcf$m_post[,order(perm)]

        Tm = fitbcf$b_post[,order(perm)] * (1.0/ (fitbcf$b1 - fitbcf$b0))

        Tc = ac * (1.0/fitbcf$msd)

        tau_post = sdy*fitbcf$b_post[,order(perm)]

        mu_post  = muy + sdy*(Tc*fitbcf$msd + Tm*fitbcf$b0)

        yhat_post = muy + sdy*fitbcf$yhat_post[,order(perm)]

        u_post = sdy*fitbcf$u[,order(perm)]

        v_post = sdy*fitbcf$v[,order(perm)]

        if (abcf) {
            mu_post    <- mu_post   + u_post
            yhat_post  <- yhat_post + u_post
        } else if (ibcf) {
            tau_post   <- tau_post  + v_post
            mu_post    <- mu_post   + u_post
            yhat_post  <- yhat_post + u_post + t(t(v_post) * z)
        }

        sigma_i = sdy*fitbcf$sigma_i[,order(perm)]

        names(fitbcf$acceptance) = c('sigma_y','sigma_u','sigma_v','rho')

        list(sigma_y    = sdy*fitbcf$sigma_y,
             sigma_u    = sdy*fitbcf$sigma_u,
             sigma_v    = sdy*fitbcf$sigma_v,
             rho        = fitbcf$rho,
             sigma_i    = sigma_i,
             yhat       = yhat_post,
             sdy        = sdy,
             con_sd     = con_sd,
             mod_sd     = mod_sd,
             muy        = muy,
             mu         = mu_post,
             tau        = tau_post,
             u          = u_post,
             v          = v_post,
             mu_scale   = fitbcf$msd,
             tau_scale  = fitbcf$bsd,
             b0         = fitbcf$b0,
             b1         = fitbcf$b1,
             delta_mu   = fitbcf$delta_con,
             acceptance = fitbcf$acceptance,
             perm       = perm,
             include_pi = include_pi,
             abcf = abcf,
             ibcf = ibcf,
             ate_prior_sd = sdy*ate_prior_sd,
             random_seed=this_seed,
             con_trees = fitbcf$con_trees,
             mod_trees = fitbcf$mod_trees
        )

    }

    if (!ibcf & !abcf) {
        #If we're not running IBCF, remove all the iBCF-specific components
        chain_out <- lapply(chain_out, function(x) {
            x$sigma <- x$sigma_y
            x$sigma_y <- x$sigma_u <- x$sigma_v <- x$rho <- x$sigma_i <- x$u <- x$v <- x$acceptance <- x$ate_prior_sd <- NULL
            return(x)
        })
    } else if (!ibcf) {
        chain_out <- lapply(chain_out, function(x) {
            x$sigma_v <- x$rho <- x$v <- x$ate_prior_sd <- NULL
            x$acceptance <- x$acceptance[c('sigma_y', 'sigma_u')]
            return(x)
        })
    }

    if (!keep_trees) {
        chain_out <- lapply(chain_out, function(x) {
            x$con_trees <- x$mod_trees <- NULL
            return(x)
        })
    }

    fitObj <- list(raw_chains = chain_out)

    if (!simplified_return) {
        fitObj <- c(fitObj,
                    list(coda_chains = .extract_coda_chains(chain_out)),
                    .get_components_from_chains(chain_out))
    }

    attr(fitObj, "class") <- "bcf"

    .cleanup_after_par(do_type_config)

    return(fitObj)
}

#' Summarize method for bcf() fitted object
#'
#' Takes a fitted bcf object produced by bcf() and produces summary stats and MCMC diagnostics.
#' This function is built using the coda package and meant to mimic output from rstan::print.stanfit().
#' It includes, for key parameters, posterior summary stats, effective sample sizes,
#' and Gelman and Rubin's convergence diagnostics.
#' By default, those parameters are: sigma (the error standard deviation when the weights
#' are all equal), tau_bar (the estimated sample average treatment effect), mu_bar
#' (the average outcome under control/z=0 across all observations in the sample), and
#' yhat_bat (the average outcome under the realized treatment assignment across all
#' observations in the sample).
#'
#' We strongly suggest updating the coda package to our
#' Github version, which uses the Stan effective size computation.
#' We found the native coda effective size computation to be overly optimistic in some situations
#' and are in discussions with the coda package authors to change it on CRAN.
#' @param object output from a BCF predict run.
#' @param ... additional arguments affecting the summary produced.
#' @param params_2_summarise parameters to summarise.
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#'
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=2000, nsim=2000)
#'
#' # Get model fit diagnostics
#' summary(bcf_fit)
#'
#'}
#' @export
summary.bcf <- function(object,
                        ...,
                        params_2_summarise = NULL){

    if (!is.null(params_2_summarise)) {

    } else if (is.null(params_2_summarise) & is.null(object$include_random_effects)) {
        params_2_summarise <- c('sigma','tau_bar','mu_bar','yhat_bar')
    } else if (is.null(params_2_summarise) & !object$include_random_effects) {
        params_2_summarise <- c('sigma','tau_bar','mu_bar','yhat_bar')
    } else if (is.null(params_2_summarise) & object$include_random_effects) {
        params_2_summarise <- c('sigma_y','sigma_u','sigma_v','rho','tau_bar','mu_bar','yhat_bar')
    }

    chains_2_summarise <- object$coda_chains[,params_2_summarise]

    message("Summary statistics for each Markov Chain Monte Carlo run")
    print(summary(chains_2_summarise))

    cat("\n----\n\n")


    message("Effective sample size for summary parameters")

    ef = function(e) {
        if(e$message == "unused argument (crosschain = TRUE)") {
            cat("Reverting to coda's default ESS calculation. See ?summary.bcf for details.\n\n")
            print(coda::effectiveSize(chains_2_summarise))
        } else {
            stop(e)
        }
    }
    tryCatch(print(coda::effectiveSize(chains_2_summarise, crosschain = TRUE)),
             error = ef)
    cat("\n----\n\n")


    if (length(chains_2_summarise) > 1){
        message("Gelman and Rubin's convergence diagnostic for summary parameters")
        print(coda::gelman.diag(chains_2_summarise, autoburnin = FALSE))
        cat("\n----\n\n")

    }

}


#' Takes a fitted bcf object produced by bcf() and produces predictions for a new set of covariate values
#'
#' This function takes in an existing BCF model fit and uses it to predict estimates for new data.
#' It is important to note that this function requires that you indicate where the trees from the model fit are saved.
#' You can do so using the save_tree_directory argument in bcf(). Otherwise, they will be saved in the working directory.
#' The bcf() function automatically saves those in the same directory as the
#' @param object output from a BCF predict run
#' @param ... additional arguments affecting the predictions produced.
#' @param x_predict_control matrix of covariates for the "prognostic" function mu(x) for predictions (optional)
#' @param x_predict_moderate matrix of covariates for the covariate-dependent treatment effects tau(x) for predictions (optional)
#' @param z_pred Treatment variable for predictions (optional except if x_pre is not empty)
#' @param pi_pred propensity score for prediction
#' @param save_tree_directory directory where the trees have been saved
#' @param n_cores An optional integer of the number of cores to run your MCMC chains on
#' @param skip_con Whether to skip prediction of mu, and only report tau
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y               = y,
#'               z               = z,
#'               x_control       = x,
#'               x_moderate      = x,
#'               pihat           = pihat,
#'               nburn           = n_burn,
#'               nsim            = n_sim,
#'               n_chains        = 2,
#'               update_interval = 1,
#'               save_tree_directory = './trees')
#'
#' # Predict using new data
#'
#' x_pred = matrix(rnorm(n*p), nrow=n)
#'
#' pred_out = predict(bcf_out=bcf_fit,
#'                    x_predict_control=x_pred,
#'                    x_predict_moderate=x_pred,
#'                    pi_pred=pihat,
#'                    z_pred=z,
#'                    save_tree_directory = './trees')
#'
#'}
#' @export
predict.bcf <- function(object,
                        x_predict_control=NULL,
                        x_predict_moderate,
                        pi_pred,
                        z_pred,
                        save_tree_directory=NULL,
                        n_cores=2,
                        skip_con=FALSE,
                        ...) {

    if (skip_con) {
        x_predict_control <- matrix(1, nrow(x_predict_moderate))
    }

    if (is.null(x_predict_control) & !skip_con) {
        stop('If you want to predict mu, you need to add values to x_pred_control')
    }

    if (any(is.na(x_predict_moderate)))      stop("Missing values in x_predict_moderate")
    if (any(is.na(x_predict_control)))       stop("Missing values in x_predict_control")
    if (any(is.na(z_pred)))                  stop("Missing values in z_pred")
    if (any(!is.finite(x_predict_moderate))) stop("Non-numeric values in x_pred_moderate")
    if (any(!is.finite(x_predict_control)))  stop("Non-numeric values in x_pred_control")
    if (any(!is.finite(pi_pred)))            stop("Non-numeric values in pi_pred")
    if (any(!(z_pred %in% c(0,1))))          stop('z_pred must be 0s and 1s only')

    if((is.null(x_predict_moderate) & !is.null(x_predict_control)) | (!is.null(x_predict_moderate) & is.null(x_predict_control))) {
        stop("If you want to predict, you need to add values to both x_pred_control and x_pred_moderate")
    }

    pi_pred = as.matrix(pi_pred)
    if(!.ident(length(z_pred),
               nrow(x_predict_moderate),
               nrow(x_predict_control),
               nrow(pi_pred))
    ) {
        stop("Data size mismatch. The following should all be equal:
            length(z_pred): ", length(z_pred), "\n",
             "nrow(x_pred_moderate): ", nrow(x_predict_moderate), "\n",
             "nrow(x_pred_control): ", nrow(x_predict_control), "\n",
             "nrow(pi_pred): ", nrow(pi_pred), "\n"
        )
    }

    cat("Initializing BCF Prediction\n")
    x_pm = matrix(x_predict_moderate, ncol=ncol(x_predict_moderate))
    x_pc = matrix(x_predict_control, ncol=ncol(x_predict_control))

    if(object$raw_chains[[1]]$include_pi=="both" | object$raw_chains[[1]]$include_pi=="control") {
        x_pc = cbind(x_predict_control, pi_pred)
    }
    if(object$raw_chains[[1]]$include_pi=="both" | object$raw_chains[[1]]$include_pi=="moderate") {
        x_pm = cbind(x_predict_moderate, pi_pred)
    }

    n_chains = length(object$raw_chains)

    if (is.null(save_tree_directory)) {
        if (is.null(object$raw_chains[[1]]$con_trees) | is.null(object$raw_chains[[1]]$mod_trees)) {
            stop('Must providee tree directory if trees not included in fit')
        }
        save_tree_directory <- tempdir()
        for (i in 1:n_chains) {
            writeLines(object$raw_chains[[i]]$con_trees, file.path(save_tree_directory, paste('con_trees', i, 'txt', sep='.')))
            writeLines(object$raw_chains[[i]]$mod_trees, file.path(save_tree_directory, paste('mod_trees', i, 'txt', sep='.')))
        }
    }

    files_not_found <- c()
    for (i in 1:n_chains) {
        cfnm <- file.path(save_tree_directory, paste('con_trees', i, 'txt', sep='.'))
        if (!skip_con & !file.exists(cfnm)) {
            files_not_found <- c(files_not_found, cfnm)
        }

        mfnm <- file.path(save_tree_directory, paste('mod_trees', i, 'txt', sep='.'))
        if (!file.exists(mfnm)) {
            files_not_found <- c(files_not_found, mfnm)
        }
    }

    if (length(files_not_found)>0) {
        stop(paste('Could not find expected files ', paste(files_not_found, collapse=', ')))
    }

    cat("Starting Prediction \n")

    templog = tempfile()
    do_type_config <- .get_do_type(n_cores,templog)
    `%doType%` <- do_type_config$doType

    chain_out <- foreach::foreach(iChain=1:n_chains) %doType% {

        tree_files = .get_chain_tree_files(save_tree_directory, iChain)

        cat("Starting to Predict Chain ", iChain, "\n")

        mods = TreeSamples$new()
        mods$load(tree_files$mod_trees)
        Tm = mods$predict(t(x_pm))

        if (!skip_con) {
            cons = TreeSamples$new()
            cons$load(tree_files$con_trees)
            Tc = cons$predict(t(x_pc))
        } else {
            Tc = matrix(0, nrow(Tm), ncol(Tm))
        }

        list(Tm = Tm,
             Tc = Tc)
    }

    all_yhat = c()
    all_mu   = c()
    all_tau  = c()

    chain_list=list()

    muy = object$raw_chains[[1]]$muy

    sdy = object$raw_chains[[1]]$sdy

    for (iChain in 1:n_chains){
        # Extract Chain Specific Information
        Tm = chain_out[[iChain]]$Tm
        Tc = chain_out[[iChain]]$Tc

        this_chain_bcf_out = object$raw_chains[[iChain]]

        b1 = this_chain_bcf_out$b1
        b0 = this_chain_bcf_out$b0
        mu_scale = this_chain_bcf_out$mu_scale

        # Calculate, tau, y, and mu
        mu  = muy + sdy*(Tc*mu_scale + Tm*b0)
        tau = sdy*(b1 - b0)*Tm
        yhat = mu + t(t(tau)*z_pred)

        # Package Output up
        all_yhat = rbind(all_yhat, yhat)
        all_mu   = rbind(all_mu,   mu)
        all_tau  = rbind(all_tau,  tau)

        if (!skip_con) {
            scalar_df <- data.frame("tau_bar"   = matrixStats::rowWeightedMeans(tau, w=NULL),
                                    "mu_bar"    = matrixStats::rowWeightedMeans(mu, w=NULL),
                                    "yhat_bar"  = matrixStats::rowWeightedMeans(yhat, w=NULL))
        } else {
            scalar_df <- data.frame("tau_bar"   = matrixStats::rowWeightedMeans(tau, w=NULL))
        }

        chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
    }

    .cleanup_after_par(do_type_config)

    if (!skip_con) {
        result <- list(tau = all_tau,
                       mu = all_mu,
                       yhat = all_yhat,
                       coda_chains = coda::as.mcmc.list(chain_list))
    } else {
        result <- list(tau = all_tau,
                       coda_chains = coda::as.mcmc.list(chain_list))
    }

    return(result)
}
