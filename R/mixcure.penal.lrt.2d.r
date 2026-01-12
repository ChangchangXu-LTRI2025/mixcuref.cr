#######################################
### Function for mixcure model under ##
#### penalized loglikelihoods or ML  ##
#######################################################################
#### PROFILE LIKELIHOOD CONFIDENCE REGION FOR BIVARIATE PARAMETERS ####
#### REGION ESTIMATION FOR CURE PART OF MC                         ####
#######################################################################

########### CREATED ON 2021-10-26;
########### LAST MODIFIED 2021-10-26:

mixcure.penal.lrt.2d <- function(formula, data, init, pl, est.tru, iterlim = 200) {
require(splines)
require(survival)
require(abind)

  #########################################################################################

  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];

  design.matrix <- model.matrix(formula, data = design.matrix);

  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix);
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;

  samp.s <- nrow(design.matrix)

  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl, PLCI=F) {

    ####  parameter and variable dep parameters;
    t <- survt[, 1];  status <- survt[, 2];  event <- (status == 1L);  cens  <- !event;  logt  <- log(t)
    # Sub-matrices for the two linear predictors
    Xc <- design.matrix[, index.cure.var, drop = FALSE]
    k.vec <- ncol(Xc)                       # number of covariates in each part
    Xs <- design.matrix[, (index.surv.var - k.vec), drop = FALSE]
    index.gamma <- 2 * k.vec + 1            # matches parameter layout
    gamma <- p[index.gamma]

    # Guard (optional but helps nlm avoid NaNs)
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    # Cure part: theta
    lp_cure <- drop(Xc %*% p[index.cure.var])
    theta <- plogis(lp_cure)
    # Survival part: eps = t^gamma * exp(Xs beta) = exp(gamma*logt + Xs beta)
    lp_surv <- drop(Xs %*% p[index.surv.var])
    logeps <- gamma * logt + lp_surv
    eps <- exp(logeps)
    # Common denominator: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    # Negative log-likelihood
    # event term: log(1-theta) + log(gamma) - log(t) + log(eps) - eps
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))

    if (!pl) return(loglikelihood)

    # Quantities for penalty
    eta   <- emeps / D
    delta <- (1 - theta) * eta
    kap   <- (1 - eta) * (1 - theta) * (theta + eta)
    pi    <- eps * emeps / (D * D)      # equals exp(eps)*eps*eta^2 but overflow-safe
    # Build Xt once: survival covariates + log(t) for gamma column
    Xt <- cbind(Xs, logt)
    # Convert event/cens to 0/1 for weighting without subsetting
    e <- as.numeric(event)
    c <- 1 - e

    # --- Block A (k x k): for cure by cure parameters
    w1 <- theta * (1 - theta)
    wA <- e * w1 + c * kap
    info.a <- crossprod(Xc, Xc * wA)

    # --- Block B (k x (k+1)): for cure by surv parameters, only censored contribute
    wB <- c * (w1 * pi)
    info.b <- -crossprod(Xc, Xt * wB)

    # --- Block D ((k+1) x (k+1)): for surv by surv parameters
    wd2 <- eps * delta - eps^2 * delta + eps^2 * delta^2
    wD <- e * eps + c * wd2
    info.d <- crossprod(Xt, Xt * wD)

    # Add your extra gamma-gamma term
    info.d[k.vec + 1, k.vec + 1] <- info.d[k.vec + 1, k.vec + 1] + sum(event) / (gamma * gamma)

    # Schur complement without explicit inverse:
    # S = A - B %*% D^{-1} %*% t(B)
    tmp <- solve(info.d, t(info.b))     # (k+1) x k
    S <- info.a - info.b %*% tmp

    # log|det| in a stable way
    detS <- determinant(S, logarithm = TRUE)
    detD <- determinant(info.d, logarithm = TRUE)

    if (!PLCI) {if (detS$sign * detD$sign <= 0) return(Inf)}
    logdet <- as.numeric(detS$modulus) + as.numeric(detD$modulus)
    loglikelihood - 0.5 * logdet
  }

  ######END of loglik.mixture####################################


  # Parameter estimation under Ha (non-restricted likelihood)
  # maximize penalized or unpenalized loglikelihood by nlm;
  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix,
    pl = pl,
    iterlim = iterlim, hessian=TRUE);

  loglik0 <- -maximizer0$minimum

  loglik.mixture.profile <- function(p, survt,k.cur, k.sur, param.est.cur, param.est.sur,
    design.matrix1 = design.matrix, design.matrix0 = design.matrix,
    index.cure.var = index.cure.v, index.surv.var = index.surv.v, pl) {
    t      <- survt[, 1]
    status <- survt[, 2]
    event  <- (status == 1L)
    cens   <- !event
    logt   <- log(t)

    design.mtx.comb <- cbind(design.matrix0, design.matrix1)

    # gamma (assumes index.gamma exists in caller env)
    gamma <- p[index.gamma - 2]
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    X0 <- as.matrix(design.matrix0)
    X1 <- as.matrix(design.matrix1)

    # ---- theta & eps (same structure; stabilise plogis and avoid exp(eps) overflow later) ----
    lp_cure <- drop(X1[, index.cure.var[-k.cur], drop = FALSE] %*%
                      as.matrix(p[index.cure.var[-length(index.cure.var)]])) -
      design.mtx.comb[, k.cur] * param.est.cur
    theta <- plogis(lp_cure)

    lp_surv <- drop(design.mtx.comb[, index.surv.var[-k.cur], drop = FALSE] %*%
                      as.matrix(p[index.surv.var[-length(index.surv.var)] - 1])) +
      design.mtx.comb[, k.sur] * param.est.sur
    logeps <- gamma * logt + lp_surv
    eps <- exp(logeps)

    # Stable rewrite: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    eta   <- emeps / D
    delta <- (1 - theta) * eta

    # IMPORTANT: keep your kap for est/PLCI
    kap <- (1 - eta) * (1 - theta) * (theta + eta)

    # pi = exp(eps)*eps*eta^2 (rewrite safely)
    pi <- eps * emeps / (D * D)

    # ---- negative log-likelihood ----
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))
    if (!isTRUE(pl)) return(loglikelihood)

    # =========================================================================================
    # Fast Fisher blocks
    # =========================================================================================
    max.len <- max(length(index.cure.var), length(index.surv.var))

    # Mirror the loop indexing behaviour: i,j treated as 1..max.len columns
    X0A <- X0[, 1:max.len, drop = FALSE]
    X1A <- X1[, 1:max.len, drop = FALSE]

    e <- as.numeric(event)
    c <- 1 - e

    # --- Block A: event weight theta(1-theta), cens weight kap; subset to index.cure.var
    wA <- e * (theta * (1 - theta)) + c * kap
    A_full <- crossprod(X0A, X0A * wA)
    info.a <- A_full[index.cure.var, index.cure.var, drop = FALSE]

    # --- Block B: -sum( design.matrix1[,i] * design.xt0[,j] * eps*(1-delta)*delta ) over cens
    Xt0 <- cbind(X0A, logt)  # (n x (max.len+1))
    wB <- c * (eps * (1 - delta) * delta)
    B_full <- -crossprod(X1A, Xt0 * wB)  # (max.len x (max.len+1))

    cols_B <- c(index.surv.var - max.len, index.gamma - max.len)
    info.b <- B_full[index.cure.var, cols_B, drop = FALSE]

    # --- Block D: built from design.xt1 = [design.matrix1, logt], then subset
    Xt1 <- cbind(X1A, logt)  # (n x (max.len+1))
    wd2 <- eps * delta - (eps^2) * delta + (eps^2) * (delta^2)
    wD <- e * eps + c * wd2

    D_full <- crossprod(Xt1, Xt1 * wD)
    D_full[max.len + 1, max.len + 1] <- D_full[max.len + 1, max.len + 1] + sum(event) / (gamma * gamma)

    rowscols_D <- c(index.surv.var - max.len, index.gamma - max.len)
    info.d <- D_full[rowscols_D, rowscols_D, drop = FALSE]

    # Schur complement without explicit inverse
    tmp <- solve(info.d, t(info.b))
    info.set0 <- info.a - info.b %*% tmp

    # Stable log(det(info.set0)*det(info.d))
    detS <- determinant(info.set0, logarithm = TRUE)
    detD <- determinant(info.d, logarithm = TRUE)
    if (detS$sign * detD$sign <= 0) return(Inf)

    logdet <- as.numeric(detS$modulus) + as.numeric(detD$modulus)
    loglikelihood - 0.5 * logdet
  }

  #################################################################
  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for parameters of a single variable;

  dim.v <- ncol(design.matrix)

    ll.compl3 <- matrix(0, nrow=dim.v,ncol = 1)
    ll.pval <- matrix(0, nrow=dim.v,ncol = 1)

    for (k in index.cure.v[-1]) {

      ik = k + length(index.cure.v)
      maximizer <- nlm(
        f = loglik.mixture.profile, p = init[-c(k,ik)], k.cur = k, k.sur=ik,
        param.est.cur=est.tru[k], param.est.sur=est.tru[ik],
        survt = survt, design.matrix0 = design.matrix,
        design.matrix1=design.matrix,
        pl=pl, iterlim = iterlim, hessian=F);

      loglik.part = -maximizer$minimum;
      ll.compl3[k,1]<- loglik.part
      ll.pval[k,1] <- pchisq(2*abs(loglik0 - loglik.part), df=2, lower.tail = F)
    }

  ll.table <- cbind.data.frame(llr=ll.compl3,pval=ll.pval);
  rownames(ll.table) <- colnames(design.matrix);


#run.time = proc.time() - init.time


#######################################
## Output tables from either method; ##
#######################################


out <- list(ll.table);
class(out) <- c('mixcure', 'list');

return(out);

}

