##################################
### mixcure.penal.wald.cr      ###
### Creates:  estimates and se ###
### Wald 2d test p, and CR area ##
##################################


mixcure.penal.wald.cr.area <- function(formula, data, apct = 0.05,init, true.val, pl){

  require(splines)
  require(survival)
  require(abind)

  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];

  design.matrix <- model.matrix(formula, data = design.matrix);

  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix);
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;

  samp.s <- nrow(design.matrix)

  # p is a vector of values for parameter vector r (radius) for all covariates X and gamma (shape)

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

  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix,
    pl = pl,
    iterlim = 100, hessian=F);

  var.est <- maximizer0$estimate
  var.mat <- solve(maximizer0$hessmat)

  wald.2d.test <- rep(NA, length(ncol(design.matrix)))
  wald.2d.test.p <- rep(NA, length(ncol(design.matrix)))
  wald.2d.cr.area <- rep(NA, length(ncol(design.matrix)))
  corr.var <- rep(NA, length(ncol(design.matrix)))

  for (i in 1: ncol(design.matrix) ) {
    mat.index = c(i, (i+ncol(design.matrix)) )
    test.2d = (var.est-true.val)[mat.index]%*%solve(vcov.mat[mat.index, mat.index])%*%(var.est-true.val)[mat.index]
    wald.2d.test.p[i] = pchisq(test.2d, df=2,lower.tail = F)
    wald.2d.cr.area[i] = sqrt(diag(var.mat[i]))*sqrt(diag(var.mat[i+ncol(design.matrix)]))*pi
    wald.2d.test[i] <- test.2d
    corr.var[i] <- vcov.mat[mat.index, mat.index][2,1]
      }

  coef.table <- cbind(
    'coef.cure'   = var.est[index.cure.v],
    'coef.surv'   = exp(var.est[index.surv.v]),
    'se.cure'    = sqrt(diag(var.mat)[index.cure.v]),
    'se.surv'    = sqrt(diag(var.mat)[index.surv.v]),
    'corr'        = corr.var[index.surv.v],
    '2dchi'       = wald.2d.test[index.cure.v],
    'Pval'    = wald.2d.test.p[index.cure.v],
    'area_cr' = wald.2d.cr.area[index.cure.v]
   );
  rownames(coef.table) <- colnames(design.matrix);

  return(coef.table);
}
