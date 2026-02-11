################################
### mixcure.penal.polar.plcr ###
################################


mixcure.penal.polar.plcr <- function(formula, data, apct = 0.05,init, pl){

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
    iterlim = 100, hessian=T);

  hessmat <- maximizer0$hessian

  var.mat <- solve(hessmat)

  alpha.hat <- maximizer0$estimate[index.gamma];

  loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value


  loglik.mixture.profile <- function(p, survt, k = k, design.matrix,
    ka.hat, kb.hat, r.est, phi.est,
    index.cure.var = index.cure.v, index.surv.var = index.surv.v, pl) {
    t      <- survt[, 1]
    status <- survt[, 2]
    event  <- (status == 1L)
    cens   <- !event
    logt   <- log(t)

    ik <- k - length(index.cure.var)

    # gamma (assumes index.gamma exists in caller env)
    gamma <- p[index.gamma - 2]
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    # ---- theta & eps ----
    X <- as.matrix(design.matrix)

    adj_cure <- (ka.hat + r.est * cos(phi.est))
    adj_surv <- (kb.hat + r.est * sin(phi.est))

    # if (length(index.cure.var) < 3) {
    #   # scalar multiplication path
    #   lp_cure <- drop(X[, index.cure.var[-k], drop = FALSE] * as.matrix(p[index.cure.var[-length(index.cure.var)]])) -
    #     X[, k] * adj_cure
    #   theta <- plogis(lp_cure)
    #
    #   lp_surv <- drop(X[, index.cure.var[-k], drop = FALSE] * as.matrix(p[(index.surv.var[-length(index.cure.var)] - 1)])) +
    #     X[, k] * adj_surv
    #   logeps <- gamma * logt + lp_surv
    #   eps <- exp(logeps)
    # } else {
      lp_cure <- drop(X[, index.cure.var[-k], drop = FALSE] %*% as.matrix(p[index.cure.var[-length(index.cure.var)]])) -
        X[, k] * adj_cure
      theta <- plogis(lp_cure)

      lp_surv <- drop(X[, index.cure.var[-k], drop = FALSE] %*% as.matrix(p[(index.surv.var[-length(index.cure.var)] - 1)])) +
        X[, k] * adj_surv
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
  #  }

    # Stable rewrite: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    eta   <- emeps / D
    delta <- (1 - theta) * eta

    # IMPORTANT: keep your "est/PLCI" kap (as in this function)
    kap <- (1 - eta) * (1 - theta) * (theta + eta)

    # pic = exp(eps)*eps*eta^2 (rewrite safely)
    pic <- eps * emeps / (D * D)

    # ---- negative log-likelihood ----
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))

    if (!isTRUE(pl)) return(loglikelihood)

    # ====================
    # Fast Fisher blocks
    # ====================
    max.len <- max(length(index.cure.var), length(index.surv.var))

    # Mirror loop indexing: i,j assumed to index first max.len cols
    XA <- X[, 1:max.len, drop = FALSE]

    e <- as.numeric(event)
    c <- 1 - e

    # --- Block A (subset to index.cure.var)
    wA <- e * (theta * (1 - theta)) + c * kap
    A_full <- crossprod(XA, XA * wA)
    info.a <- A_full[index.cure.var, index.cure.var, drop = FALSE]

    # --- Block B: uses design.xt = [design.matrix, logt]
    Xt <- cbind(XA, logt)  # (n x (max.len+1))
    wB <- c * (eps * (1 - delta) * delta)
    B_full <- -crossprod(XA, Xt * wB)   # left matrix is design.matrix[,i]
    cols_B <- c(index.surv.var - max.len, index.gamma - max.len)
    info.b <- B_full[index.cure.var, cols_B, drop = FALSE]

    # --- Block D: uses same design.xt on both sides
    wd2 <- eps * delta - (eps^2) * delta + (eps^2) * (delta^2)
    wD <- e * eps + c * wd2
    D_full <- crossprod(Xt, Xt * wD)
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
  #########################################################
  l.null = loglik - 0.5 * qchisq(1-apct,df=2,ncp = 0,lower.tail=T)

  phi_value = seq(0, pi*2, length=100)

  ni = 1
  tol = 0.4
  r_phi <- matrix(0, nrow = 100, ncol = (length(index.cure.v)+1))

  #ncores=detectCores()-1
 # cl<-makeCluster(10) #change the 2 to your number of CPU cores
 # registerDoParallel(cl)

  #r_phi <-rep(0, (length(index.cure.v)+1))

  #  for(i.phi in 1:100){
  for (i.phi in 1:100){
    require(R.utils)

      for (k in index.cure.v[-1]) {
        ik = k+length(index.cure.v);
        max.est = maximizer0$estimate
        n = ni + 1
        r.est=0.25
        param.est.up = max.est
        sign.delta =1

      converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; delta.up = 0
      while (!converge & iter1 <= 25 & !EXIT1 & !is.nan(delta.up)) {

        maximizer.temp <-  nlm(
          f = loglik.mixture, p = max.est, design.matrix=design.matrix,
          survt=survt, pl=pl,
          iterlim = 60, hessian=TRUE)
        score.temp = maximizer.temp$gradient
        hessian.temp = maximizer.temp$hessian
        if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06

        l0.b.up <- -maximizer.temp$minimum
        inv.hessian.temp <- solve(hessian.temp)


        lambda <- (2*(l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)

        # maximizing loop for unpenalized estimates;
        #if (pl == F) {
        inside <- FALSE; iter2 <- 1;
        while (!inside & iter2 <= 100 & !is.nan(delta.up)) {

          # add increment to stepwise parameter value;
          param.est.temp.up <- param.est.up[-c(k,ik)]
          r.est <- r.est + delta.up *sign.delta
          # if (k==2) {param.est.temp.up[1] <- -1;param.est.temp.up[9] <- 0.1}

          #compute loglikelihood function using updated parameter values;


          maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.up,
                                  survt=survt, k = k, design.matrix=design.matrix,
                                  ka.hat=max.est[k], kb.hat=max.est[ik], r.est=r.est, phi.est=phi_value[i.phi],
                                  pl = pl, iterlim = 100, hessian=F)

          l.temp.up = -maximizer.temp1$minimum

          #if (!is.nan(l.temp.up))

          #compare to see if updated l is still
          inside <- (l.temp.up < (l.null - 0.05)) #l.null - 0.05 for all others, 0.2 for k=3 of high rate H0
          alevel.up <- pchisq(2*(loglik - l.temp.up), df=2, ncp=0,lower.tail = T)
          #print(c(delta.up, alevel.up, n,l.temp.up,k,iter1,iter2, i.phi))
          if (!inside) {delta.up <- delta.up/((n+4)/n); sign.delta=1; iter2 <- iter2 + 1}  #(n+0.1)/n for low rate H0;
          # if (is.nan(delta.up)) {param.est.temp.up[k] <- NA}
        } #for iter2

        #}
        #Using converged increment for parameter to get corresponding score and variance expressions;

        param.est.up <- insert(maximizer.temp1$estimate, ats=k, values=r.est[,1]*cos(phi_value[i.phi]))
        param.est.up <- insert(param.est.up, ats=ik, values=r.est[,1]*sin(phi_value[i.phi]))

        l0.b.up = l.temp.up

        diff.up = l0.b.up - l.null
        converge <- (abs(diff.up) <= tol)
        if ((!converge| is.nan(l0.b.up)) & !is.nan(delta.up)) {iter1 <- iter1 + 1; n = n + 1; sign.delta=-1} else {EXIT1 = T}
        if (is.nan(delta.up)==T) {param.est.up[k] <- NA}
      } #for iter1

      r_phi[i.phi, k] <- r.est
     # print(c(i.phi, k))
      }
      r_phi[i.phi, 1] <- i.phi
      r_phi[i.phi, (length(index.cure.v)+1)] <- phi_value[i.phi]
  #    print(r_phi)
  #    return(r_phi)
  }

 # r_phi[, (length(index.cure.v)+1)] <- phi_value
  colnames(r_phi) <- c(colnames(design.matrix), "phi")

  z.score <- maximizer0$estimate / sqrt(diag(var.mat));

  coef.table.cure <- cbind(
    'coef'        = maximizer0$estimate[index.cure.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.cure.v]),
    'se(coef)'    = sqrt(diag(var.mat)[index.cure.v]),
    'z'           = z.score[index.cure.v],
    'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.cure.v]))),
    'LCI.95%' = maximizer0$estimate[index.cure.v] - 1.96 * sqrt(diag(var.mat)[index.cure.v]),
    'UCI.95%' = maximizer0$estimate[index.cure.v] + 1.96 * sqrt(diag(var.mat)[index.cure.v])
  );
  rownames(coef.table.cure) <- colnames(design.matrix);

  coef.table.surv <- cbind(
    'coef'        = maximizer0$estimate[index.surv.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.surv.v]),
    'se(coef)'    = sqrt(diag(var.mat)[index.surv.v]),
    'z'           = z.score[index.surv.v],
    'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.surv.v]))),
    'LCI.95%' = maximizer0$estimate[index.surv.v] - 1.96 * sqrt(diag(var.mat)[index.surv.v]),
    'UCI.95%' = maximizer0$estimate[index.surv.v] + 1.96 * sqrt(diag(var.mat)[index.surv.v])
  );
  rownames(coef.table.surv) <- colnames(design.matrix);

  coef.table.alpha <- cbind(
    'coef'     = alpha.hat,
    'se(coef)' = sqrt(diag(var.mat)[index.gamma]),
    'z'        = z.score[index.gamma],
    'Pr(>|z|)' = 2 * (1 - pnorm(abs(z.score[index.gamma]))),
    'LCI.95%'  = maximizer0$estimate[index.gamma] - 1.96 * sqrt(diag(var.mat)[index.gamma]),
    'UCI.95%'  = maximizer0$estimate[index.gamma] + 1.96 * sqrt(diag(var.mat)[index.gamma]),
    'loglik' = -maximizer0$minimum
  );
  rownames(coef.table.alpha) <- 'alpha';

  out <- list(
    coefficients = list(
      cure = coef.table.cure,
      surv = coef.table.surv,
      alpha = coef.table.alpha
      #   run.time
    ),
    plcr = as.data.frame(r_phi),
    cov = var.mat
  );

  return(out)
}
