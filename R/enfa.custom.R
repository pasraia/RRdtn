enfa.custom<-function (s.dat, x, obs_col, geoID_col, time_col = NULL){
  x[, !grepl(paste(c(geoID_col, time_col), collapse = "|"),
             colnames(x))] <- as.data.frame(scale(x[, !grepl(paste(c(geoID_col,
                                                                     time_col), collapse = "|"), colnames(x))]))
  if (is.null(time_col)) {
    pres <- s.dat[, geoID_col]
    S <- x[x[, geoID_col] %in% pres, !grepl(geoID_col, colnames(x))]
  }
  if (!is.null(time_col)) {
    x_time <- split(x, x[, time_col])
    S_time <- lapply(names(x_time), function(u) {
      TT <- x_time[u][[1]]
      pres <- s.dat[s.dat[, time_col] == u, geoID_col]
      S_time <- TT[TT[, geoID_col] %in% pres, !grepl(paste(c(geoID_col,
                                                             time_col), collapse = "|"), colnames(TT))]
      S_time <- as.matrix(S_time)
      S_time
    })
    S <- do.call(rbind, S_time)
  }
  S <- as.matrix(S)
  nS <- nrow(S)
  Rg <- cov(x[, !grepl(paste(c(geoID_col, time_col),
                             collapse = "|"), colnames(x))])
  p <- s.dat[, obs_col]
  p.sum <- sum(p)
  mar <- apply(S, 2, function(x) sum(x * p))/p.sum
  Sm <- sweep(S, 2, mar)
  DpSm <- apply(Sm, 2, function(x) x * p)
  Rs <- crossprod(Sm, DpSm)/(p.sum - 1)
  x <- x[, !grepl(paste(c(geoID_col, time_col), collapse = "|"),
                  colnames(x))]
  cZ <- ncol(x)
  m <- sqrt(as.numeric(t(mar) %*% mar))
  if (max(Im(eigen(Rs)$values)) > 1e-05)
    stop("complex eigenvalues. Try removing correlated variables.")
  eigRs <- lapply(eigen(Rs), Re)
  keep <- (eigRs$values > 1e-09)
  Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*%
    t(eigRs$vectors[, keep])
  W <- Rs12 %*% Rg %*% Rs12
  z <- Rs12 %*% mar
  y <- z/sqrt(sum(z^2))
  H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*%
                                            t(y))
  sf <- eigen(H)$values[-cZ]
  s.p <- (t(mar) %*% Rg %*% mar)/(t(mar) %*% Rs %*% mar)
  s <- c(s.p, sf)
  spec <- sqrt(mean(s))
  s.p <- abs(s)/sum(abs(s))
  v <- Re(eigen(H)$vectors)
  U <- matrix(nrow = cZ, ncol = cZ)
  u <- as.matrix((Rs12 %*% v)[, 1:(cZ - 1)])
  norw <- sqrt(diag(t(u) %*% u))
  U[, -1] <- sweep(u, 2, norw, "/")
  U[, 1] <- mar
  nm <- c("Marg", paste0("Spec", (1:(cZ - 1))))
  colnames(U) <- names(s.p) <- names(s) <- nm
  rownames(U) <- names(mar) <- names(x)
  mod <- list(call = "enfa", mf = mar, marginality = m,
              sf = s, specialization = spec, sf.prop = s.p, co = U,
              cov = Rs, weights = s.dat)
  return(mod)
}
