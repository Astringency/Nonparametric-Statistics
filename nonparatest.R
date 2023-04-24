one.sample.sign.test <- function(x, p, q0) {
  # one sample sign test
  s1 <- sum(x < q0)
  s2 <- sum(x > q0)
  n <- s1 + s2
  p1 <- pbinom(s1, n, p)
  p2 <- 1 - pbinom(s1 - 1, n, p)
  if (p1 > p2) {
    m1 <- "One tail test:H1:Q<q0"
  } else {
    m1 <- "One tail test: H1: Q>q0"
  }
  p.value <- min(p1, p2)
  m2 <- "Two tails test"
  p.value2 <- 2 * p.value
  if (q0 == median(c(t(x)))) {
    p.value <- 0.5
    p.value2 <- 1
  }

  list(
    Sign.test1 = m1, p.values.of.one.tail.test = p.value,
    p.value.of.two.tai.test = p.value2
  )
}

mci <- function(x, alpha = 0.05) {
  # Median confidence interval 1
  x <- sort(x)
  n <- length(x)
  b <- 0
  i <- 0
  while (b <= alpha / 2 & i <= floor(n / 2)) {
    b <- pbinom(i, n, 0.5)
    k1 <- i
    k2 <- n - i + 1
    a <- 2 * pbinom(k1 - 1, n, 0.5)
    i <- i + 1
  }
  z <- c(k1, k2, a, 1 - a)
  z2 <- "Entire range!"

  if (k1 >= 1) {
    out <- list(Confidence.level = 1 - a, CI = c(x[k1], x[k2]))
  } else {
    out <- list(Confidence.level = 1 - 2 * pbinom(0, n, 0.5), CI = z2)
  }
  out
}

mci2 <- function(x, alpha = 0) {
  # Median confidence interval 2
  x <- sort(x)
  n <- length(x)
  q <- 0.5
  m <- floor(n * q)
  s1 <- pbinom(0:m, n, q)
  s2 <- pbinom(m:(n - 1), n, q, low = F)
  ss <- c(s1, s2)
  nn <- length(ss)
  a <- NULL
  for (i in 0:m) {
    b1 <- ss[i + 1]
    b2 <- ss[nn - i]
    b <- b1 + b2
    d <- 1 - b
    if ((b) > 1) break
    a <- rbind(a, c(b, d, x[i + 1], x[n - i]))
  }
  colnames(a) <- c(
    "p-value", "Confidence level",
    "lower limit", "Higher limit"
  )
  if (a[1, 1] > alpha) {
    out <- a
  } else {
    for (i in 1:nrow(a)) {
      if (a[i, 1] > alpha) {
        out <- list(a[i - 1, 1], a[i - 1, 2], a[i - 1, 3:4])
        break
      }
    }
  }
  out
}

qci <- function(x, alpha = 0.05, q = 0.25) {
  # Quantile confidence interval
  x <- sort(x)
  n <- length(x)
  a <- alpha / 2
  r <- qbinom(a, n, q)
  s <- qbinom(1 - a, n, q)
  CL <- pbinom(s, n, q) - pbinom(r - 1, n, q)

  if (r == 0) {
    lo <- NA
  } else {
    lo <- x[r]
  }

  if (s == n) {
    up <- NA
  } else {
    up <- x[s + 1]
  }

  list(c("lower limit" = lo, "upper limit" = up),
    "1-alpha" = 1 - alpha, "true confidence level" = CL
  )
}


wilcoxon.test <- function(x, m, a = "greater") {
  # Wilcoxon Sign Test
  x_ <- c()
  for (i in x) {
    if (i == m) {
      x_ <- x_
    } else {
      x_ <- append(x_, i)
    }
  }
  x_ <- sort(x_)

  su <- cbind(x_, abs(x_ - m), rank(abs(x_ - m)), x_ < m)
  colnames(su) <- c("X_Real", "Abs_Diff", "Rank+", "Sign")

  W.neg <- sum(su[, 3] * (su[, 4] == 1))
  W.pos <- sum(su[, 3] * (su[, 4] == 0))

  list(
    "Summary Result" = su,
    "Statistics W-" = W.neg,
    "Statistics W+" = W.pos,
    "Test Result" = wilcox.test(x_ - m, alt = a)
    )
}

wilcoxon.estimate <- function(x, alpha = 0.05){
  # Wilcoxon Sign Rank Test Estimate
  n <- length(x)
  walsh <- NULL
  for (i in 1:n){
      for (j in i:length(x)){
            walsh <- c(walsh, (x[i] + x[j])/2)
      }
  }

  walsh <- sort(walsh)

  theta.hat <- median(walsh)

  N <- length(walsh)
  k <- qsignrank(alpha / 2, n)

  list(
    "Walsh Table" = walsh,
    "Statistics" = k,
    "Point Estimate" = theta.hat,
    "Lower Limit" = walsh[k[1] + 1],
    "Upper Limit" = walsh[N - k[1]],
    "Confidence Level" = 1 - alpha
  )
}

cox.stuart.test <- function(x, alpha=0.05) {
  # Cox-Stuart Test
  c <- NULL
  x1 <- NULL
  x2 <- NULL

  if (length(x) %% 2 == 0) {
    c <- length(x) / 2
  } else {
    c <- (length(x) + 1) / 2
  }

  if (length(x) %% 2 == 0) {
    x1 <- x[1:c]
    x2 <- x[(c + 1):length(x)]
  } else {
    x1 <- x[1:(c - 1)]
    x2 <- x[(c + 1):length(x)]
  }

  D <- x1 - x2

  s1 <- sum(sign(D) == 1)
  s2 <- sum(sign(D) == -1)
  n <- s1 + s2
  k <- min(s1, s2)
  p <- pbinom(k, n, 0.5)

  if (p < alpha) {
    h <- "Reject H0"
  } else {
    h <- "Fail to reject H0"
  }

  list(
    "D table" = D,
    "Statistics" = k,
    "p-value" = p,
    "p-value.two.tails" = 2 * p,
    "Hypothesis" = h
  )
}

run.test <- function(y, cut = 0) {
  # Run Test
  if (cut != 0) {
    x <- (y > cut) * 1
  } else {
    x <- y
  }

  N <- length(x)
  k <- 1
  for (i in 1:(N - 1)) {
    if (x[i] != x[i + 1]){
      k <- k + 1
    }
  }

  r <- k
  m <- sum(1 - x)
  n <- N - m

  P1 <- function(m, n, k) {
    2 * choose(m - 1, k - 1) / choose(m + n, n) * choose(n - 1, k - 1)
  }

  P2 <- function(m, n, k) {
    choose(m - 1, k - 1) * choose(n - 1, k) / choose(m + n, n) + choose(m - 1, k) * choose(n - 1, k - 1) / choose(m + n, n)
  }
  
  r2 <- floor(r/2)
  if (r2 == r/2) {
    pv <- 0
    for (i in 1:r2) {
      pv <- pv + P1(m, n, i)
    }
    for (i in 1:(r2 - 1)) {
      pv <- pv + P2(m, n, i)
    }
  } else {
    pv <- 0
    for (i in 1:r2) {
      pv <- pv + P1(m, n, i)
    }
    for (i in 1:r2) {
      pv <- pv + P2(m, n, i)
    }
  }
  if (r2 == r/2) {
    pv1 <- 1 - pv + P1(m, n, r2)
  } else {
    pv1 <- 1 - pv + P2(m, n, r2)
  }
  
  z <- (r - 2 * m * n / N - 1) / sqrt(2 * m * n * (2 * m * n - m - n) / (m + n)^2 / (m + n - 1))
  ap1 <- pnorm(z)
  ap2 <- 1 - ap1
  tpv <- min(pv, pv1) * 2
  
  list(
    m = m,
    n = n,
    N = N,
    R = r,
    Exact.pvalue1 = pv,
    Exact.pvalue2 = pv1,
    Aprox.pvalue1 = ap1,
    Aprox.pvalue2 = ap2,
    Exact.2sided.pvalue = tpv,
    Approx.2sided.pvalue = min(ap1, ap2) * 2
    )
}

Brown.Mood.test <- function(data, alt = "less") {
  k <- unique(data[, 2])
  m <- median(data[, 1])
  m1 <- NULL
  m2 <- NULL

  for (i in k) {
      m1 <- c(m1, sum(data[data[, 2] == i, 1] > m))
      m2 <- c(m2, sum(data[data[, 2] == i, 1] < m))
  }

  C <- rbind(m1, m2)

  a <- C[1, 1]
  b <- C[1, 2]
  m <- sum(C[, 1])
  n <- sum(C[, 2])
  t <- sum(C[1, ])
  N <- m + n

  if (a <= b) {
    z <- (a + 0.5 - m * t/N) / sqrt(m * n * t * (N-t)/N^3)
  } else {
    z <- (a - 0.5 - m * t/N) / sqrt(m * n * t * (N-t)/N^3)
  }

  chi <- (2 * a - m)^2 * (m + n) / m / n
  
  p.hyper.one.tail <- min(phyper(a, m, n, min(m, t)), 1 - phyper(a, m, n, min(m, t)))
  p.hyper.two.tail <- 2 * min(phyper(a, m, n, t), 1 - phyper(a, m, n, t))
  p.norm <- pnorm(z)
  p.chi <- 1 - pchisq(chi, 1)
  
  list(
    "Exact p-value one way" = p.hyper.one.tail, 
    "Exact p-value two way" = p.hyper.two.tail,
    "Approx p-value Normal" = p.norm, 
    "Approx p-value Chi2" = p.chi,
    "Fisher test" = fisher.test(C, alt = alt)
  )
}