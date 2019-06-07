
library('rmutil')


maxi <- function(v) {
  l <- 1
  for(i in 2:length(v)) {
    if(v[l] < v[i]) {
      l <- i
    }
  }
  return(l)
}

variance <- function(x, m) {
  return((1/(length(x) - 1))*sum((x - m)^2))
}



pub_interval <- function(db, a) {
  m <- mean(db)
  radius <- sqrt(variance(db, m)/length(db))*qnorm(1 - a/2)
  return(c(m - radius, m + radius))
}

pub_range <- function(db) {
  return(c(min(db), max(db)))
}


pub_histogram_learner <- function(db, bins) {
  #outputs a normalized histogram of db separated by the intervals in bins
  db <- sort(db)
  probs <- rep(0, length(bins) - 1)
  db_i <- 1
  while(db[db_i] < bins[1]) {
    db_i <- db_i + 1
    if(db_i > length(db)) {
      return(probs/sum(probs))
    }
  }
  for(i in 1:length(probs)) {
    while(db[db_i] < bins[i + 1]) {
      probs[i] <- probs[i] + 1
      db_i <- db_i + 1
      if(db_i > length(db)) {
        return(probs/sum(probs))
      }
    }
  }
  return(probs/sum(probs))
}

priv_std <- function(db, a, e, d, stdmin, stdmax) {
  bins_base <- floor(log2(stdmin) - 2)
  bins <- 2^(bins_base:ceiling(log2(stdmax) + 2))
  y <- 1:floor(length(db)/2)
  for(i in 1:length(y)) {
    y[i] <- abs(db[2*i] - db[2*i - 1])
  }

  l <- maxi(pub_histogram_learner(y, bins))
  return(2^((l + bins_base - 1) + 2))
}

priv_mean <- function(db, e, d, std, r) {
  rn <- ceiling(r/std)
  bins_base <- -rn
  bins <- ((bins_base:(rn + 1)) - .5)*std

  l <- maxi(pub_histogram_learner(db, bins))
  return((l + bins_base - 1)*std)
}

priv_range <- function(db, a1, a2, e1, e2, d1, d2, stdmin, stdmax, r) {
 priv_std <- priv_std(db, a1, e1, d1, stdmin, stdmax)
 priv_mean <- priv_mean(db, e2, d2, priv_std, r)

 radius <- 4*priv_std*sqrt(log(length(db)/a2))
 return(c(priv_mean - radius, priv_mean + radius))
}


priv_vadhan <- function(db, a0, a1, a2, a3, e1, e2, e3, d, stdmin, stdmax, r) {
  n <- length(db)
  xrange <- priv_range(db, a3/2, a3/2, e3/2, e3/2, d/2, d/2, stdmin, stdmax, r)
  xmin <- xrange[1]
  xmax <- xrange[2]
  xdist <- xmax - xmin

  for(i in 1:n) {
    #clamp x
    x <- db[i]
    if(x < xmin) {
      x <- xmin
    } else if(x > xmax) {
      x <- xmax
    }
    db[i] <- x
  }

  mean_var <- xdist/(e1*n)
  priv_mean <- mean(db) + rlaplace(1, 0, mean_var)

  var_var <- xdist^2/(e2*(n - 1))
  #priv_var <- public_var + extra_var + lap_noise
  priv_var <- variance(db, priv_mean) + var_var*log(1/a2) + rlaplace(1, 0, var_var)

  priv_radius <- sqrt(priv_var/n)*qt(1 - a0/2, n - 1) + mean_var*log(1/a1)

  return(c(priv_mean - priv_radius, priv_mean + priv_radius))
}


