require(optimbase)

# Bart’s density ----------------------------------------------------------

colos <- c(rgb(32/255, 74/255, 135/255, 0.7),
           rgb(204/255, 0, 0, 0.7),
           rgb(200/255, 141/255, 0, 0.7),
           rgb(78/255, 153/255, 6/255, 0.7),
           rgb(32/255, 74/255, 135/255, 0.3),
           rgb(204/255, 0, 0, 0.3),
           rgb(200/255, 141/255, 0, 0.3),
           rgb(78/255, 153/255, 6/255, 0.3))

# Package and function
suppressMessages(require(mixtools, quietly = T))


n <- 500
set.seed(1314)
XX <- rnormmix(n, 
               lambda = c(0.5, rep(0.1,5)), 
               mu     = c(0, ((0:4)/2)-1), 
               sigma  = c(1, rep(0.1,5)) 
)

# Make an histogram of the data
hist(XX, prob = T, col = gray(.8), border = NA, xlab = "x",
     main = paste("Data from Bart's density",sep=""),
     sub = paste("n = ", n, sep = ""),
     breaks = 50)
# Show the data points
rug(XX, col = rgb(0,0,0,.5))

# Plot the true density
true.den = function(x) 0.5*dnorm(x, 0, 1) + 
  0.1*dnorm(x,-1.0, 0.1) + 0.1*dnorm(x, -0.5, 0.1) +
  0.1*dnorm(x, 0.0, 0.1) + 0.1*dnorm(x,  0.5, 0.1) +
  0.1*dnorm(x, 1.0, 0.1)
curve(true.den, col = rgb(1,0,0,0.4), lwd = 3, n = n , add = TRUE)

# Kernel density estimate
lines(density(XX),            col = colos[3], lwd = 3)   # Oversmoothing
lines(density(XX, bw = .08),  col = colos[4], lwd = 3)   # Just Right
lines(density(XX, bw = .008), col = colos[5], lwd = 3)   # Undersmoothing

# Add a legend
legend("topright", c("True","Over", "Just right", "Under"), lwd = 5,
       col = c(rgb(1,0,0,0.4), colos[3], colos[4],colos[5]), cex = 0.8, bty = "n")

# Question 1: Extended version of handmade.em ---------------------------------------------------------

handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T)
{

  like <- p[1]*dnorm(y, mu[1], sigma[1])
  for (i in 2:length(p)){
    like <- like + p[i]*dnorm(y, mu[i], sigma[i])
  }
  
  deviance <- -2*sum(log(like))
  n <- length(p)
  l <- (2+(n*3))
  res      <- matrix(NA , n_iter + 1, l)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  for (iter in 1:n_iter) {
    # E step
    d <- c(p[1]*dnorm(y, mu[1], sigma[1]))
    for (i in 2: length(p)){
      d <- c(d, p[i]*dnorm(y, mu[i], sigma[i]))
    }
    d <- matrix(d, nrow = length(p), byrow = TRUE)
    
    r <- d[1,]/colSums(d)
    for (i in 2: length(p)){
      r <- array(c(r, d[i,]/colSums(d)))
    }
    r <- matrix(r, nrow = length(p), byrow = TRUE)
    
    # M step
    
    for (i in 1: length(p)){
      p[i] <- mean(r[i,])
      mu[i] <- sum(r[i,]*y)/sum(r[i,])
      sigma[i] <- sqrt(sum(r[i,]*(y^2))/sum(r[i,]) - (mu[i])^2)
    }
    
    
    # -2 x log-likelihood (a.k.a. deviance)
    like <- p[1]*dnorm(y, mu[1], sigma[1])
    for (i in 2:length(p)){
      like <- like + p[i]*dnorm(y, mu[i], sigma[i])
    }
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
  }
  res <- data.frame(res)
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}

# Question 2: Two different sample sizes ----------------------------------


# Question 3 --------------------------------------------------------------


# 3.a AIC ------------------------------------------------------------------

aic <- function(kmax, xx, plot_  = FALSE){
  n <- length(xx)
  k <- kmax
  aic_scores <- c()
  check_aic <- 99999999
  a=as.numeric(Sys.time())
  for (i in 2:kmax){
    set.seed(a)
    p <- ones(2,10)[1,]
    mu <- c(0, ((0:8)/2)-1)
    sigma <- runif(10, 0.0, 2.5)
    hem_fit <- handmade.em(xx,
                           p = p[1:i],
                           mu = mu[1:i],
                           sigma = sigma[1:i],
                           n_iter = 1000)
    
    MLE <- hem_fit$deviance/n
    loglikelihood <- -n/2*(log(2*pi) + 1 +log(MLE))
    q <- 3*i
    AIC <- -2*loglikelihood + 2*q
    aic_scores <- c(aic_scores, AIC)
    if (AIC < check_aic){
      check_aic <- AIC
      k <- i
    }
  }
  if (plot_ == TRUE){plot(c(2:kmax), aic_scores, type = "l",
       ylab = 'AIC Scores ', xlab = 'k-values', col = 'blue')
  text(c(2:kmax), aic_scores, round(aic_scores, 2), cex=0.6)
  sprintf("The optimal k value is found as: %i", k)}
  return(k)
}

aic(9, XX, plot_ = TRUE)


# 3. b BIC ----------------------------------------------------------------


bic <- function(kmax, xx, plot_ = FALSE){
  n <- length(xx)
  k <- kmax
  bic_scores <- c()
  check_bic <- 99999999
  a=as.numeric(Sys.time())
  for (i in 2:kmax){
    set.seed(a)
    p <- ones(2,10)[1,]
    mu <- c(0, ((0:8)/2)-1)
    sigma <- runif(10, 0.0, 2.5)
    hem_fit <- handmade.em(xx,
                           p = p[1:i],
                           mu = mu[1:i],
                           sigma = sigma[1:i],
                           n_iter = 1000)
    
    MLE <- hem_fit$deviance/n
    loglikelihood <- -n/2*(log(2*pi) + 1 +log(MLE))
    q <- 3*i
    BIC <- -2*loglikelihood + log(n)*q
    bic_scores <- c(bic_scores, BIC)
    if (BIC < check_bic){
      check_bic <- BIC
      k <- i
    }
  }
  
  if(plot_ == TRUE) {plot(c(2:kmax), bic_scores, type = "l",
       ylab = 'BIC Scores ', xlab = 'k-values', col = 'red')
  text(c(2:kmax), bic_scores, round(bic_scores, 2), cex=0.6)
  sprintf("The optimal k value is found as: %i", k)}
  
  return(k)
}

bic(9, XX, plot_ = TRUE)

# 3. c/d/e Sample splitting --------------------------------------------------------

sample_splitting <- function(size, xx, kmax, plot_ = FALSE){
  ## size% of the sample size
  size <- size/100
  smp_size <- floor(size * length(xx))
  
  ## set the seed to make your partition reproducible
  train_ind <- sample(seq_len(length(xx)), size = smp_size)
  train <- xx[train_ind]
  test <- xx[-train_ind]
  
  check <- 99999999
  k <- 0
  deviances <- c()
  a=as.numeric(Sys.time())
  for (i in 2:kmax){
  set.seed(a)
  p <- ones(2,10)[1,]
  mu <- c(0, ((0:8)/2)-1)
  sigma <- runif(10, 0.000001, (max(XX)- min(XX))/2)
  #now that we have train and test, start by finding the best parameters based on train
  hem_fit <- handmade.em(train,
                         p = p[1:i],
                         mu = mu[1:i],
                         sigma = sigma[1:i],
                         n_iter = 1000)
  
  par <- hem_fit$parameters
  #finally test those parameters on test
  p <- par[1:i]
  mu <- par[i+1:i]
  sigma <- par[2*i+1:i]
  
  like <- p[1]*dnorm(test, mu[1], sigma[1])
  for (j in 2:length(p)){
    like <- like + p[j]*dnorm(test, mu[j], sigma[j])
  }
  deviance <- -2*sum(log(like))
  deviances <- c(deviances, deviance)
  if (deviance < check){
    check <- deviance
    k <- i
  }
  }


  
  if(plot_ == TRUE) {plot(c(2:kmax), deviances, type = "l",
       ylab = 'Deviance ', xlab = 'k-values', col = 'green')
  text(c(2:kmax), deviances, round(deviances, 2), cex=0.6)
  sprintf("The optimal k value is found as: %i", k)}
  
  return(k)
  
  }

sample_splitting(30, XX, 8)


# 3. f/g n-fold Cross Validation ------------------------------------------

cross_validation <- function(n_fold, xx, kmax, plot_ = FALSE){
  
  folds <- createFolds(xx, k = n_fold, list = TRUE, returnTrain = FALSE)
  deviances <- c()
  a=as.numeric(Sys.time())
  for (k in 2:kmax){
    
    set.seed(a)
    p <- ones(2,10)[1,]
    mu <- c(0, ((0:8)/2)-1)
    sigma <- runif(10, 0.0, 2.5)
    
    avg_dev <- 0
    for (f in 1:n_fold){
      test <- xx[folds[[f]]]
      train <- xx[-folds[[f]]]
      
      
      #now that we have train and test, start by finding the best parameters based on train
      hem_fit <- handmade.em(train,
                             p = p[1:k],
                             mu = mu[1:k],
                             sigma = sigma[1:k],
                             n_iter = 1000)
      par <- hem_fit$parameters
      #finally test those parameters on test
      p <- par[1:k]
      mu <- par[k+1:k]
      sigma <- par[2*k+1:k]
      
      like <- p[1]*dnorm(test, mu[1], sigma[1])
      for (j in 2:length(p)){
        like <- like + p[j]*dnorm(test, mu[j], sigma[j])
      }
      deviance <- -2*sum(log(like))
      avg_dev <- avg_dev + deviance
      
      
    }
    avg_dev <- avg_dev/n_fold
    deviances <- c(deviances, avg_dev)
  }
  if (plot_ == TRUE) {plot(c(2:kmax), deviances, type = "l",
       ylab = 'Deviance ', xlab = 'k-values', col = 'green')
  text(c(2:kmax), deviances, round(deviances, 2), cex=0.6)
  sprintf("The optimal k value is found as: %i", k)}
  
  return(which.min(deviances) + 1)
  
}

cross_validation(5, XX, 8, plot_ = TRUE)


# Deciding on k* ----------------------------------------------------------

select_k <- function(M, n, kmax){
  store_k = c()
  for (i in 1:M){
    XX <- rnormmix(n, 
                   lambda = c(0.5, rep(0.1,5)), 
                   mu     = c(0, ((0:4)/2)-1), 
                   sigma  = c(1, rep(0.1,5)) 
    )
    options(show.error.messages = FALSE)
    try(store_k <- c(store_k, aic(kmax, XX)))
    try(store_k <- c(store_k, bic(kmax, XX)))
    try(store_k <- c(store_k, sample_splitting(50, XX, kmax)))
    try(store_k <- c(store_k, sample_splitting(70, XX, kmax)))
    try(store_k <- c(store_k, sample_splitting(30, XX, kmax)))
    try(store_k <- c(store_k, cross_validation(5, XX, kmax)))
    try(store_k <- c(store_k, cross_validation(10, XX, kmax)))
  }
  
  data=data.frame(value=store_k)
  p <- ggplot(data.frame(value=data), aes(x=value)) + 
    geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9)
  p
  print(store_k)
}

select_k(20, 500, 8)



# 3. h --------------------------------------------------------------------

