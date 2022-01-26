#' @import mvtnorm
#' @import glmnet
#' @import progress

aux <- new.env()

#' @export
get_aux<-function()
{
  return(aux)
}




# bootstrap: 0=no (parametric); 1=ordinary; 2=Bayesian
#' @export
voi.glm<-function(reg_obj, n_sim=1000, bootstrap=0, lambdas=(1:99)/100)
{
  sample_size <- dim(reg_obj$data)[1]
  mu <- coefficients(reg_obj)
  sigma <- vcov(reg_obj)

  pi<-predict(reg_obj,type="response")
  aux$pi <- pi

  NB_test <-  NB_all <-  NB_max  <- rep(0, length(lambdas))

  for(i in 1:n_sim)
  {
    if(bootstrap>0)
    {
      ws <- bootstrap(sample_size, bootstrap-1)
      bs_data <- cbind(reg_obj$data, ws=ws)
      bs_reg <- glm(data=bs_data, formula = reg_obj$formula, family=reg_obj$family, weights = ws)
      p <- as.vector(predict(bs_reg, newdata=reg_obj$data, type="response"))
    }
    else
    {
      new_betas <- rmvnorm(1, mean=mu, sigma=sigma)
      p <- as.vector(1/(1+exp(-(new_betas %*% t(model.matrix(reg_obj))))))
    }

    if(i==1) plot(pi,p)

    for(j in 1:length(lambdas))
    {
      NB_test[j] <- NB_test[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))

      NB_all[j] <- NB_all[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)

      NB_max[j] <- NB_max[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
    }
    #cat('.')
  }

  voi <- (NB_max-pmax(0,NB_test,NB_all))/n_sim

  res <-cbind(lambda=lambdas, voi=voi, NB_all=NB_all/n_sim, NB_test=NB_test/n_sim, NB_max=NB_max/n_sim)

  return(res)
}





#' @export
#' @param reg_obj: any object that you can apply predict with new data to get predictions. The structure should represent the correct model
#' @param x: The model matrix of predictors
#' @param y: The vector of responses
#' @param pi: optional. Predictions from the current model. If not supplied, the predictions from reg_object will be used
voi.glmnet <- function(reg_obj, x, y, pi=NULL, n_sim=1000, lambdas=(1:99)/100, Bayesian_bootstrap=F, empirical=F)
{
  aux$coeffs <- t(as.matrix(coefficients(reg_obj)))
  aux$x <- x
  aux$y <- y
  aux$reg_obj <- reg_obj

  sample_size <- dim(x)[1]

  if(is.null(pi)) pi <- predict(reg_obj, type="response", newx=x)

  aux$pi <- pi

  NB_model <- rep(0, length(lambdas))
  NB_all <- NB_model
  NB_max <- NB_model

  NB_model_s2 <- rep(0, length(lambdas))
  NB_all_s2 <- NB_model_s2
  NB_max_s2 <- NB_model_s2
  NB_model_all_s2 <- NB_model_s2
  NB_model_max_s2 <- NB_model_s2
  NB_all_max_s2 <- NB_model_s2

  p_win_model <- p_win_all <- p_win_none <- rep(0, length(lambdas))

  dc_model <- NB_model
  dc_all <- NB_model

  optimism <- NB_model

  aux$bs_coeffs <- matrix(NA,nrow=n_sim, ncol=dim(coefficients(reg_obj)))

  colnames(aux$bs_coeffs) <- rownames(coefficients(reg_obj))

  for(j in 1:length(lambdas))
  {
    dc_model[j] <- dc_model[j] + mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
    dc_all[j] <- dc_all[j] + mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * 1)
  }

  pb <- progress::progress_bar$new(total=n_sim)

  for(i in 1:n_sim)
  {
    pb$tick()
    repeat
    {
      weights <- bootstrap(sample_size, Bayesian = Bayesian_bootstrap)
      issues <- F
      tryCatch(
      {
        tmp <- cv.glmnet(x=x, y=y, family="binomial",weights = as.vector(weights))
        bs_reg <- glmnet(x=x, y=y, family="binomial", lambda=tmp$lambda.min, weights=as.vector(weights))
      },warning=function(cond) {message("Warning occured! repeating with a new sample"); issues<- T; })
      if(!issues) break
    }

    aux$bs_coeffs[i,] <- t(as.matrix(coefficients(bs_reg)))

    p <- as.vector(predict(bs_reg, newx=x, type="response"))

    if(i==1) plot(pi,p)

    for(j in 1:length(lambdas))
    {
      tmp1 <- mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
      NB_model[j] <- NB_model[j] + tmp1
      NB_model_s2[j] <- NB_model_s2[j] + tmp1^2
      tmp2 <- mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)
      NB_all[j] <- NB_all[j] + tmp2
      NB_all_s2[j] <- NB_all_s2[j] + tmp2^2
      tmp3 <- mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
      NB_max[j] <- NB_max[j] + tmp3
      NB_max_s2[j] <- NB_max_s2[j] + tmp3^2
      NB_model_all_s2[j] <- NB_model_all_s2[j] + (tmp1-tmp2)^2
      NB_model_max_s2[j] <- NB_model_max_s2[j] + (tmp1-tmp3)^2
      NB_all_max_s2[j] <- NB_all_max_s2[j] + (tmp2-tmp3)^2

      winner <- which.max(c(tmp1,tmp2,0))
      if(winner==1) p_win_model[j] <- p_win_model[j]+1 else if(winner==2) p_win_all[j] <- p_win_all[j]+1 else p_win_none[j]<-p_win_none[j]+1

      dc_model_int <- sum((y + (1 - y) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]) * weights) / sum(weights)
      dc_model_ext <- mean((y + (1 - y) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
      optimism[j] <- optimism[j] + dc_model_int - dc_model_ext
    }
  }

  NB_model <- NB_model / n_sim
  NB_all <- NB_all / n_sim
  NB_max <- NB_max / n_sim
  NB_model_s2 <- NB_model_s2 / n_sim
  NB_all_s2 <- NB_all_s2 / n_sim
  NB_max_s2 <- NB_max_s2 / n_sim
  NB_model_all_s2 <- NB_model_all_s2 / n_sim
  NB_model_max_s2 <- NB_model_max_s2 / n_sim
  NB_all_max_s2 <- NB_all_max_s2 / n_sim
  p_win_model <- p_win_model/n_sim
  p_win_all <- p_win_all/n_sim
  p_win_none <- p_win_none/n_sim

  optimism <- optimism / n_sim

  voi <- (NB_max-pmax(0,NB_model,NB_all))

  res <-cbind(lambda=lambdas, voi=voi, NB_all=NB_all, NB_model=NB_model, NB_max=NB_max, p_win_model=p_win_model, p_win_all=p_win_all, p_win_none=p_win_none, dc_model=dc_model, dc_all=dc_all, optimism=optimism, NB_all_s2=NB_all_s2, NB_model_s2=NB_model_s2, NB_max_s2=NB_max_s2, NB_model_all_s2=NB_model_all_s2, NB_model_max_s2=NB_model_max_s2, NB_all_max_s2=NB_all_max_s2)

  return(res)
}





#' @export
voi.glmnet2 <- function(formula, data, pi, n_sim=1000, lambdas=(1:99)/100, Bayesian_bootstrap=F, empirical=F, weights=NULL)
{
  x <- model.matrix(formula,data)
  y <- data[,all.vars(formula)[1]]

  aux$coeffs <- colnames(x)
  aux$x <- x
  aux$y <- y
  aux$reg_obj <- NULL

  sample_size <- dim(x)[1]

  aux$pi <- pi

  NB_model <- rep(0, length(lambdas))
  NB_all <- NB_model
  NB_max <- NB_model

  NB_model_s2 <- rep(0, length(lambdas))
  NB_all_s2 <- NB_model_s2
  NB_max_s2 <- NB_model_s2
  NB_model_all_s2 <- NB_model_s2
  NB_model_max_s2 <- NB_model_s2
  NB_all_max_s2 <- NB_model_s2

  p_win_model <- p_win_all <- p_win_none <- rep(0, length(lambdas))

  dc_model <- NB_model
  dc_all <- NB_model

  optimism <- NB_model

  aux$bs_coeffs <- matrix(NA,nrow=n_sim, ncol=dim(x)[2])

  colnames(aux$bs_coeffs) <- colnames(x)

  for(j in 1:length(lambdas))
  {
    dc_model[j] <- dc_model[j] + mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
    dc_all[j] <- dc_all[j] + mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * 1)
  }

  pb <- progress::progress_bar$new(total=n_sim)

  for(i in 1:n_sim)
  {
    pb$tick()
    repeat
    {
      weights2 <- bootstrap(sample_size, Bayesian = Bayesian_bootstrap, weights=weights)
      issues <- F
      tryCatch(
        {
          tmp <- cv.glmnet(x=x, y=y, family="binomial",weights = as.vector(weights2))
          bs_reg <- glmnet(x=x, y=y, family="binomial", lambda=tmp$lambda.min, weights=as.vector(weights2))
        },warning=function(cond) {message("Warning occured! repeating with a new sample"); issues<- T; })
      if(!issues) break
    }

    aux$bs_coeffs[i,] <- t(as.matrix(coefficients(bs_reg)))[-2]

    p <- as.vector(predict(bs_reg, newx=x, type="response"))

    if(i==1) plot(pi,p)

    for(j in 1:length(lambdas))
    {
      tmp1 <- mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
      NB_model[j] <- NB_model[j] + tmp1
      NB_model_s2[j] <- NB_model_s2[j] + tmp1^2
      tmp2 <- mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)
      NB_all[j] <- NB_all[j] + tmp2
      NB_all_s2[j] <- NB_all_s2[j] + tmp2^2
      tmp3 <- mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
      NB_max[j] <- NB_max[j] + tmp3
      NB_max_s2[j] <- NB_max_s2[j] + tmp3^2
      NB_model_all_s2[j] <- NB_model_all_s2[j] + (tmp1-tmp2)^2
      NB_model_max_s2[j] <- NB_model_max_s2[j] + (tmp1-tmp3)^2
      NB_all_max_s2[j] <- NB_all_max_s2[j] + (tmp2-tmp3)^2

      winner <- which.max(c(tmp1,tmp2,0))
      if(winner==1) p_win_model[j] <- p_win_model[j]+1 else if(winner==2) p_win_all[j] <- p_win_all[j]+1 else p_win_none[j]<-p_win_none[j]+1

      dc_model_int <- sum((y + (1 - y) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]) * weights) / sum(weights)
      dc_model_ext <- mean((y + (1 - y) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
      optimism[j] <- optimism[j] + dc_model_int - dc_model_ext
    }
  }

  NB_model <- NB_model / n_sim
  NB_all <- NB_all / n_sim
  NB_max <- NB_max / n_sim
  NB_model_s2 <- NB_model_s2 / n_sim
  NB_all_s2 <- NB_all_s2 / n_sim
  NB_max_s2 <- NB_max_s2 / n_sim
  NB_model_all_s2 <- NB_model_all_s2 / n_sim
  NB_model_max_s2 <- NB_model_max_s2 / n_sim
  NB_all_max_s2 <- NB_all_max_s2 / n_sim
  p_win_model <- p_win_model/n_sim
  p_win_all <- p_win_all/n_sim
  p_win_none <- p_win_none/n_sim

  optimism <- optimism / n_sim

  voi <- (NB_max-pmax(0,NB_model,NB_all))
  voi_r <- (NB_max-pmax(0,NB_all))/(NB_model-pmax(0,NB_all))

  res <-cbind(lambda=lambdas, voi=voi, voi_r=voi_r, NB_all=NB_all, NB_model=NB_model, NB_max=NB_max, p_win_model=p_win_model, p_win_all=p_win_all, p_win_none=p_win_none, dc_model=dc_model, dc_all=dc_all, optimism=optimism, NB_all_s2=NB_all_s2, NB_model_s2=NB_model_s2, NB_max_s2=NB_max_s2, NB_model_all_s2=NB_model_all_s2, NB_model_max_s2=NB_model_max_s2, NB_all_max_s2=NB_all_max_s2)

  return(res)
}







bootstrap <- function (n, Bayesian=F, weights=NULL)
{
  if(Bayesian)
  {
    if(!is.null(weights)) stop("BAyesian bootstrap currently does not work with weighted samples.")
    u <- c(0,sort(runif(n-1)),1)
    return((u[-1] - u[-length(u)])*n)
  }
  else
  {
    if(is.null(weights)) weights <-rep(1/n,n)
    u <- rmultinom(1,n,weights)
    return(u)
  }
}








#' @export
process_results <- function(res, graphs=c("voi","summit","dc"),th=NULL)
{
  if(is.null(th)){
  out <- list()
  out$inb_current <- mean(pmax(0,res[,'NB_model'],res[,'NB_all'])-pmax(0,res[,'NB_all']))
  out$inb_perfect <- mean(res[,'NB_max']-pmax(0,res[,'NB_all']))
  out$voi_r <- out$inb_perfect/out$inb_current
  } else{
    index <- which(res[,'lambda'] %in% th)
    out <- list()
    out$inb_current <- pmax(0,res[,'NB_model'],res[,'NB_all'])-pmax(0,res[,'NB_all'])
    out$inb_perfect <- res[,'NB_max']-pmax(0,res[,'NB_all'])
    out$voi_r <- (out$inb_perfect/out$inb_current)[index]
  }

  if(!is.na(match("voi",graphs)))
  {
    plot(res[,'lambda'], res[,'NB_max']-pmax(0,res[,'NB_model'],res[,'NB_all']),type='l', lwd=2, col="red", xlab="Threshold", ylab="EVPI")
  }

  if(!is.na(match("summit",graphs)))
  {
    plot(res[,'lambda'],pmax(0,res[,'NB_model'],res[,'NB_all'])-pmax(0,res[,'NB_all']),type='l', xlab='Threshold', ylab='Incremental net benefit', lwd=2)
    lines(res[,'lambda'],res[,'NB_max']-pmax(0,res[,'NB_all']),type='l',col="red", lwd=2)
  }

  if(!is.na(match("dc",graphs)))
  {
    max_y <- max(res[,'NB_all'],res[,'NB_model'])
    plot(res[,'lambda'],res[,'NB_model'],type='l', xlab='Threshold', ylab='Net benefit', lwd=2, xlim=c(0,1), ylim=c(0,max_y), col="red")
    lines(res[,'lambda'],res[,'lambda']*0,type='l', lwd=1, col="gray")
    lines(res[,'lambda'],res[,'NB_all'],type='l', lwd=1, col="black")
    lines(res[,'lambda'],res[,'NB_max'],type='l',col="blue", lwd=2)
  }

  return(out)

}




#' @export
decision_curve <- function(y, pi, lambdas=(1:99)/100)
{
  NB_model <- rep(0, length(lambdas))
  NB_all <- NB_model

  for(j in 1:length(lambdas))
  {
    NB_model[j] <- mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
    NB_all[j] <- mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * 1)
  }

  return(cbind(lambda=lambdas, NB_none=NB_all*0, NB_model=NB_model, NB_all=NB_all))
}



#' @export
plot_decision_curve <- function(dc_data)
{
  max_y <-max(dc_data[,c('NB_none','NB_model','NB_all')])
  plot(dc_data[,'lambda'],dc_data[,'NB_none'],type='l',ylim=c(0,max_y),xlim=c(0,1),col='gray', xlab="Threshold", ylab="Net benefit")
  lines(dc_data[,'lambda'],dc_data[,'NB_all'],type='l',ylim=c(0,max_y),col='black')
  lines(dc_data[,'lambda'],dc_data[,'NB_model'],type='l',ylim=c(0,max_y),col='red',lw=2)
}


