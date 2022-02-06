#' @import mvtnorm
#' @import glmnet
#' @import progress

aux <- new.env()

#' @export
get_aux<-function()
{
  return(aux)
}



#' @title  voi_glmn function
#' @param n_sim Number of simulations required. 
#' @param glm_object The GLM object representing the proposed model
#' @param thresholds thresholds at which EVPI is calculated
#' @param mc_type any of "bootstrap", "Bayesian_bootstrap", or "likelihood"
#' @param data for EVPI calculations. 
#' @param pi A vector of predicted values. Required only if data is provided. If a string or single number it indicates the corresponding columns in data
#' @param truth_formula Formula of the correct model. Its parameters will be estimated repeatedly using the data
#' @param family GLM family and link function 
#' @export
voi_glm <- function(n_sim, thresholds=(0:99)/100, glm_object=NULL, mc_type="bootstrap", data=NULL, pi=NULL, truth_formula=NULL, family=binomial(link="logit"))
{
  if(is.null(data))
    if(is.null(glm_object))
      stop("both glm_object and data were NULL. One is needed.\n")
    else
      data <- glm_object$data
    
  
  if(is.null(pi))
    if(is.null(glm_object))
      stop("both glm_object and model_formula were NULL. One is needed.\n")
    else
      pi <- predict(glm_object,newdata = data, type="response")
  else
  {
    if(length(pi)==1) pi <- data[,pi]
  }
  
  if(is.null(truth_formula))
    if(is.null(glm_object))
      stop("both glm_object and turth_formula were NULL. One is needed.\n")
  else
    truth_formula <- glm_object$call$formula
  
  if(is.null(family))
    if(is.null(glm_object))
      stop("both glm_object and turth_formula were NULL. One is needed.\n")
  else
    family <- glm_object$call$family
  
  n <- dim(data)[1]
  
  if(mc_type=="bootstrap")  
  {
    NB_all <- NB_model <- NB_max <- rep(0,length(thresholds))
    for(i in 1:n_sim)
    {
      bs_data <- data[sample(1:n,n,replace = T),]
      bs_model <- glm(formula = truth_formula, family=family, data = bs_data)
      p <- predict(bs_model, newdata=data, type="response")
      
      for(j in 1:length(thresholds))
      {
        NB_all[j] <- NB_all[j] + mean((p-(1-p)*thresholds[j]/(1-thresholds[j])))
        NB_model[j] <- NB_model[j] + mean((pi>thresholds[j])*(p-(1-p)*thresholds[j]/(1-thresholds[j])))
        NB_max[j] <- NB_max[j] + mean((p>thresholds[j])*(p-(1-p)*thresholds[j]/(1-thresholds[j])))
      }
    }
  }
  
  
  if(mc_type=="Bayesian_bootstrap")  
  {
    NB_all <- NB_model <- NB_max <- rep(0,length(thresholds))
    for(i in 1:n_sim)
    {
      w <- c(0,sort(runif(n-1)),1)
      data$w <- w[-1]-w[-length(w)]
      bs_model <- glm(formula = truth_formula, family=family, data = data, weights = w)
      p <- predict(bs_model, newdata=data, type="response")
      
      for(j in 1:length(thresholds))
      {
        NB_all[j] <- NB_all[j] + mean((p-(1-p)*thresholds[j]/(1-thresholds[j])))
        NB_model[j] <- NB_model[j] + mean((pi>thresholds[j])*(p-(1-p)*thresholds[j]/(1-thresholds[j])))
        NB_max[j] <- NB_max[j] + mean((p>thresholds[j])*(p-(1-p)*thresholds[j]/(1-thresholds[j])))
      }
    }
  }
  
  
  if(mc_type=="likelihood")  
  {
    NB_all <- NB_model <- NB_max <- rep(0,length(thresholds))
    mu <- coefficients(glm_object)
    covmat <- vcov(glm_object)
    for(i in 1:n_sim)
    {
      betas <- rmvnorm(1,mu,covmat)
      model$coefficients <- betas
      p <- predict(model, newdata=data, type="response")
      
      for(j in 1:length(thresholds))
      {
        NB_all[j] <- NB_all[j] + mean((p-(1-p)*thresholds[j]/(1-thresholds[j])))
        NB_model[j] <- NB_model[j] + mean((pi>thresholds[j])*(p-(1-p)*thresholds[j]/(1-thresholds[j])))
        NB_max[j] <- NB_max[j] + mean((p>thresholds[j])*(p-(1-p)*thresholds[j]/(1-thresholds[j])))
      }
    }
  }
  
  NB_all <- NB_all / n_sim
  NB_model <- NB_model / n_sim
  NB_max <- NB_max / n_sim
  
  INB_current <- pmax(0,NB_all,NB_model) - pmax(0,NB_all)
  INB_perfect <- NB_max - pmax(0,NB_all)
  
  EVPI <- INB_perfect - INB_current
  EVPIr <- INB_perfect / INB_current
  
  return(data.frame(threshold=thresholds, EVPI=EVPI, EVPIr=EVPIr, INB_perfect=INB_perfect, INB_current=INB_current, NB_all=NB_all, NB_model=NB_model, NB_max=NB_max))
}





#  voi.glmnet function
# @param reg_obj: any object that you can apply predict with new data to get predictions. The structure should represent the correct model
# @param x: The model matrix of predictors
# @param y: The vector of responses
# @param pi: optional. Predictions from the current model. If not supplied, the predictions from reg_object will be used
#' @export
voi_glmnet <- function(reg_obj, x, y, pi=NULL, n_sim=1000, lambdas=(1:99)/100, Bayesian_bootstrap=F, empirical=F)
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
voi.glmnet2 <- function(formula, data, pi, n_sim=1000, lambdas=(1:99)/100, Bayesian_bootstrap=F, weights=NULL)
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
plot_evpir <- function(EVPIr, lambdas=(0:99)/100, max_y=10, ...)
{
  args <- list(...)
  if(!is.null(args$col))
  {
    if(length(args$col)<3) stop("If providing color, 3 should be mentioned (for the curve, 00, and inf).")
    cols <- args$col
  }
  else
    cols <- c("red","grey","black")
  
  max_y <- min(max(EVPIr,na.rm = T), max_y)
  yNaN <- rep(0,length(lambdas))
  yNaN[which(is.nan(EVPIr))] <- 1
  yInf <- rep(0,length(lambdas))
  yInf[which(EVPIr>10)] <- 1
  w <- rep(1/length(lambdas),length(lambdas))
  plot(z, EVPIr, type='l', col=cols[1], ylim=c(0,max_y), xlab="Threshold", ylab="Relative EVPI")
  par(new=T)
  barplot(yNaN, w, border=cols[2], col=cols[2], xlim=c(0,1), ylim=c(0,max_y), xlab=NULL, ylab=NULL, space=0, axes=FALSE)
  par(new=T)
  barplot(yInf, w, border=cols[3], col=cols[3], xlim=c(0,1), ylim=c(0,max_y), xlab=NULL, ylab=NULL, space=0, axes=FALSE)
  legend(0.8, max_y-1, legend=c("0/0",paste0(">",max_y)),col=c(cols[2],cols[3]), lty=c(1,1), lwd=c(10,10), border=NA)
}


