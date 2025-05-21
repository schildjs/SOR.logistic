#' expit a log odds to obtain a probability
#' 
#' expit a log odds to obtain a probability.  This is the antilogit function
#' @param x a log odds
#' @return a probability
#' @export
expit <- function(x){ex <- exp(x)
                     return(ex/(1+ex))}

#' logit a probability to obtain a log odds
#' 
#' logit a probability to obtain a log odds
#' @param p a probability
#' @return a log odds
#' @export
logit <- function(p){return(log(p) - log(1-p))}

#' Sequential Offsetted (logistic) Regression
#' 
#' Conducts sequential offsetted logistic regression analyses
#' @param ymod The target response model
#' @param zmod The auxiliary variable model
#' @param id Subject identifier
#' @param pi1.pi0.ratio The ratio of sampling probabilities pi(z=1, X)/pi(Z=0, X)
#' @param DAT.ods The dataframe
#' @param CORSTR The correlation structure.  Presently only allows "independence", "exchangeabel" and "ar1".
#'              Note that one should only used "independence" if using observation level sampling.
#' @return beta: Parameter estimates for the target model.
#' @return robvarbeta: Robust covariance matrix for the target model
#' @return gamma: Parameter estimates for the auxiliary variable model
#' @return robvargamma: Robust covariance matrix for the auxiliary variable model
#' @export 
SOLR <-
  function(ymod, 
           zmod, 
           id,
           pi1.pi0.ratio,  
           DAT.ods, 
           CORSTR="independence"){
    DAT.ods$id       <- id
    DAT.ods$offset.z <- log(pi1.pi0.ratio) 
    #print(zmod)
    mod.z            <- glm(formula(paste(zmod[2],zmod[1],zmod[3], "+offset(offset.z)", sep="")), 
                            x=TRUE, data=DAT.ods, family="binomial")
    
    DAT.ods1     <- DAT.ods
    DAT.ods1[,paste(ymod[2])] <- 1
    DAT.ods0     <- DAT.ods
    DAT.ods0[,paste(ymod[2])] <- 0
    
    #DAT.ods1[,"Y"] <- 1
    #DAT.ods0     <- DAT.ods
    #DAT.ods0[,"Y"] <- 0
    z.formula <- zmod[c(1,3)]
    #print(z.formula)
    W0           <- model.matrix(z.formula, data=DAT.ods0)
    W1           <- model.matrix(z.formula, data=DAT.ods1)
    
    ## Observed situations where coefficient estimated for Y interactions in model for Z
    ## were very large so that when I calculated linear predictors, they were over 709
    ## which seems to be R's limit for exponentiation
    #print("blah0")
    lp.lambda1.s <- predict.glm(mod.z, DAT.ods1)
    lp.lambda0.s <- predict.glm(mod.z, DAT.ods0)
    #print("blah0")
    lambda1.ok <- lp.lambda1.s<710
    lambda0.ok <- lp.lambda0.s<710
    
    lambda1.p <- ifelse(lambda1.ok, expit(-1*DAT.ods$offset.z + lp.lambda1.s), 1)
    lambda0.p <- ifelse(lambda0.ok, expit(-1*DAT.ods$offset.z + lp.lambda0.s), 1)
    lambda1.s <- ifelse(lambda1.ok, expit(lp.lambda1.s ), 1)
    lambda0.s <- ifelse(lambda0.ok, expit(lp.lambda0.s ), 1)
    
    lp.lambda.s <- predict(mod.z) 
    lambda.ok   <- lp.lambda.s<710
    lambda.s    <- ifelse(lambda.ok, expit( lp.lambda.s) , 1)
    lambda.p    <- ifelse(lambda.ok, expit( -1*DAT.ods$offset.z + lp.lambda.s) , 1)
        
    rho1.scaled   <- (1-lambda1.p) + pi1.pi0.ratio*lambda1.p
    rho0.scaled   <- (1-lambda0.p) + pi1.pi0.ratio*lambda0.p
    DAT.ods$offset.y <- log(rho1.scaled/rho0.scaled)
    
    ## Calculate F_i (y) for y=1 and y=0
    Const <- 1- pi1.pi0.ratio
    F1    <- lambda1.p*(1-lambda1.p)* Const/(1-(Const*lambda1.p))
    F0    <- lambda0.p*(1-lambda0.p)* Const/(1-(Const*lambda0.p))
    
    ## Model for Y in pseudopopulation
    mod.y       <- geeglm(formula(paste(ymod[2],ymod[1],ymod[3], "+offset(offset.y)", sep="")), id=id, 
                          data=DAT.ods, x=TRUE, corstr=CORSTR, family="binomial")
    
    mu.s      <- mod.y$fitted.values
    cor.param <- mod.y$geese$alpha
    sca.param <- mod.y$geese$gamma
    
    p.z <- length(mod.z$coef)
    p.y <- length(mod.y$coef)
    p.all <- p.z+p.y
    
    W  <- mod.z$x
    Z  <- mod.z$y
    X  <- mod.y$x
    Y  <- mod.y$y
    
    ## Subject specific calculation
    I   <- matrix(0, p.all, p.all)
    Q   <- matrix(0, p.all, p.all)
    Itt <- matrix(0, p.z, p.z)
    Iuu <- matrix(0, p.y, p.y)
    T   <- rep(0, p.z)
    U   <- rep(0, p.y)
    Iut <- matrix(0, p.y, p.z)
    
    corstr <- CORSTR
    
    for (i in unique(DAT.ods$id)){ 
      
      keep <- (DAT.ods$id==i)    
      ni   <- sum(keep)
      Wi   <- W[keep,]
      Xi   <- X[keep,]
      Zi   <- Z[keep]
      Yi   <- Y[keep]
      
      F1i <- diag(F1[keep], ni, ni)
      F0i <- diag(F0[keep], ni, ni)
      
      W1i <- W1[keep,]
      W0i <- W0[keep,]
      
      lambda.si <- lambda.s[keep]
      lambda.pi <- lambda.p[keep]
      mu.si     <- mu.s[keep]
      
      ## Calculate working correlation and inverse of working covariance
      if (ni>1){ if (corstr=="independence"){          Ci      <- diag(rep(1, ni))
      }else if (corstr=="exchangeable"){ Ci      <- as.matrix(diag(rep(1-cor.param, ni), ni, ni)) + matrix(cor.param, ni, ni)
      }else if (corstr=="ar1"){          tmp.cor <- expand.grid(WAVE, WAVE)
                                         Ci      <- cor.param^abs(matrix(tmp.cor[,2]-tmp.cor[,1], length(WAVE)))
      }else{ stop("Need to use independence, exchangeable, or ar1 working covariance structure") }
      }else{ Ci <- 1}
      #print(mu.si)
      #print(Ci)
      
      Ai   <- diag((mu.si*(1-mu.si)), ni, ni)
      Vinv <- solve(sqrt(Ai) %*% Ci %*% sqrt(Ai))
      
      if (ni==1){ Wi  <- t(as.matrix(Wi)); Xi  <- t(as.matrix(Xi))
                  W0i <- t(as.matrix(W0i));W1i <- t(as.matrix(W1i)) }
      Di <- Ai %*% Xi
      Ti <- t(Wi) %*% (Zi - lambda.si)
      Ui <- t(Di) %*% Vinv %*% (Yi - mu.si)
      Qi <- outer(c(Ti,Ui), c(Ti,Ui))
      
      Itti <- matrix(0,p.z, p.z)
      for (j in 1:ni){ Itti <- Itti + (outer(Wi[j,], Wi[j,]) * (lambda.si[j]*(1-lambda.si[j]))) }
      
      Iuui <- t(Di) %*% Vinv %*% Di
      
      dBi.dgam <- t(t(W0i)%*%F0i - t(W1i)%*%F1i)
      Iuti     <- t(Di) %*% Vinv %*% Ai %*% dBi.dgam
      
      Iuu <- Iuu+Iuui
      Itt <- Itt+Itti
      Iut <- Iut+Iuti
      
      T <- T + Ti
      U <- U + Ui
      Q <- Q + Qi
    }
    
    Itu  <- matrix(0, p.z, p.y)
    I    <- rbind(cbind(Itt, Itu), cbind(Iut, Iuu))
    Iinv <- solve(I)
    AVAR <- Iinv %*% Q %*% t(Iinv)
    
    out <- list(beta=mod.y$coef, robvarbeta=AVAR[(p.z+1):(p.z+p.y),(p.z+1):(p.z+p.y)] ,
                gamma=mod.z$coef, robvargamma = AVAR[c(1:p.z), c(1:p.z)])
    out
  }
