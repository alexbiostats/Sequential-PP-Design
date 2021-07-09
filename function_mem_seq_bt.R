### File with collection of functions used in various iterations to implement simulations for sequential predictive probability monitoring for a basket trial with 10 baskets

MEM.mat <- function(Indices, mod.mat, H){
  M <- matrix(NA, H, H ); diag(M) <- rep(1,dim(M)[1])
  for(i in 1:(length(Indices)-1)){ M[(i+1):dim(M)[2], i] <- M[i, (i+1):dim(M)[2]] <- mod.mat[[i]][Indices[i],] }
  M[dim(M)[2], i+1] <- M[i+1, dim(M)[2]] <- (0:1)[ Indices[length(Indices)] ]
  return(M)
}

rdirichlet <- function (n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

logMarg.Dens <- function(I, mod.mat, xvec, nvec, avec, bvec){
  M <- MEM.mat(I, mod.mat, length(xvec))
  marg.vec <- rep(NA, dim(M)[1])
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  for(i in 1:dim(M)[1]){ p.vec <- prod( prod.vec^(1-M[i,])  )
                         marg.vec[i] <- (beta(avec[i] + M[i,]%*%xvec, bvec[i] + M[i,]%*%(nvec-xvec)) / beta(avec[i],bvec[i]) ) * p.vec
  } #return( list(mod = M, marg = prod(marg.vec)) )
  return( sum(log(marg.vec)) )
}

logMarg.DensSA <- function(M, mod.mat, xvec, nvec, avec, bvec){
  #M <- MEM.mat(I, mod.mat, length(xvec))
  marg.vec <- rep(NA, dim(M)[1])
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  for(i in 1:dim(M)[1]){ p.vec <- prod( prod.vec^(1-M[i,])  )
                         marg.vec[i] <- (beta(avec[i] + M[i,]%*%xvec, bvec[i] + M[i,]%*%(nvec-xvec)) / beta(avec[i],bvec[i]) ) * p.vec
  } #return( list(mod = M, marg = prod(marg.vec)) )
  return( sum(log(marg.vec)) )
}

logMarg.DensINIT <- function(mod, xvec, nvec, a, b){
  M <- c(1,mod)
  prod.vec <- beta(xvec + a, nvec + b - xvec) / beta(a, b) #calculate the product portion of integrated marginal likelihood
  p.vec <- prod( prod.vec^(1-M)  )
  marg.vec <- (beta(a + M%*%xvec, b + M%*%(nvec-xvec)) / beta(a,b) ) * p.vec
  return( log(marg.vec) )
}

MEM_marginalINIT <- function(xvec, nvec, avec, bvec){
  ###Function to calculate the MEM marginal densities for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  
  mod.mat <- list(); k <- length(xvec)-1; j <- 1
  while(k > 0){ mod.mat[[j]] <- ( as.matrix(expand.grid( rep(list(c(0,1)), k) )) )[order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
                k <- k - 1; j <- j + 1 }
  #H <- length(mod.mat)
  #temp <- list(); for(h in 1:(H-1)){ temp[[h]] <- 1:dim(mod.mat[[h]])[1] }; temp[[H]] <- 1:2
  #Mod.I <- temp[[1]] #as.matrix(expand.grid( temp ))
  
  ## identify maximizer of marginal density ##
  MAX <- matrix(NA, length(xvec), length(xvec) ); diag(MAX) <- rep(1,dim(MAX)[1])
  log.Marg <- apply( mod.mat[[1]], MARGIN=1, FUN=logMarg.DensINIT, xvec, nvec, avec[1], bvec[1])
  o <- order(log.Marg,decreasing=TRUE)[1]
  MAX[1,] <- c(1, mod.mat[[1]][o[1],])
  K <- length(xvec)
  for(j in 2:(K-1)){ Ii <- c(j,1:(j-1),(j+1):K)
                     log.Marg <- apply( mod.mat[[1]], MARGIN=1, FUN=logMarg.DensINIT, xvec[Ii], nvec[Ii], avec[j], bvec[j])
                     o <- order(log.Marg,decreasing=TRUE)[1]
                     MAX[j,Ii] <- c(1, mod.mat[[1]][o[1],])
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))
  log.Marg <- apply( mod.mat[[1]], MARGIN=1, FUN=logMarg.DensINIT, xvec[Ii], nvec[Ii], avec[j], bvec[j])
  o <- order(log.Marg,decreasing=TRUE)[1]
  MAX[j,Ii] <- c(1, mod.mat[[1]][o[1],])
  colnames(MAX) <- rownames(MAX) <- xvec
  
  ###Create list to return matrix which identifies which sources are included in each model and the marginal densities
  ret <- list(mod.mat = mod.mat, maximizer=MAX)
  return(ret)
}


# Fast optimizer
MEM_marginalLoop <- function(xvec, nvec, avec, bvec, INITIAL){
  ###Function to calculate the MEM marginal densities for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  
  mod.mat <- list(); k <- length(xvec)-1; j <- 1
  while(k > 0){ mod.mat[[j]] <- ( as.matrix(expand.grid( rep(list(c(0,1)), k) )) )[order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
                k <- k - 1; j <- j + 1 }

 
  if( is.na(INITIAL[1]) ){ M.init <-matrix(1,length(xvec),length(xvec))}

 
  #logDenOld<-logMarg.DensSA(M.init, mod.mat, xvec, nvec, avec, bvec)
  #M[1,2]<-1-M[1,2]
  MOld<-M.init  
  d<-length(xvec)   
  for( i in 1:(d-1))
    for(k in (i+1):d)
    {
      if (MOld[i,k]!=MOld[k,i])
      {
        MOld[i,k]<-MOld[k,i]<-rbinom(1,1,0.5)
      }
    }
  logDenOld<-logMarg.DensSA(MOld, mod.mat, xvec, nvec, avec, bvec)
  for(iter in 1:1000)
  {
    cat(iter,"\n")
    changed<-FALSE
    for( i in 1:(d-1))
      for(k in (i+1):d)
      {
        loc<-c(k,i) 
        MNew<-MOld
        MNew[loc[1],loc[2]]<-MNew[loc[2],loc[1]]<-1-MOld[loc[1],loc[2]]
        logDenNew<-logMarg.DensSA(MNew, mod.mat, xvec, nvec, avec, bvec)
        if (logDenNew > logDenOld)
        {
          MOld<-MNew
          changed<-TRUE
          logDenOld<-logDenNew
          cat(logDenOld, "\n")
        }
      }
    if (!changed)
    {
      break
    }
  }
  
  ## identify maximizer of marginal density ##
  #log.Marg <- apply( Mod.I, MARGIN=1, FUN=logMarg.Dens, mod.mat, xvec, nvec, avec, bvec)
  
  #browser()
  
  colnames(M.init) <- rownames(M.init) <- xvec
  colnames(MOld) <- rownames(MOld) <- xvec
  cat("final density: ", logDenOld, "\n")
#browser()
  return(list(mod.mat = mod.mat, init=M.init, maximizer=MOld))
}

MEM_marginal <- function(xvec, nvec, avec, bvec){
  ###Function to calculate the MEM marginal densities for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  
  mod.mat <- list(); k <- length(xvec)-1; j <- 1
  while(k > 0){ mod.mat[[j]] <- ( as.matrix(expand.grid( rep(list(c(0,1)), k) )) )[order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
                k <- k - 1; j <- j + 1 }
  H <- length(mod.mat)
  temp <- list(); for(h in 1:(H-1)){ temp[[h]] <- 1:dim(mod.mat[[h]])[1] }; temp[[H]] <- 1:2
  Mod.I <- as.matrix(expand.grid( temp ))
  #browser()
  ## identify maximizer of marginal density ##
  log.Marg <- apply( Mod.I, MARGIN=1, FUN=logMarg.Dens, mod.mat, xvec, nvec, avec, bvec)
  o <- order(log.Marg,decreasing=TRUE)[1]
  print(log.Marg)
  MAX <- matrix(NA, length(xvec), length(xvec) ); diag(MAX) <- rep(1,dim(MAX)[1])
  for(i in 1:(length(Mod.I[o,])-1)){ MAX[(i+1):dim(MAX)[2], i] <- MAX[i, (i+1):dim(MAX)[2]] <- mod.mat[[i]][Mod.I[o,i],] }
  MAX[dim(MAX)[2], i+1] <- MAX[i+1, dim(MAX)[2]] <- (0:1)[ Mod.I[o,length(Mod.I[o,])] ]
  colnames(MAX) <- rownames(MAX) <- xvec
  
  ###Create list to return matrix which identifies which sources are included in each model and the marginal densities
  ret <- list(mod.mat = mod.mat, maximizer=MAX)
  return(ret)
}

MEM_modweight <- function(mod.mat, source.vec){
  ###Function to calculate model weights for each MEM given source inclusion probabilities
  # Applied to each row of the symmetric source inclusion prior probabilities
  #mod.mat: matrix which specifies which supplemental sources are included for analysis of each subtype
  #source.mat: rows of symmetric matrix of source inclusion probabilties for MEM
  
  s.in <- source.vec #source inclusion probability
  s.ex <- 1-source.vec #source exclusion probability
  
  ###Calculate model weights
  q.mat <- sapply( 1:dim(mod.mat)[1], function(x) s.in^(mod.mat[x,]) * s.ex^(1-mod.mat[x,]))
  
  if(length(s.in)==1){ q.vec <- q.mat #if only one supplemental sources, q.mat is equal to q.vec for model weights
  }else{ q.vec <- apply( q.mat, 2, prod ) } #if more than one supplemental source, take product of source probabilities for each model
  
  return(q.vec)
}

eval.Post <- function(p0,X,N,Omega,w){
  a=b=0.5 
  alph <- a + Omega%*%X
  beta <- b + (Omega%*%N-Omega%*%X)
  return( sum( (1-pbeta(p0, alph, beta))*w ) )
}

MEM.cdf <- function(Data,pars,p0,marg.M){
  #marg.M <- MEM_marginal(xvec=Data$X, nvec=Data$N, avec=rep(0.5,length(Data$N)), bvec=rep(0.5,length(Data$N)))
  pr.Inclus <- pars$UB*marg.M$maximizer; pr.Inclus[which(pr.Inclus==0)] <- pars$LB
  weights <- list() 
  for(i in 1:nrow(pr.Inclus)){ 
    weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], source.vec=pr.Inclus[i,][-i]) 
  }
  models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
  out <- eval.Post(p0, Data$X, Data$N, models, weights[[1]])
  K <- length(Data$X)
  for(j in 2:(K-1)){ Ii <- c(j,1:(j-1),(j+1):K)
                     out <- c(out, eval.Post(p0, Data$X[Ii], Data$N[Ii], models, weights[[j]]))
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))
  out <- c(out, eval.Post(p0, Data$X[Ii], Data$N[Ii], models, weights[[j]])) #; print(pr.Inclus)
  return(out)
}

PostProb.mlab <- function(Data,pars,p0){
  pp.mem <- MEM.cdf(Data,pars,p0)
  if(is.na(pars$ESS)){ out <- sum( pars$m.BasketWeights*pp.mem )
  } else{ out <- rdirichlet(1, pars$m.BasketWeights*pars$ESS)*pp.mem }
  return(out)
}

genObs.mlab <- function(hyp,N,pars,p0){
  D <- list(X=rbinom(rep(1,length(N)),N,hyp), N=N)
  return( PostProb.mlab(D,pars,p0) )
}

checkThres.mlab <- function(THRES,null,N,pars,p0){ N.rep <- 5000
                                                   return( sum( replicate(N.rep, genObs.mlab(null,N,pars,p0) ) > THRES )/N.rep )
}

CalibrateTypeIError.mlab <- function(null,Grid,N,pars,p0){
  typeIerror <- sapply(X=Grid,FUN=checkThres.mlab,null,N,pars,p0)
  return(list(Threshold=Grid, TypeIerror=typeIerror)) 
}

Power.mlab <- function(ALT,thres,N,pars,p0){ N.rep <- 5000
                                             return( list( PostProbThreshold=thres, power=sum( replicate(N.rep, genObs.mlab(ALT,N,pars,p0) ) > thres )/N.rep ) )
}

genObs.bask <- function(THRES,hyp,N,pars,p0){
  D <- list(X=rbinom(rep(1,length(N)),N,hyp), N=N)
  return( as.numeric( MEM.cdf(D,pars,p0) > THRES ) )
}

checkThres.bask <- function(THRES,null,N,pars,p0){ N.rep <- 5000 ## Compute FWER ##
                                                   return(  1 - ( sum( rowSums( t( replicate(N.rep, genObs.bask(THRES,null,N,pars,p0) ) ) ) == 0 ) / N.rep ) )
}

CalibrateFWER.bask <- function(null,Grid,N,pars,p0){
  FWER <- sapply(X=Grid,FUN=checkThres.bask,null,N,pars,p0)
  return(list(Threshold=Grid, FWER=FWER)) 
}

Power.bask <- function(ALT,thres,N,pars,p0){ N.rep <- 5000
                                             return( list( PostProbThreshold=thres, 
                                                           power=colSums( t( replicate(N.rep, genObs.bask(thres,ALT,N,pars,p0) ) ) )/N.rep ) )
}

# Return the predicted final alpha and beta values of the basket with the index, including the prediction data
getFinal<-function(Data, model, weights, index)
{
  Omega<-model
  if (index==1)
  {
    a=b=0.5 
    alph <- a + Omega%*%Data$X
    beta <- b + (Omega%*%Data$N-Omega%*%Data$X)
    return(list(alph=alph, beta=beta))
  }
  
  K <- length(Data$X)
  if (index < K)
  {  
    j <- index
    Ii <- c(j,1:(j-1),(j+1):K)
    a=b=0.5 
    alph <- a + Omega%*%Data$X[Ii]
    beta <- b + (Omega%*%Data$N[Ii]-Omega%*%Data$X[Ii])
    return(list(alph=alph, beta=beta))
  }
  j <- K
  Ii <- c(j,1:(j-1))
  a=b=0.5 
  alph <- a + Omega%*%Data$X[Ii]
  beta <- b + (Omega%*%Data$N[Ii]-Omega%*%Data$X[Ii])
  return(list(alph=alph, beta=beta))
}

# Return the alpha and beta values of the basket with the index based on the current data
getAlphaBeta<- function(Data,pars,index){
  marg.M <- MEM_marginalLoop(xvec=Data$X, nvec=Data$N, avec=rep(0.5,length(Data$N)), bvec=rep(0.5,length(Data$N)),INITIAL=NA)
  pr.Inclus <- pars$UB*marg.M$maximizer; pr.Inclus[which(pr.Inclus==0)] <- pars$LB
  weights <- list() 
  for(i in 1:nrow(pr.Inclus)){ 
    weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], source.vec=pr.Inclus[i,][-i]) 
  }
  
  model <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
  Omega<-model
  if (index==1)
  {
    a=b=0.5 
    alph <- a + model%*%Data$X
    beta <- b + (Omega%*%Data$N-Omega%*%Data$X)
    return(list(alph=alph, beta=beta))
  }
  
  K <- length(Data$X)
  if (index < K)
  {  
    j <- index
    Ii <- c(j,1:(j-1),(j+1):K)
    a=b=0.5 
    alph <- a + model%*%Data$X[Ii]
    beta <- b + (Omega%*%Data$N[Ii]-Omega%*%Data$X[Ii])
    return(list(alph=alph, beta=beta))
  }
  j <- K
  Ii <- c(j,1:(j-1))
  a=b=0.5 
  alph <- a + model%*%Data$X[Ii]
  beta <- b + (Omega%*%Data$N[Ii]-Omega%*%Data$X[Ii])
  return(list(alph=alph, beta=beta))
}


# Predictive probability
predProb <- function(Data,pars, remain, p0, pt){
##### predProb function modified by Alex to account for arms that have stopped enrolling
### Data: X (events), N (sample size)
### pars: LB, UB parameters for MEM prior
### remain: number to still enroll in a given basket
### p0: null value
### pt: threshold

	marg.M <- MEM_marginalLoop(xvec=Data$X, nvec=Data$N, avec=rep(0.5,length(Data$N)), bvec=rep(0.5,length(Data$N)),INITIAL=NA)
	print(marg.M)
	pr.Inclus <- pars$UB*marg.M$maximizer; pr.Inclus[which(pr.Inclus==0)] <- pars$LB
	weights <- list() 

	for(i in 1:nrow(pr.Inclus)){ 
		weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], source.vec=pr.Inclus[i,][-i]) 
	}

	models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
	remain<- remain #old version: round(totalRemain * Data$N / sum(Data$N))
	prob<-c()
	K <- length(Data$X)

	# All alpha beta values
	allAlphBeta<-vector(mode = "list", length = K)
	for (i in 1:K){
		ab<-getAlphaBeta(Data, pars, i)
		allAlphBeta[[i]]<-ab
	}

	# Generate prediction data for K baskets that are still enrolling
	rep<-5000
	count<-rep(0, K)
	allAB<-c()
	allM<-c()
	allP<-c()

	for (j in 1:rep){
		allSucc<-rep(0,K)

		# Successes for all baskets
		for (i in which(remain!=0) ){ #replace 1:K with "which(remain!=0)"
			ab<-allAlphBeta[[i]]
			allR<-c()
			sel<-rmultinom(1, 1, weights[[i]])
			k<-which(sel==1)[1]
			rate<-rbeta(1,ab$alph[k], ab$beta[k])
      
			# Prediction
			allSucc[i] <- rbinom(1,remain[i],rate)
		}

		newData<-Data
		for (i in which(remain!=0) ){
			newData$N[i]<-newData$N[i]+remain[i]
			newData$X[i]<-newData$X[i]+allSucc[i]
		}

		# Final inference
		for (i in which(remain!=0) ){
			abNew<-getFinal(newData, models, weights, i) # weights appears to be an unnecessary term not used in the function
			alph<-abNew$alph
			beta<-abNew$beta
			pp<-sum((1-pbeta(p0, alph, beta))*weights[[i]] )

			# Sample size
			count[i]<-count[i]+(pp>pt)
		} 
	}
	prob<-count/rep
	return(prob)
}



# Posterior probability
postProb <- function(Data,pars, p0, pt){
##### posterior probability calculation for after all enrollment is complete, modification of predProb
### Data: X (events), N (sample size)
### pars: LB, UB parameters for MEM prior
### p0: null value
### pt: threshold

	marg.M <- MEM_marginalLoop(xvec=Data$X, nvec=Data$N, avec=rep(0.5,length(Data$N)), bvec=rep(0.5,length(Data$N)),INITIAL=NA)
	print(marg.M)
	pr.Inclus <- pars$UB*marg.M$maximizer; pr.Inclus[which(pr.Inclus==0)] <- pars$LB
	weights <- list() 

	for(i in 1:nrow(pr.Inclus)){ 
		weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], source.vec=pr.Inclus[i,][-i]) 
	}

	models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
	prob<-c()
	K <- length(Data$X)

	# All alpha beta values
	allAlphBeta<-vector(mode = "list", length = K)
	for (i in 1:K){
		ab<-getAlphaBeta(Data, pars, i)
		allAlphBeta[[i]]<-ab
	}

	# Final inference
	pp_vec <- rep(NA,K)

	for (i in 1:K ){
		abNew<-getFinal(Data, models, weights, i)
		alph<-abNew$alph
		beta<-abNew$beta
		pp_vec[i] <-sum((1-pbeta(p0, alph, beta))*weights[[i]] )
	} 
	return(pp_vec)
}


trial_sim <- function(rates, numPat=25, start=5, UB, p0=0.1, postp, prep, seed, filename){
### Function to simulate trial based on given inputs
# rates: vector of rates for baskets
# numPat: max number of patients
# start: number of patients to wait for before starting monitoring
# UB: upper boundary borrowing parameter
# p0: null response
# postp: posterior probability threshold to use for efficacy
# prep: predictive probability threshold to use for futility monitoring
# seed: seed to set for simulation, default is 515

	set.seed(seed)

	numArm <- length(rates)
	dat <- sapply(1:numArm, function(x) rbinom(n=numPat, size=1, prob=rates[x])) # simulate trial data for all baskets

	pars <- list(LB=0, UB=UB) ## MMLE lower, upper boundary
  
	if(prep==0){ 
		X <- colSums(dat)
		N <- rep(numPat,numArm)
		Data <- list( X=X,N=N )
	}else{

		# Calculate number of events and sample size at start
		X <- colSums( dat[1:(start-1),] )
		N <- rep( (start-1), numArm )
		be <- rep(TRUE, numArm) # indicator for if a basket is still enrolling
  
		# Loop through each enrollment and stop baskets as needed
		for (j in start:numPat){

			# Update date for jth enrollment if basket still enrolling
			X[which(be==T)] <- sapply(which(be==T), function(z) if(be[z]==T){ X[z] <- X[z] + dat[j,z]})
			N[which(be==T)] <- sapply(which(be==T), function(z) if(be[z]==T){ N[z] <- N[z] + 1})
			Data <- list( X=X, N=N ) ## Vectors of Observed number of Responses (X) and Patients (N)

			if( j == numPat ){ break } # break to calculate posterior probability after loop
			Nremain <- sapply(1:numArm, function(z) if(be[z]==T){ numPat-N[z] }else{ 0 }) # max number left to enroll
			r <- predProb(Data=Data, pars=pars, remain=Nremain, p0=p0, pt=postp)

			be[which(be==T)] <- sapply(which(be==T), function(z) if(be[z]==T){ be[z] <- r[z] >= prep }) # if pred. prob >= prep (threshold), keep enrolling
		
			if( sum(be==F)==numArm ){ break } # break if all arms have stopped enrollment
		}
	}

	post_prob <- postProb(Data=Data, pars=pars, p0=p0, pt=postp)
	reject_ind <- post_prob >= postp
	reject_ind[ which(N < numPat) ] <- 0 #fix added after seeing results where early stop could also result in reject pending postp
	earlystop_ind <- (N < numPat)

	ret_vec <- c('seed'=seed, 'UB'=UB, 'p0'=p0, 'prep'=prep, 'postp'=postp, 'rate'=rates,'X'=X, 'N'=N, 'PP'=post_prob, 'reject'=reject_ind, 'earlystop'=earlystop_ind) 
	write.table(matrix(ret_vec,nrow=1), file=filename, col.names = F, row.names = F, append = TRUE)
}


trial_simon_sim <- function(rates, numPat=25, start=5, UB='s', p0=0.1, postp=NA, prep=NA, seed, filename){
### Function to simulate trial based on given inputs
# rates: vector of rates for baskets
# numPat: max number of patients
# start: number of patients to wait for before starting monitoring
# UB: upper boundary borrowing parameter
# p0: null response
# postp: posterior probability threshold to use for efficacy
# prep: predictive probability threshold to use for futility monitoring
# seed: seed to set for simulation, default is 515

	set.seed(seed)

	numArm <- length(rates)
	dat <- sapply(1:numArm, function(x) rbinom(n=numPat, size=1, prob=rates[x])) # simulate trial data for all baskets

	X1 <- colSums(dat[1:16,])
	be <- (X1 > 1) #keep enrolling (<=1 indicates stop enrollment)

	X <- sapply(1:numArm, function(x) if(be[x]==T){sum(dat[,x])}else{X1[x]})
	N <- sapply(1:numArm, function(x) if(be[x]==T){25}else{16})

	post_prob <- rep(NA,numArm)
	reject_ind <- as.numeric( sapply(1:numArm, function(x) if(be[x]==T){X[x] > 4}else{FALSE}) )
	earlystop_ind <- as.numeric(N==16)

	ret_vec <- c('seed'=seed, 'UB'=UB, 'p0'=p0, 'prep'=prep, 'postp'=postp, 'rate'=rates,'X'=X, 'N'=N, 'PP'=post_prob, 'reject'=reject_ind, 'earlystop'=earlystop_ind) 
	write.table(matrix(ret_vec,nrow=1), file=filename, col.names = F, row.names = F, append = TRUE)
}

