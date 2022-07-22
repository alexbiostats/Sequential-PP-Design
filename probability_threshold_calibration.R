### File with code for calibration of decision rules based on a global null scenario (all 10 baskets with 10% response rate) with NO interim monitoring

source("~/function_mem_seq_bt.R")

rates <- rep(0.1,10) # set scenario to be null for all baskets with 10% response rate
numArm <- length(rates) # define number of baskets based on response rate


### Bayesian design with information sharing (UB=0.1)

p0 <- 0.1
pars <- list(LB=0, UB=0.1) #UB=0.1 to indicate potential for some borrowing if exchangeable
pp_mat_ub1 <- matrix(nrow=1000, ncol=10)

for(i in 1:1000){
	set.seed(i)

	dat <- sapply(1:numArm, function(x) rbinom(n=25, size=1, prob=rates[x])) # simulate trial data for all baskets
	X <- colSums(dat)
	N <- rep(25,numArm)
	Data <- list( X=X,N=N )

	post_prob <- postProb(Data=Data, pars=pars, p0=p0, pt=NULL)
	pp_mat_ub1[i,] <- post_prob	
}

# Alter values to achieve ~10% type I error rates (or whatever desired level is)

table(pp_mat_ub1 < 0.838) # marginal (i.e., basket-wise) type I error rate
table( rowSums(pp_mat_ub1 < 0.979) == 10) # family-wise type I error rate


### Bayesian design without information sharing (UB=0)

p0 <- 0.1
pars <- list(LB=0, UB=0) #UB=0 to indicate no information sharing
pp_mat_ub0 <- matrix(nrow=1000, ncol=10)

for(i in 1:1000){
	set.seed(i)

	dat <- sapply(1:numArm, function(x) rbinom(n=25, size=1, prob=rates[x])) # simulate trial data for all baskets
	X <- colSums(dat)
	N <- rep(25,numArm)
	Data <- list( X=X,N=N )

	post_prob <- postProb(Data=Data, pars=pars, p0=p0, pt=NULL)
	pp_mat_ub0[i,] <- post_prob	
}

# Alter values to achieve ~10% type I error rates (or whatever desired level is)

table(pp_mat_ub0 < 0.9) # marginal (i.e., basket-wise) type I error rate
table( rowSums(pp_mat_ub0 < 0.990) == 10) # family-wise type I error rate


### Simon's 2016 Posterior Probability Monitoring Design

p0 <- 0.1
p1 <- 0.3
pp_mat_simon <- matrix(nrow=1000, ncol=10)

for(i in 1:1000){
	set.seed(i)

	dat <- sapply(1:numArm, function(x) rbinom(n=25, size=1, prob=rates[x])) # simulate trial data for all baskets
	
	post_prob <- find_postk(lambda=0.1, gamma=0.33, r=colSums(dat), n=rep(25,10), plo=p0, phi=p1)
	pp_mat_simon[i,] <- post_prob	
}

# Alter values to achieve ~10% type I error rates (or whatever desired level is)

table(pp_mat_simon < 0.154) # marginal (i.e., basket-wise) type I error rate
table( rowSums(pp_mat_simon < 0.665) == 10) # family-wise type I error rate

