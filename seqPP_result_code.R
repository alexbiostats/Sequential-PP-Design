### Code to summarize results and create figures for the sequential PP paper based on simulation results from snowfall_sim_memseqbt.R

library(xtable)

dat <- read.table('~/mem_seq_bt_results.txt', header=T, stringsAsFactors=F)
dat$scen <- apply(dat[,paste0('rate',1:10)] ,1,paste0,collapse="/")


### Functions to synthesize the results

sim_sum <- function(res, nullp=0.1, altp=0.3, scen, postp0, postp1){
### Function to create table summarizing results for B=0 and B=0.1
# res: file with results
# nullp: response rate in null baskets (to be used in identifying match to null_ind)
# altp: response rate in alternative baskets (to be used in identifying match to alt_ind)
# scen: vector of 0.1 or 0.3 values for each basket (use all 0.1 for global null/alternative)
# postp0: posterior probability threshold used for UB=0
# postp1: posterior probability threshold used for UB=0.1

	#Extract null and alternative data sets (written 
	scen_null <- paste0(scen,collapse='/')
	null_basket <- which(scen == nullp)

	scen_alt <- if( sum(scen==nullp) == length(scen) ){ paste0(rep(altp,length(scen)), collapse='/') }else{ scen_null }
	alt_basket <- if( sum(scen==nullp) == length(scen) ){ null_basket }else{ which(scen == altp) }
	
	asd0n <- res[which(res$postp==postp0 & res$scen==scen_null & res$UB==0),]
	asd1n <- res[which(res$postp==postp1 & res$scen==scen_null & res$UB==0.1),]
	asd0a <- res[which(res$postp==postp0 & res$scen==scen_alt & res$UB==0),]
	asd1a <- res[which(res$postp==postp1 & res$scen==scen_alt & res$UB==0.1),]

	#Create object to store results
	col_names <- c('rr_n_0','rr_n_1','fw_n_0','fw_n_1','es_n_0','es_n_1','sr_n_0','sr_n_1','as_n_0','as_n_1','rr_a_0','rr_a_1','es_a_0','es_a_1','sr_a_0','sr_a_1','as_a_0','as_a_1')
	mat_res <- matrix( nrow=11, ncol=18, dimnames=list(seq(0,0.5,by=0.05),col_names)) 

	for( i in c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)){ 
		rr_n_0 <- mean( unlist(asd0n[which(asd0n$prep==i),paste0('reject',null_basket)]) )
		fw_n_0 <- mean( rowSums(asd0n[which(asd0n$prep==i),paste0('reject',null_basket)]) > 0 )
		es_n_0 <- mean( unlist(asd0n[which(asd0n$prep==i),paste0('N',null_basket)]) )
		sr_n_0 <- mean( unlist(asd0n[which(asd0n$prep==i),paste0('earlystop',null_basket)]) )
		as_n_0 <- mean( rowSums(asd0n[which(asd0n$prep==i),paste0('earlystop',null_basket)])==length(null_basket) )

		rr_n_1 <- mean( unlist(asd1n[which(asd1n$prep==i),paste0('reject',null_basket)]) )
		fw_n_1 <- mean( rowSums(asd1n[which(asd1n$prep==i),paste0('reject',null_basket)]) > 0 )
		es_n_1 <- mean( unlist(asd1n[which(asd1n$prep==i),paste0('N',null_basket)]) )
		sr_n_1 <- mean( unlist(asd1n[which(asd1n$prep==i),paste0('earlystop',null_basket)]) )
		as_n_1 <- mean( rowSums(asd1n[which(asd1n$prep==i),paste0('earlystop',null_basket)])==length(null_basket) )

		rr_a_0 <- mean( unlist(asd0a[which(asd0a$prep==i),paste0('reject',alt_basket)]) )
		es_a_0 <- mean( unlist(asd0a[which(asd0a$prep==i),paste0('N',alt_basket)]) )
		sr_a_0 <- mean( unlist(asd0a[which(asd0a$prep==i),paste0('earlystop',alt_basket)]) )
		as_a_0 <- mean( rowSums(asd0a[which(asd0a$prep==i),paste0('earlystop',alt_basket)])==length(alt_basket) )

		rr_a_1 <- mean( unlist(asd1a[which(asd1a$prep==i),paste0('reject',alt_basket)]) )
		es_a_1 <- mean( unlist(asd1a[which(asd1a$prep==i),paste0('N',alt_basket)]) )
		sr_a_1 <- mean( unlist(asd1a[which(asd1a$prep==i),paste0('earlystop',alt_basket)]) )
		as_a_1 <- mean( rowSums(asd1a[which(asd1a$prep==i),paste0('earlystop',alt_basket)])==length(alt_basket) )

		mat_res[paste0(i),] <- c(rr_n_0,rr_n_1,fw_n_0,fw_n_1,es_n_0,es_n_1,sr_n_0,sr_n_1,as_n_0,as_n_1,rr_a_0,rr_a_1,es_a_0,es_a_1,sr_a_0,sr_a_1,as_a_0,as_a_1)
	}

	return(mat_res)

}

simon_sum <- function(res){
### Function to create table summarizing results Simon two-stage design
# res: file with results
# nullp: response rate in null baskets (to be used in identifying match to null_ind)
# altp: response rate in alternative baskets (to be used in identifying match to alt_ind)
# scen: vector of 0.1 or 0.3 values for each basket (use all 0.1 for global null/alternative)

	asd1 <- res[which(res$scen=='0.1/0.1/0.1/0.1/0.1/0.1/0.1/0.1/0.1/0.1' & res$UB=='s'),]
	asd2 <- res[which(res$scen=='0.3/0.3/0.3/0.3/0.3/0.3/0.3/0.3/0.3/0.3' & res$UB=='s'),]
	asd3 <- res[which(res$scen=='0.1/0.1/0.1/0.1/0.1/0.1/0.1/0.1/0.3/0.3' & res$UB=='s'),]
	asd4 <- res[which(res$scen=='0.1/0.1/0.1/0.1/0.1/0.3/0.3/0.3/0.3/0.3' & res$UB=='s'),]

	#Create object to store results
	col_names <- c('rr_n','fw_n','es_n','sr_n','as_n','rr_a','es_a','sr_a','as_a')
	mat_res <- matrix( nrow=3, ncol=9, dimnames=list(c('global','mixed_8null_2alt','mixed_5null_5alt'),col_names)) 

	# Global scenarios
	rr_n_0 <- mean( unlist(asd1[,paste0('reject',1:10)]) )
	fw_n_0 <- mean( rowSums(asd1[,paste0('reject',1:10)]) > 0 )
	es_n_0 <- mean( unlist(asd1[,paste0('N',1:10)]) )
	sr_n_0 <- mean( unlist(asd1[,paste0('earlystop',1:10)]) )
	as_n_0 <- mean( rowSums(asd1[,paste0('earlystop',1:10)])==10 )

	rr_a_0 <- mean( unlist(asd2[,paste0('reject',1:10)]) )
	es_a_0 <- mean( unlist(asd2[,paste0('N',1:10)]) )
	sr_a_0 <- mean( unlist(asd2[,paste0('earlystop',1:10)]) )
	as_a_0 <- mean( rowSums(asd2[,paste0('earlystop',1:10)])==10 )

	mat_res['global',] <- c(rr_n_0,fw_n_0,es_n_0,sr_n_0,as_n_0,rr_a_0,es_a_0,sr_a_0,as_a_0)

	# Mixed scenario: 8 null, 2 alt baskets
	rr_n_m <- mean( unlist(asd3[,paste0('reject',1:8)]) )
	fw_n_m <- mean( rowSums(asd3[,paste0('reject',1:8)]) > 0 )
	es_n_m <- mean( unlist(asd3[,paste0('N',1:8)]) )
	sr_n_m <- mean( unlist(asd3[,paste0('earlystop',1:8)]) )
	as_n_m <- mean( rowSums(asd3[,paste0('earlystop',1:8)])==8 )

	rr_a_m <- mean( unlist(asd3[,paste0('reject',9:10)]) )
	es_a_m <- mean( unlist(asd3[,paste0('N',9:10)]) )
	sr_a_m <- mean( unlist(asd3[,paste0('earlystop',9:10)]) )
	as_a_m <- mean( rowSums(asd3[,paste0('earlystop',9:10)])==2 )

	mat_res['mixed_8null_2alt',] <- c(rr_n_m,fw_n_m,es_n_m,sr_n_m,as_n_m,rr_a_m,es_a_m,sr_a_m,as_a_m)

	# Mixed scenario: 5 null, 5 alt baskets
	rr_n_m <- mean( unlist(asd4[,paste0('reject',1:5)]) )
	fw_n_m <- mean( rowSums(asd4[,paste0('reject',1:5)]) > 0 )
	es_n_m <- mean( unlist(asd4[,paste0('N',1:5)]) )
	sr_n_m <- mean( unlist(asd4[,paste0('earlystop',1:5)]) )
	as_n_m <- mean( rowSums(asd4[,paste0('earlystop',1:5)])==8 )

	rr_a_m <- mean( unlist(asd4[,paste0('reject',6:10)]) )
	es_a_m <- mean( unlist(asd4[,paste0('N',6:10)]) )
	sr_a_m <- mean( unlist(asd4[,paste0('earlystop',6:10)]) )
	as_a_m <- mean( rowSums(asd4[,paste0('earlystop',6:10)])==2 )

	mat_res['mixed_5null_5alt',] <- c(rr_n_m,fw_n_m,es_n_m,sr_n_m,as_n_m,rr_a_m,es_a_m,sr_a_m,as_a_m)

	return(mat_res)
}


simon16_sum <- function(res, nullp=0.1, altp=0.3, scen, postp){
### Function to create table summarizing results for Simon 2016 results
# res: file with results
# nullp: response rate in null baskets (to be used in identifying match to null_ind)
# altp: response rate in alternative baskets (to be used in identifying match to alt_ind)
# scen: vector of 0.1 or 0.3 values for each basket (use all 0.1 for global null/alternative)
# postp: posterior probability threshold used to declare efficacy

	#Extract null and alternative data sets (written 
	scen_null <- paste0(scen,collapse='/')
	null_basket <- which(scen == nullp)

	scen_alt <- if( sum(scen==nullp) == length(scen) ){ paste0(rep(altp,length(scen)), collapse='/') }else{ scen_null }
	alt_basket <- if( sum(scen==nullp) == length(scen) ){ null_basket }else{ which(scen == altp) }
	
	asd0n <- res[which(res$postp==postp & res$scen==scen_null & res$UB=='simon'),]
	asd0a <- res[which(res$postp==postp & res$scen==scen_alt & res$UB=='simon'),]

	#Create object to store results
	col_names <- c('rr_n','fw_n','es_n','sr_n','as_n','rr_a','es_a','sr_a','as_a')
	mat_res <- matrix( nrow=11, ncol=9, dimnames=list(seq(0,0.5,by=0.05),col_names)) 

	for( i in c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)){ 
		rr_n <- mean( unlist(asd0n[which(asd0n$prep==i),paste0('reject',null_basket)]) )
		fw_n <- mean( rowSums(asd0n[which(asd0n$prep==i),paste0('reject',null_basket)]) > 0 )
		es_n <- mean( unlist(asd0n[which(asd0n$prep==i),paste0('N',null_basket)]) )
		sr_n <- mean( unlist(asd0n[which(asd0n$prep==i),paste0('earlystop',null_basket)]) )
		as_n <- mean( rowSums(asd0n[which(asd0n$prep==i),paste0('earlystop',null_basket)])==length(null_basket) )

		rr_a <- mean( unlist(asd0a[which(asd0a$prep==i),paste0('reject',alt_basket)]) )
		es_a <- mean( unlist(asd0a[which(asd0a$prep==i),paste0('N',alt_basket)]) )
		sr_a <- mean( unlist(asd0a[which(asd0a$prep==i),paste0('earlystop',alt_basket)]) )
		as_a <- mean( rowSums(asd0a[which(asd0a$prep==i),paste0('earlystop',alt_basket)])==length(alt_basket) )

		mat_res[paste0(i),] <- c(rr_n,fw_n,es_n,sr_n,as_n,rr_a,es_a,sr_a,as_a)
	}

	return(mat_res)

}

### Summarize the results with above functions

# MEM approaches
marg_global <- sim_sum(res=dat, nullp=0.1, altp=0.3, scen=rep(0.1,10), postp0=0.9, postp1=0.848)
fw_global <- sim_sum(res=dat, nullp=0.1, altp=0.3, scen=rep(0.1,10), postp0=0.99, postp1=0.979)

marg_mix <- sim_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,8),rep(0.3,2)), postp0=0.9, postp1=0.848)
fw_mix <- sim_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,8),rep(0.3,2)), postp0=0.99, postp1=0.979)

marg_mix5 <- sim_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,5),rep(0.3,5)), postp0=0.9, postp1=0.848)
fw_mix5 <- sim_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,5),rep(0.3,5)), postp0=0.99, postp1=0.979)

# Simon 2-stage
simon <- simon_sum(res=dat)

# Simon 2016 Posterior Probability Monitoring Approach
marg16_global <- simon16_sum(res=dat, nullp=0.1, altp=0.3, scen=rep(0.1,10), postp=0.154)
fw16_global <- simon16_sum(res=dat, nullp=0.1, altp=0.3, scen=rep(0.1,10), postp=0.665)

marg16_mix <- simon16_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,8),rep(0.3,2)), postp=0.154)
fw16_mix <- simon16_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,8),rep(0.3,2)), postp=0.665)

marg16_mix5 <- simon16_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,5),rep(0.3,5)), postp=0.154)
fw16_mix5 <- simon16_sum(res=dat, nullp=0.1, altp=0.3, scen=c(rep(0.1,5),rep(0.3,5)), postp=0.665)



######################################
### TeX tables

xtable(marg_global[,1:10], digits=c(0,3,3,3,3,1,1,3,3,3,3), align='l|cc|cc|cc|cc|cc')
xtable(marg_global[,11:18], digits=c(0,3,3,1,1,3,3,3,3), align='l|cc|cc|cc|cc')

xtable(marg_global[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')
xtable(fw_global[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')
xtable(marg_mix[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')
xtable(fw_mix[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')

xtable(simon[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')


######################################
### TeX tables

xtable(marg_global[,1:10], digits=c(0,3,3,3,3,1,1,3,3,3,3), align='l|cc|cc|cc|cc|cc')
xtable(marg_global[,11:18], digits=c(0,3,3,1,1,3,3,3,3), align='l|cc|cc|cc|cc')

xtable(marg_global[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')
xtable(marg_mix[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')
xtable(marg_mix5[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')

xtable(fw_global[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')
xtable(fw_mix[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')
xtable(fw_mix5[,1:16], digits=c(0,3,3,3,3,1,1,3,3,3,3,3,3,1,1,3,3), align='l|cc|cc|cc|cc|cc|cc|cc|cc')

xtable(marg16_global[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')
xtable(marg16_mix[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')
xtable(marg16_mix5[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')

xtable(fw16_global[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')
xtable(fw16_mix[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')
xtable(fw16_mix5[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')

xtable(simon[,1:8], digits=c(0,3,3,1,3,3,3,1,3), align='l|c|c|c|c|c|c|c|c')


######################################
### Plot Results

## Select scenario and type I error rate control to plot (i.e., global (all null or all alternative), mixed with 8 null/2 alternative, mixed with 5 null/5 alternative; marginal (basket-wise) or family-wise type I error rate control)

# Global, Marginal
res_use <- marg_global; simon_use <- 'global'; s16_use <- marg16_global

# Global, Familywise
res_use <- fw_global; simon_use <- 'global'; s16_use <- fw16_global


# Mixed 8 null/2 alt, Marginal
res_use <- marg_mix; simon_use <- 'mixed_8null_2alt'; s16_use <- marg16_mix

# Mixed 8 null/2 alt, Familywise
res_use <- fw_mix; simon_use <- 'mixed_8null_2alt'; s16_use <- fw16_mix


# Mixed 5 null/5 alt, Marginal
res_use <- marg_mix5; simon_use <- 'mixed_5null_5alt'; s16_use <- marg16_mix5

# Mixed 5 null/5 alt, Familywise
res_use <- fw_mix5; simon_use <- 'mixed_5null_5alt'; s16_use <- fw16_mix5



### Panel set
par(mfrow=c(2,2))

### Rejection rate
plot( x=-10, y=-10, xlab='Probability Threshold', ylab='Marginal Rejection Rate', xlim=c(0,0.5), ylim=c(0,1) )
title('Rejection Rate')

abline(h=simon[simon_use,'rr_n'], lty=3, lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'rr_n_1'], lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'rr_n_0'], lty=2, lwd=2)
lines(x=seq(0,0.5,by=0.05), y=s16_use[,'rr_n'], lty=4, lwd=2, col='black')

abline(h=simon[simon_use,'rr_a'], lty=3, col='gray65', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'rr_a_1'], col='gray65', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'rr_a_0'], lty=2, col='gray85', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=s16_use[,'rr_a'], lty=4, lwd=2, col='gray65')

### Expected size
plot( x=-10, y=-10, xlab='Probability Threshold', ylab='Expected Size', xlim=c(0,0.5), ylim=c(0,25) )
title('Expected Size')

abline(h=simon[simon_use,'es_n'], lty=3, lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'es_n_1'], lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'es_n_0'], lty=2, lwd=2)
lines(x=seq(0,0.5,by=0.05), y=s16_use[,'es_n'], lty=4, lwd=2, col='black')

abline(h=simon[simon_use,'es_a'], lty=3, col='gray65', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'es_a_1'], col='gray65', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'es_a_1'], lty=2, col='gray85', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=s16_use[,'es_a'], lty=4, lwd=2, col='gray65')

### Stop rates
plot( x=-10, y=-10, xlab='Probability Threshold', ylab='Stop Rate', xlim=c(0,0.5), ylim=c(0,1) )
title('Stop Rates')

abline(h=simon[simon_use,'sr_n'], lty=3, lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'sr_n_1'], lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'sr_n_0'], lty=2, lwd=2)
lines(x=seq(0,0.5,by=0.05), y=s16_use[,'sr_n'], lty=4, lwd=2, col='black')

abline(h=simon[simon_use,'sr_a'], lty=3, col='gray65', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'sr_a_1'], col='gray65', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=res_use[,'sr_a_0'], lty=2, col='gray85', lwd=2)
lines(x=seq(0,0.5,by=0.05), y=s16_use[,'sr_a'], lty=4, lwd=2, col='gray65')

### Legend
plot( x=-10, y=-10, xlim=c(0,1), ylim=c(0,1), xlab='', ylab='', xaxt='n', yaxt='n', frame.plot=F)
leg_text <- c('Null Baskets, Information Sharing',
	'Alternative Baskets, Information Sharing',
	'Null Baskets, Bayesian','Alternative Baskets, Bayesian',
	'Null Baskets, Simon 2016','Alternative Baskets, Simon 2016',
	'Null Baskets, Simon Two-Stage','Alternative Baskets, Simon Two-Stage')
legend('left', inset=-0.25, lty=c(1,1,2,2,4,4,3,3), col=rep(c('black','gray65','black','gray85','black','gray65','black','gray65')), lwd=c(2,2,2,2,2,2,2,2), legend=leg_text, bty='n', cex=1, xpd=T)




