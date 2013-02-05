set.seed(1234)

#simsize <- 1000
n <- 300
err.tol <- 0.001
#N <- n*simsize


# Tianxiao: the following 4 examples are included in the paper

# Ex 1 (psi0 = 0 = psi1):
# gamma <- c(0,0,0,0,0,0,0)
# delta <- c(.5,.5)

#Ex 2 (psi0 = 0.5, psi1 = 0):
# gamma <- c(0,0,0,0,0.5,0,0.5)
# delta <- c(.5,.5)

#Ex 3 (psi0 = 1, psi1 = 0.5):
# gamma <- c(0,0,0.5,0.5,1,0.5,0.5)
# delta <- c(1,0)

#Ex 4 (psi0 = 1, psi1 = 1):
# gamma <- c(0,0,0.5,0.77,1,1,1)
# delta <- c(1,0)

expit <- function(x) { exp(x)/(1+exp(x)) }


# Tianxiao: the following 'para.generation' function is for generating some additional scenarios, not included in the original paper I sent you

para.generation <- function(){
gamma<-numeric(7)
f<-k<-numeric(4)
delta<-runif(2)
delta<-round(delta,digits=4)
cat("delta = ", delta, "\n")
k[1]<-0.25*expit(delta[1]+delta[2])
k[2]<-0.25*expit(-delta[1]+delta[2])
k[3]<-0.25*expit(delta[1]-delta[2])
k[4]<-0.25*expit(-delta[1]-delta[2])
k<-round(k,digits=4)
cat("k = ", k, "\n")
gamma[5:7]<- rnorm(3)
f[1]<-gamma[5]+gamma[6]+gamma[7]
f[2]<-gamma[5]+gamma[6]-gamma[7]
f[3]<-gamma[5]-gamma[6]+gamma[7]
f[4]<-gamma[5]-gamma[6]-gamma[7]
f<-round(f,digits=4)
cat("f = ", f, "\n")
gamma[3]<-gamma[5]-(k[1]+k[2])*(abs(f[1])-abs(f[4]))+(k[3]+k[4])*(abs(f[2])-abs(f[3]))
gamma[4]<-gamma[6]-(k[1]-k[2])*(abs(f[1])-abs(f[3]))+(k[3]-k[4])*(abs(f[2])-abs(f[4]))
gamma<-round(gamma,digits=4)
cat("gamma = ", gamma, "\n")
gamma.delta<-c(gamma,delta)
return (gamma.delta)
}


# Tianxiao: the following function 'data.generation' is for generating data (O's, A's and Y) based on the above parameter setting
 
data.generation <- function(N, gamma, delta){ 
	
	full.data <- NULL
      O <- matrix(0,N,2)
	temp <- 2*matrix(rbinom(N*3, 1, 0.5),N,3)-1
	A <- temp[,c(1:2)]
	O[,1] <- temp[,3]
	p <- expit(delta[1]*O[,1]+delta[2]*A[,1])	
	O[,2] <- 2*matrix(rbinom(N,1,p),N,1)-1
	
	
	mu <- gamma[1] + gamma[2]*O[,1] + gamma[3]*A[,1] + gamma[4]*O[,1]*A[,1] + gamma[5]*A[,2] + gamma[6]*O[,2]*A[,2] + gamma[7]*A[,1]*A[,2]
	Y <- rnorm(N,0,1)
	
	full.data$Y <- Y + mu
	full.data$A <- A
	full.data$O <- O
	
	return(full.data)
}


# Tianxiao: the following 'calc.pseudooutcome' is a sub-routine to calculate the (hard-max) pseudo-outcome that is to be used in the main program
 
calc.pseudooutcome <- function(sim.data,Beta.Psi=NULL){
	
	Y <- sim.data$Y 
	A <- sim.data$A 
	O <- sim.data$O	
		
	beta.psi.2 <- NULL
	if (!is.null(Beta.Psi)){beta.psi.2 <- Beta.Psi[1:7]}
	else{
		X2 <- cbind(1, O[,1], A[,1], O[,1]*A[,1], A[,2], O[,2]*A[,2], A[,1]*A[,2])
		beta.psi.2 <- as.numeric(coef(lm(Y~X2-1)))
		}
	
	Yprimeprime <- beta.psi.2[5] + beta.psi.2[6]*O[,2] + beta.psi.2[7]*A[,1]
	
	stage2dim <- length(beta.psi.2)

	Sigma2 <- n*solve(t(X2)%*%X2)
	
	Z2 <- (diag(array(Y - X2%*%matrix(beta.psi.2,stage2dim,1)))%*%X2%*%Sigma2)/sqrt(n-7)

    Cov2 <- t(Z2)%*%Z2
    Sigma.stage2 = Cov2[c(5:7),c(5:7)]
    h<-cbind(1, O[,2], A[,1])
    k<-diag(h%*%Cov2[c(5:7),c(5:7)]%*%t(h))
    n <- length(Y)
    lambda <- 3*k/n
    soft.th <- (1-lambda/Yprimeprime/Yprimeprime)
    soft.factor <- soft.th*(soft.th > 0)
	
	Y.tilde <- beta.psi.2[1] +  beta.psi.2[2]*O[,1] +  beta.psi.2[3]*A[,1] +  beta.psi.2[4]*A[,1]*O[,1] + abs(Yprimeprime)*soft.factor
	
    return(Y.tilde)	
	}



# Tianxiao: the following 'q.hard.shared.est' is the shared hardmax procedure (uses 'qr.solve' function)

q.hard.shared.est <- function(sim.data, err.tol, maximum.iter=500){
	
	# array of psi values at each iteration
	results <- list()
	
	Y <- sim.data$Y 
	A <- sim.data$A 
	O <- sim.data$O
	
	Y.tilde <- calc.pseudooutcome(sim.data,Beta.Psi=NULL)
	
	Y.star <- c(Y, Y.tilde)
	Z <- rbind(cbind(1, O[,1], A[,1], O[,1]*A[,1], A[,2], O[,2]*A[,2], A[,1]*A[,2], 0, 0), cbind(0, 0, 0, 0, A[,1], A[,1]*O[,1], 0, 1, O[,1]))  # initial value of the matrix
	
     	Beta.Psi <- qr.solve(Z,Y.star)
	Y.star <- matrix(Y.star,,1)
	
	output.vec <- matrix(Beta.Psi[c(5:6)],1,2)
	
	output.vec.prev <- output.vec*1000  # just for initialization
	total.iter <- 1
	
	results$Psi0[[total.iter]] <- output.vec[1]
	results$Psi1[[total.iter]] <- output.vec[2]
	
	while ((max(abs(output.vec-output.vec.prev)) > err.tol) & (total.iter<=maximum.iter)){
		total.iter <- total.iter + 1
		output.vec.prev <- output.vec
		Y.tilde <- calc.pseudooutcome(sim.data,Beta.Psi)
		Y.star <- matrix(c(Y, Y.tilde),,1)
		Beta.Psi <- qr.solve(Z,Y.star)
		output.vec <- matrix(Beta.Psi[c(5:6)],1,2)
		
		results$Psi0[[total.iter]] <- output.vec[1]
	    results$Psi1[[total.iter]] <- output.vec[2]
		}
	if (total.iter > maximum.iter) {print('Maximum iterations reached in q.hard.shared')}
	
	results$total.iter <- total.iter
	results$output.vec <- output.vec
	
	return(results)
}


# Tianxiao: the following 'nonsmooth.min.est' is a variation of the above 'q.hard.shared.est' procedure (uses 'nlminb' function instead of 'qr.solve') -- the goal is to investigate if these various solving/optimization procedures make any difference

nonsmooth.min.est <- function(sim.data, err.tol, maximum.iter=500){
	
	# array of psi values at each iteration
	results <- list()
	
	Y <- sim.data$Y 
	A <- sim.data$A 
	O <- sim.data$O
	
	Y.tilde <- calc.pseudooutcome(sim.data,Beta.Psi=NULL)
	
	Y.star <- c(Y, Y.tilde)
	Z <- rbind(cbind(1, O[,1], A[,1], O[,1]*A[,1], A[,2], O[,2]*A[,2], A[,1]*A[,2], 0, 0), cbind(0, 0, 0, 0, A[,1], A[,1]*O[,1], 0, 1, O[,1]))  # initial value of the matrix
	
      residual<-function(theta){
      r<-Y.star-Z%*%theta
      t(r)%*%r
      }      

	Beta.Psi <- nlminb(c(0,0,0,0,0,0,0,0,0),residual)$par
	Y.star <- matrix(Y.star,,1)
	
	output.vec <- matrix(Beta.Psi[c(5:6)],1,2)
	
	output.vec.prev <- output.vec*1000  # just for initialization
	total.iter <- 1
	
	results$Psi0[[total.iter]] <- output.vec[1]
	results$Psi1[[total.iter]] <- output.vec[2]
	
	while ((max(abs(output.vec-output.vec.prev)) > err.tol) & (total.iter<=maximum.iter)){
		total.iter <- total.iter + 1
		output.vec.prev <- output.vec
		Y.tilde <- calc.pseudooutcome(sim.data,Beta.Psi)
		Y.star <- matrix(c(Y, Y.tilde),,1)
		Beta.Psi <- nlminb(Beta.Psi,residual)$par
		output.vec <- matrix(Beta.Psi[c(5:6)],1,2)
		
		results$Psi0[[total.iter]] <- output.vec[1]
	    results$Psi1[[total.iter]] <- output.vec[2]
		}
	if (total.iter > maximum.iter) {print('Maximum iterations reached in q.hard.shared')}
	
	results$total.iter <- total.iter
	results$output.vec <- output.vec
	
	return(results)
}


# Tianxiao: the following 'nonlinear.eqa.est' is a third variation (uses 'multiroot' function) -- the goal is to investigate if these various solving/optimization procedures make any difference

library(rootSolve)
nonlinear.eqa.est <- function(sim.data, err.tol, maximum.iter=500){
     
	# array of psi values at each iteration
	results <- list()
	
	Y <- sim.data$Y 
	A <- sim.data$A 
	O <- sim.data$O
	
	Y.tilde <- calc.pseudooutcome(sim.data,Beta.Psi=NULL)
	
	Y.star <- c(Y, Y.tilde)
	Z <- rbind(cbind(1, O[,1], A[,1], O[,1]*A[,1], A[,2], O[,2]*A[,2], A[,1]*A[,2], 0, 0), cbind(0, 0, 0, 0, A[,1], A[,1]*O[,1], 0, 1, O[,1]))  # initial value of the matrix
	
      equation<-function(theta){
      r<-Y.star-Z%*%theta
      t(Z)%*%r
      }      

	Beta.Psi <- multiroot(equation,start=c(0,0,0,0,0,0,0,0,0))$root
	Y.star <- matrix(Y.star,,1)
	
	output.vec <- matrix(Beta.Psi[c(5:6)],1,2)
	
	output.vec.prev <- output.vec*1000  # just for initialization
	total.iter <- 1
	
	results$Psi0[[total.iter]] <- output.vec[1]
	results$Psi1[[total.iter]] <- output.vec[2]
	
	while ((max(abs(output.vec-output.vec.prev)) > err.tol) & (total.iter<=maximum.iter)){
		total.iter <- total.iter + 1
		output.vec.prev <- output.vec
		Y.tilde <- calc.pseudooutcome(sim.data,Beta.Psi)
		Y.star <- matrix(c(Y, Y.tilde),,1)
		Beta.Psi <- multiroot(equation,start=Beta.Psi)$root
		output.vec <- matrix(Beta.Psi[c(5:6)],1,2)
		
		results$Psi0[[total.iter]] <- output.vec[1]
	    results$Psi1[[total.iter]] <- output.vec[2]
		}
	if (total.iter > maximum.iter) {print('Maximum iterations reached in q.hard.shared')}
	
	results$total.iter <- total.iter
	results$output.vec <- output.vec
	
	return(results)
}




# Actually run the function
para<-para.generation()
gamma<-para[1:7]
delta<-para[8:9]
sim.data <- data.generation(n, gamma, delta)
results.qr <- q.hard.shared.est(sim.data, err.tol, maximum.iter=500)
results.nm <- nonsmooth.min.est(sim.data, err.tol, maximum.iter=500)
results.ne <- nonlinear.eqa.est(sim.data, err.tol, maximum.iter=500)

par(mfrow=c(1,3))
plot(results.qr$Psi0,type="b",ylim=c(min(results.qr$Psi0,results.qr$Psi1,gamma[5:6]),max(results.qr$Psi0,results.qr$Psi1,gamma[5:6])),ylab="Example 14 qr.solve")
points(results.qr$Psi1,type="b",lty=2)
legend("right",legend=c("Psi0","Psi1"),lty=1:2)
abline(h=gamma[5])
abline(h=gamma[6],lty=2)
plot(results.nm$Psi0,type="b",ylim=c(min(results.nm$Psi0,results.nm$Psi1,gamma[5:6]),max(results.nm$Psi0,results.nm$Psi1,gamma[5:6])),ylab="Example 14 nlminb")
points(results.nm$Psi1,type="b",lty=2)
legend("right",legend=c("Psi0","Psi1"),lty=1:2)
abline(h=gamma[5])
abline(h=gamma[6],lty=2)
plot(results.ne$Psi0,type="b",ylim=c(min(results.ne$Psi0,results.ne$Psi1,gamma[5:6]),max(results.ne$Psi0,results.ne$Psi1,gamma[5:6])),ylab="Example 14 multiroot")
points(results.ne$Psi1,type="b",lty=2)
legend("right",legend=c("Psi0","Psi1"),lty=1:2)
abline(h=gamma[5])
abline(h=gamma[6],lty=2)

results.qr
results.nm
results.ne



# Other functions to use in the iterative scheme:

# "nlminb" for nonlinear optimization
# replace the 'max' operator by sigmoid approximation
# replace the 'max' operator by soft-threshold approximation
