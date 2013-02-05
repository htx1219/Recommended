all<-matrix(0, 12000, 12)
for (i in 1:10){
	all[,i]<-rnorm(12000)
}
all[,11]<-all[,1]^2+all[,2]^2+all[,3]^2+all[,4]^2+all[,5]^2+all[,6]^2+all[,7]^2+all[,8]^2+all[,9]^2+all[,10]^2
all[,12]<-2*(all[,11]>qchisq(0.5,10))-1
all<-data.frame(all)
train<-all[1:2000,]
test<-all[2001:12000,]

library(ada)
library(rpart)

maxIter = 400

adahtx<-function(train, maxIter, print = 0, ...){
N = dim(train)[1]
weight = rep(1/N, N)
trees = list()
weights<-matrix(0, maxIter, N)
err<-rep(0, maxIter)
alpha <- rep(0, maxIter)
for (i in 1:maxIter){
	if (print > 0){
		if (i %% print  == 0){
			print(i)
		}
	}
	weights[i, ]<-weight
	trees[[i]]<-rpart(X12~.-X11-X12, data=train, weights = weight, method = "class", ...) #, maxdepth=1,cp=-1,minsplit=0,xval=0)
	pred <- predict(trees[[i]])[,2]
	errcase = 1-(train$X12 ==(2*(pred>0.5)-1))
	err[i] <- sum(weight*errcase)/sum(weight)
	alpha[i]<-log((1-err[i])/err[i])
	weight = weight*exp(alpha[i]*errcase)
}
htxada <- list()
htxada$trees <- trees
htxada$alpha <- alpha
return(htxada)
}

adapredhtx<-function(htxada, iter, data){
	trees <- htxada$trees
	alpha <- htxada$alpha
	allpred = matrix(0, iter, dim(data)[1])
	finalpred = rep(0, dim(data)[1])
	for (i in 1:iter){
		pred = predict(trees[[i]], data)[,2]
		pred <- 2*(pred>0.5)-1
		finalpred = finalpred + alpha[i]*pred
		allpred[i, ]<-2*(finalpred>0)-1
	}
	return(allpred)
}

checkerr<-function(allpred, pred){
	iter = dim(allpred)[1]
	err <- rep(0, iter)
	for (i in 1:iter){
		err[i] = 1-mean(pred == allpred[i, ])
	}
	return(err)
}

# Adaboost Using Tree Iter=200
htxada1<-adahtx(train, 200, 20)
allpred1 = adapredhtx(htxada1, 200, train)
err1<-checkerr(allpred1, train$X12)

allpred2 = adapredhtx(htxada1, 200, test)
err2<-checkerr(allpred2, test$X12)
plot(c(0,200), c(0, 0.35), ylab="Error", xlab="Iteration 1 to 200", main = "Training And Testing Error", type='n')
lines(err1, type = 'l', col=2)
lines(err2, type = 'l', col=4)
legend(x="topright", legend=c("Training Error", "Test Error"), lty=1, col=c(2,4), lwd=2)

# Adaboost Using Tree Iter=1000
htxada2<-adahtx(train, 1000, 50)
allpred2.1 = adapredhtx(htxada2, 1000, train)
err2.1<-checkerr(allpred2.1, train$X12)

allpred2.2 = adapredhtx(htxada2, 1000, test)
err2.2<-checkerr(allpred2.2, test$X12)
plot(c(0,1000), c(0, 0.35), ylab="Error", xlab="Iteration 1 to 1000", main = "Training And Testing Error", type='n')
lines(err2.1, type = 'l', col=2)
lines(err2.2, type = 'l', col=4)
legend(x="topright", legend=c("Training Error", "Test Error"), lty=1, col=c(2,4), lwd=2)

# Adaboost Using Stumps Iter=400
htxada3<-adahtx(train, 400, 40, maxdepth=1,cp=-1,minsplit=0,xval=0)
allpred3.1 = adapredhtx(htxada3, 400, train)
err3.1<-checkerr(allpred3.1, train$X12)

allpred3.2 = adapredhtx(htxada3, 400, test)
err3.2<-checkerr(allpred3.2, test$X12)
plot(c(0,400), c(0, 0.35), ylab="Error", xlab="Iteration 1 to 400", main = "Training And Testing Error\n Using Stumps", type='n')
lines(err3.1, type = 'l', col=2)
lines(err3.2, type = 'l', col=4)
legend(x="topright", legend=c("Training Error", "Test Error"), lty=1, col=c(2,4), lwd=2)

# Question 4
part1<-all[1:6000,]
part1[, 12]<--1
all2<-matrix(0, 30000, 12)
for (i in 1:10){
	all2[,i]<-rnorm(30000)
}
all2[,11]<-all2[,1]^2+all2[,2]^2+all2[,3]^2+all2[,4]^2+all2[,5]^2+all2[,6]^2+all2[,7]^2+all2[,8]^2+all2[,9]^2+all2[,10]^2
all2[,12]<-1
all2<-data.frame(all2)
part2<-all2[all2$X11>12,]
part2.r<-part2[1:6000,]
d<-rbind(part1, part2.r)
d<-d[sample(12000, 12000),]
train2<-d[1:2000,]
test2<-d[2001:12000,]

htxada4<-adahtx(train2, 200, 20)
allpred4.1 = adapredhtx(htxada4, 200, train2)
err4.1<-checkerr(allpred4.1, train2$X12)

allpred4.2 = adapredhtx(htxada4, 200, test2)
err4.2<-checkerr(allpred4.2, test2$X12)
plot(c(0,200), c(0, 0.35), ylab="Error", xlab="Iteration 1 to 200", main = "Training And Testing Error", type='n')
lines(err4.1, type = 'l', col=2)
lines(err4.2, type = 'l', col=4)
legend(x="topright", legend=c("Training Error", "Test Error"), lty=1, col=c(2,4), lwd=2)

# from package "ada"
ada1<-ada(train[,1:10], train[,12], text.x<-test[,1:10], test.y<-test[,12], loss="exponential", iter=200, bag.frac=1, nu=1)
plot(c(0,200), c(0, 0.35), ylab="Error", xlab="Iteration 1 to 200", main = "Training And Testing Error", type='n')
lines(ada1$model$errs[,1], col='red')
lines(ada1$model$errs[,3], col='blue')
legend(x="topright", legend=c("Training Error", "Test Error"), lty=1, col=c(2,4), lwd=2)

ada3<-ada(train2[,1:10], train2[,12], text.x<-test2[,1:10], test.y<-test2[,12], loss="exponential", iter=200, bag.frac=1, nu=1)
plot(c(0,200), c(0, 0.35), ylab="Error", xlab="Iteration 1 to 200", main = "Training And Testing Error", type='n')
lines(ada3$model$errs[,1], col='red')
lines(ada3$model$errs[,3], col='blue')