setwd('C:/Users/topalw/Desktop/')
dat <- read.csv('contigs.lengths')


# the logic behind it 
plot(cumsum(dat$Length))
abline(h=(sum(dat$Length)/10*1:length(dat$Length)))


### FUNCTION ###
divider <- function(dist, chunks=10,reps=100){
  ##########################################
  # function divides a vector of values into k parts of ~ similar sum value
  # this is a stupid logic that works and is only a bit optimized by a randomization 
  # step and a minimize the variance of sums of all groups 
  # However there are still errors ex an equidistant point will fail or there might be 
  # more optimal ways of dividing but for rough estimation this will do.
  # also it is quite slow in reps > 1000 
  ##########################################
  # vector for sum of square deviations from target value
  ssr <- rep(0,reps)
  # make mat to save grouping
  group.mat <- matrix(data=NA,nrow=length(dist),ncol=reps)
  for(srep in 1:reps){
  # shuffle the order of the elements
  rep.order <- sample(1:length(dist), length(dist))
  dist <- dist[rep.order]
  # make cum addition
  cum.add <- cumsum(dist)
  # define targets and make result vector
  targets <- sum(dist)/chunks * seq(1,chunks)
  result <- rep(0,length(dist))
  # loop through targets and assign groups
  last <- 0
  for(target in targets){
    distances <- abs(cum.add - target)
    group.name <- which(targets == target)
    # upper bound 
    upper.bound <- which(distances == min(distances))
    lower.bound <- last + 1 
    result[lower.bound:upper.bound] <- group.name
    last <- upper.bound
  }
  # save grouping by re-ordering to match original
  group.mat[,srep] <- result[order(rep.order)]
  # save SSR of each rep
  ssr[srep] <-var(aggregate(dist, by=list(result), FUN=sum)[,2])
}
  # after reps finish 
  # choose best ssr 
  best <- which(ssr == min(ssr))[1] # might be > 1 
  # return best order 
  return(ssr)
  #return(group.mat[,best])
}


### Checking the function ###  ran with return(ssr)
#plot(ssr,ylab='variance between groups',xlab='replicate of random order')
#hist(ssr,xlab='variance between groups / replicate')

### What variances do different # of groups give? ###
#min.var <- rep(0,9)
#for(k in 2:10){
#min.var[k-1] <- min(divider(dat$Length,k,1000))[1]
#}
#plot(y=min.var, x=seq(2,10)) # obviously larger and fewer groups have smaller variance

### does the function find the smallest variance split consistently? ###
#min.var <- rep(0,100)
#for(k in 1:100){
#  min.var[k] <- min(divider(dat$Length,10,1000))[1]
#}
#plot(y=min.var, x=seq(1,100))
#hist(min.var) # range 1-5e+13 for 100 x 1000 reps

### does the function depend a lot on the number of replicates? ###
# min.var <- rep(0,length(seq(100,10000,by=100)))
# for(k in seq(100,10000,by=100)){
#   min.var[k/100] <- min(divider(dat$Length,10,k))[1]
# }
# plot(y=min.var, x=seq(100,10000,by=100))
# abline(lm(min.var~seq(100,10000,by=100)))
# summary(lm(min.var~seq(100,10000,by=100))) # -1.295e+09  2.284e+08  -5.671 1.43e-07 ***


### Running the function ###  ran with return(group.mat[,best])
dat$groups <- divider(dat$Length,10,1000)
# plot
plot(aggregate(dat$Length, by=list(dat$groups), FUN=sum)[,2])
abline(h=sum(dat$Length)/10)
#write.csv(dat,'conting_lengths_grouped.csv',quote=F, row.names = F)
