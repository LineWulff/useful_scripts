## Function converting a list of vectors of varying length
## to list of vectors of the same length - inserting NAs
## useful in creating csvs of uneven lists

#list needs indivdual names to worj, NOT numbers

vectorlist <- function(unevenlist){
  maxlen <- 0
  for (clus in names(unevenlist)){
    if (maxlen<length(unevenlist[[clus]])){
      maxlen <- length(unevenlist[[clus]])}}
  #print(maxlen)
  evenlist <- list()
  for (clus in names(unevenlist)){
    diffb <- maxlen-length(unevenlist[[clus]])
    #print(paste(clus,maxlen,diffb))
    evenlist[[clus]] <- c(unevenlist[[clus]],rep(NA,times=as.numeric(diffb)))}
  return(evenlist)
}
