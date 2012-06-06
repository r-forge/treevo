doSimulation<-function(splits, intrinsicFn, extrinsicFn, startingStates, intrinsicValues, extrinsicValues, timeStep, saveHistory=FALSE, saveRealParams=FALSE, jobName="") {

  taxIds<-function(i){
    if (id(taxa[[i]])==splits[1, 2]){
      taxa<-c(taxa, taxa[[i]], taxa[[i]])
      id(taxa[[originallength+1]])<-splits[1, 3]
      timeSinceSpeciation(taxa[[originallength+1]])<-0
      id(taxa[[originallength+2]])<-splits[1, 4]
      timeSinceSpeciation(taxa[[originallength+2]])<-0 #this is still confusing to me!
      taxontodelete<<-i
      taxa
    }
  }



  otherstatefn<-function(x){
    states(taxa[[x]])
  }

  otherMatrix<-function(i){
    taxvec<-c(1:length(taxa))
    taxvec<-taxvec[-which(taxvec==i)]
    otherstatesvector<-sapply(taxvec,otherstatefn)
    otherstatesmatrix<-matrix(otherstatesvector, ncol=length(states(taxa[[i]])), byrow=TRUE)
    newvalues<-states(taxa[[i]])+intrinsicFn(params=intrinsicValues, states=states(taxa[[i]]), timefrompresent =timefrompresent)+extrinsicFn(params=extrinsicValues, selfstates=states(taxa[[i]]), otherstates=otherstatesmatrix, timefrompresent =timefrompresent)
    nextstates(taxa[[i]])<-newvalues
    return(taxa[[i]])
  }
  stateNextState<-function(i){
    states(i)<-nextstates(i)
    return(i)
  }

  timeSinceSp<-function(i){
    timeSinceSpeciation(i)<-timeSinceSpeciation(i)+timeStep
    return(i)
  }

  if (saveRealParams){
    RealParams<-vector("list", 2)
    names(RealParams)<-c("matrix", "vector")	
    RealParams$vector<-c(startingStates, intrinsicValues, extrinsicValues)
    maxLength<-(max(length(startingStates), length(intrinsicValues), length(extrinsicValues)))
    RealParams$matrix<-matrix(ncol=maxLength, nrow=3)
    rownames(RealParams$matrix)<-c("startingStates", "intrinsicFn", "extrinsicFn")
    RealParams$matrix[1,]<-c(startingStates, rep(NA, maxLength-length(startingStates)))
    RealParams$matrix[2,]<-c(intrinsicValues, rep(NA, maxLength-length(intrinsicValues)))
    RealParams$matrix[3,]<-c(extrinsicValues, rep(NA, maxLength-length(extrinsicValues)))
    save(RealParams, file=paste("RealParams", jobName, ".Rdata", sep=""))
  }
  if (saveHistory) {
    startVector<-c()
    endVector<-c()
    startTime<-c()
    endTime<-c()
  }
	numberofsteps<-floor(splits[1, 1]/timeStep)
	mininterval<-min(splits[1:(dim(splits)[1]-1), 1]-splits[2:(dim(splits)[1]), 1])
	if (numberofsteps<1000) {
		warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", splits[1, 1]/1000))
	}
	if (floor(mininterval/timeStep)<50) {
		warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more if you expect change on this branch. I would suggest decreasing timeStep to no more than ", mininterval/50))
	}
#initial setup
	timefrompresent=splits[1, 1]
	taxa<-c(new("abctaxon", id=splits[1, 3], states=startingStates), new("abctaxon", id=splits[1, 4], states=startingStates))
	splits<-splits[2:dim(splits)[1], ] #pop off top value	
#start running
	while(timefrompresent>0) {

		while ((timefrompresent-timeStep)<=splits[1, 1]) { #do speciation. changed from if to while to deal with effectively polytomies
			originallength<-length(taxa)
                        idvec<-c(1:originallength)
                        foo<-lapply(idvec,taxIds)
                        foo<-foo[!sapply(foo,is.null)]
                        taxa<-foo[[1]]
                                
			taxa<-taxa[-taxontodelete]
			if(dim(splits)[1]>1) {
				splits<-splits[2:(dim(splits)[1]), ] #pop off top value
			}
			else {
				splits[1, ]<-c(-1, 0, 0, 0)
			}
                       
                      }



                taxavec<-c(1:length(taxa))
                taxa<-lapply(taxavec,otherMatrix)               
                taxa<-lapply(taxa,stateNextState)	
		timefrompresent<-timefrompresent-timeStep
                taxa<-lapply(taxa,timeSinceSp)

              
              }

	summarizeTaxonStates(taxa)
}
