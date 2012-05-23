library(TreEvo)
phy<-rcoal(200)
splits<-getSimulationSplits(phy)
phy$edge.length<-phy$edge.length/max(branching.times(phy))
taxframe<-doSimulation(
		splits=splits, 
		intrinsicFn= brownianIntrinsic, #simple Brownian model
		extrinsicFn= nullExtrinsic, #null
		startingStates=c(10), #root state
		intrinsicValues=c(0.05), #BM rate parameter
		extrinsicValues=c(0), #fixed at 0
		timeStep=0.001,
		saveHistory=F)
tf<-list(taxframe)
taxlist<-rep(tf,10000)
convertTaxonFrameToGeigerLapply<-function(taxonframe, phy) {
  Rnames<-function(ids){
    names<-phy$tip.label[ids]
  }
  ntax<-dim(taxonframe)[1]
  newmatrix<-matrix(data=taxonframe[, 4:(dim(taxonframe)[2])], nrow= ntax)
  newrownames<-lapply(taxonframe$taxonid,Rnames)
  geigerframe<-data.frame(newmatrix, row.names= newrownames, stringsAsFactors=F)
  geigerframe
}
convertTaxonFrameToGeigerData<-function(taxonframe, phy) {
	ntax<-dim(taxonframe)[1]
	newmatrix<-matrix(data=taxonframe[, 4:(dim(taxonframe)[2])], nrow= ntax)
	newrownames<-c(rep(0, ntax))
	for (i in 1:ntax) {
		newrownames[i]<-phy$tip.label[(taxonframe$taxonid[i])]
	}
	geigerframe<-data.frame(newmatrix, row.names= newrownames, stringsAsFactors=F)
	geigerframe
}

a<-system.time(lapply(taxlist,convertTaxonFrameToGeigerLapply,phy))
b<-system.time(lapply(taxlist,convertTaxonFrameToGeigerData,phy))
