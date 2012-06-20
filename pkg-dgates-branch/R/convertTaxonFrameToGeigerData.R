convertTaxonFrameToGeigerData<-function(taxonframe, phy) {
	ntax<-dim(taxonframe)[1]
	newmatrix<-matrix(data=taxonframe[, 4:(dim(taxonframe)[2])], nrow= ntax)
	newrownames<-c(rep(0, ntax))
        ntaxVec<-c(1:ntax)
        newRowfn<-function(i) return(phy$tip.label[(taxonframe$taxonid[i])])
        newrownames<-sapply(ntaxVec,newRowfn)
	geigerframe<-data.frame(newmatrix, row.names= newrownames, stringsAsFactors=F)
	geigerframe
}
