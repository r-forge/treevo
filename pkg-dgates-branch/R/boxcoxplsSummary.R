boxcoxplsSummary<-function(todo, summaryValues, prunedPlsResult, boxcoxLambda, boxcoxAddition) {
  sumValvec<-c(1:length(summaryValues))
  sumValsfn<-function(summaryValueIndex) return((summaryValues[summaryValueIndex] + boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex])
		#if (todo[summaryValueIndex]==1) {
			#print(paste("summaryValueIndex =", summaryValueIndex))
			#print(paste("summaryValues[summaryValueIndex] (before)=",summaryValues[summaryValueIndex]))
		#}
  summaryValues<-sapply(sumValvec,sumValsfn)
		#if (todo[summaryValueIndex]==1) {
			#print(paste(" boxcoxAddition[summaryValueIndex]=",boxcoxAddition[summaryValueIndex]))
			#print(paste(" boxcoxLambda[summaryValueIndex]=",boxcoxLambda[summaryValueIndex]))
			#print(paste(" (sum)^thingie=",(summaryValues[summaryValueIndex] + boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]))
			#print(paste("summaryValues[summaryValueIndex] (after)=",summaryValues[summaryValueIndex]))
		#}
  prunedSummaryValues<-summaryValues[which(todo>0)]
  predictResult<-(predict(prunedPlsResult, matrix(prunedSummaryValues, nrow=1))$predict[, , 1])
  predictResult
}               

