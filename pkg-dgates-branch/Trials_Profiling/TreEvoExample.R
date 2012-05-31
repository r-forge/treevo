library(TreEvo)
library(foreach)
library(doMC)
system(command=paste("rm .Rdata"))
set.seed(793450)


# phylogeny 
phy<-rcoal(20)
splits<-getSimulationSplits(phy)
phy$edge.length<-phy$edge.length/max(branching.times(phy))

#Character evolution along phylogeny using TreEvo models
char<-convertTaxonFrameToGeigerData(
	doSimulation(
		splits=splits, 
		intrinsicFn= brownianIntrinsic, #simple Brownian model
		extrinsicFn= nullExtrinsic, #null
		startingStates=c(10), #root state
		intrinsicValues=c(0.05), #BM rate parameter
		extrinsicValues=c(0), #fixed at 0
		timeStep=0.001,
		saveHistory=F), phy
)

#Calculate summary statistics on characters
summaryStatsLong(phy, char)

#ABC run on simulated phylogeny and characters
#Rprof('doRun_Omeara_2.out',interval=0.9)
a<-doRun(
	phy = phy,
	traits = char,
	intrinsicFn=brownianIntrinsic, #name which intrinsic model
	extrinsicFn=nullExtrinsic, #name which extrinsic model
	startingPriorsFns="normal", #distribution for starting state prior
	startingPriorsValues=matrix(c(mean(char[,1]), sd(char[,1]))), #prior values
	intrinsicPriorsFns=c("exponential"), #distribution for intrinsic prior
	intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE), #intrinsic prior values
	extrinsicPriorsFns=c("fixed"), #distribution for extrinsic prior
	extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE), #fixed at 0
	TreeYears=1000,
	standardDevFactor=0.2, 	
	plot=F,
	StartSims=200, #number of initial simulations
	epsilonProportion=0.7, #What proportion of initial particles to keep
	epsilonMultiplier=0.7, #What proportion of particles to keep for subsequent generations
	nStepsPRC=5, #number of generations
	maxTries=1,
	numParticles=100, #how many accepted particles to keep for each generation
	debug=F,
	whenToKill=20,
	jobName=4,
	stopRule=F,
	multicore=F, #initial sims (ie StartSims) are set up to be multithreaded if T
	coreLimit=1 
)
Rprof(NULL)
