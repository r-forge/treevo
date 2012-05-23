library(TreEvo)
library(foreach)
library(doMC)

#Create a phylogeny 
phy<-rcoal(10)
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
Rprof('doRun_labMac.out',interval=0.9)		
a<-doRun(
	phy = phy,
	traits = char,
	intrinsicFn=brownianIntrinsic, #name which intrinsic model
	extrinsicFn=nullExtrinsic, #name which extrinsic model
	startingPriorsFns="uniform", #distribution for starting state prior
	startingPriorsValues=matrix(c(min(char[,1]), max(char[,1]))), #prior values
	intrinsicPriorsFns=c("exponential"), #distribution for intrinsic prior
	intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE), #intrinsic prior values
	extrinsicPriorsFns=c("fixed"), #distribution for extrinsic prior
	extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE), #fixed at 0
	TreeYears=10000,
	standardDevFactor=0.2, 	
	plot=F,
	StartSims=10, #number of initial simulations
	epsilonProportion=0.5, #What proportion of initial particles to keep
	epsilonMultiplier=0.5, #What proportion of particles to keep for subsequent generations
	nStepsPRC=5, #number of generations
	maxTries=1,
	numParticles=25, #how many accepted particles to keep for each generation
	debug=F,
	whenToKill=20,
	jobName=6,
	stopRule=F,
	multicore=F, #initial sims (ie StartSims) are set up to be multithreaded if T
	coreLimit=1 
)
Rprof(NULL)


