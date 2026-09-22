# Code for chapter 7, table 7.6: simulation of double sampling estimator
rm(list = ls(all = TRUE))
gc()

source("http://hdl.handle.net/10079/dfn2zfw")

# binary outcomes
# 50% of observations randomly assigned to treatment
# attrition depends flexibly on the underlying "type" of the observation, formed by the 2x2 table of potential outcomes


# 18 simulation settings, corresponding to the 18 rows of table 7.6
settings <- matrix(NA,18,9)
settings[1,] <- c(.5,.5,.5,.5,.05,.05,.05,.05,100)
settings[2,] <- c(.5,.5,.5,.5,.05,.05,.05,.05,400)

settings[3,] <- c(.5,.5,.5,.5,.25,.25,.25,.25,100)
settings[4,] <- c(.5,.5,.5,.5,.25,.25,.25,.25,400)

settings[5,] <- c(.5,.5,.5,.5,.5,.5,.5,.5,100)
settings[6,] <- c(.5,.5,.5,.5,.5,.5,.5,.5,400)

settings[7,] <- c(.5,.5,.7,.3,.05,.05,.1,0,100)
settings[8,] <- c(.5,.5,.7,.3,.05,.05,.1,0,400)

settings[9,] <- c(.5,.5,.7,.3,.25,.25,.3,0,100)
settings[10,] <- c(.5,.5,.7,.3,.25,.25,.3,0,400)

settings[11,] <- c(.5,.5,.7,.3,.5,.5,.55,.45,100)
settings[12,] <- c(.5,.5,.7,.3,.5,.5,.55,.45,400)

settings[13,] <- c(.5,.5,.7,.3,.05,.05,.05,.05,100)
settings[14,] <- c(.5,.5,.7,.3,.05,.05,.05,.05,400)

settings[15,] <- c(.5,.5,.7,.3,.25,.25,.25,.25,100)
settings[16,] <- c(.5,.5,.7,.3,.25,.25,.25,.25,400)

settings[17,] <- c(.5,.5,.7,.3,.5,.5,.5,.5,100)
settings[18,] <- c(.5,.5,.7,.3,.5,.5,.5,.5,400)

final.tab <- matrix(NA,18,8)


k <- 1
for(k in 1:18)	{

	# Types are:
	# "I"   # success only under treatment
	# "II"  # success only under control
	# "III" # never success
	# "IV"  # always success 

	# define probabilities of missingness in FIRST round depending on type for CONTROL obs
	typeI.p.1.Co <- settings[k,1]
	typeII.p.1.Co <- settings[k,2]
	typeIII.p.1.Co <- settings[k,3]
	typeIV.p.1.Co <- settings[k,4]

	# define probabilities of missingness in SECOND round depending on type for CONTROL
	# obs
	typeI.p.2.Co <- settings[k,5]
	typeII.p.2.Co <- settings[k,6]
	typeIII.p.2.Co <- settings[k,7]
	typeIV.p.2.Co <- settings[k,8]

	# define probabilities of missingness in FIRST round depending on type for TREATED obs
	typeI.p.1.Tr <- 0
	typeII.p.1.Tr <- 0
	typeIII.p.1.Tr <- 0
	typeIV.p.1.Tr <- 0

	# define probabilities of missingness in SECOND round depending on type for TREATED
	#  obs
	typeI.p.2.Tr <- 0
	typeII.p.2.Tr <- 0
	typeIII.p.2.Tr <- 0
	typeIV.p.2.Tr <- 0


	sims <- 10000	# number of simulation runs
	n <- 10000      # sample size


	y0.prob <- .45  # probability of Y0 == 1
	y1.prob <- .55  # probability of Y1 == 1


	# number of observations with initial missingness that are chosen for follow-up
	follow.up.Co.no <- settings[k,9]
	follow.up.Tr.no <- settings[k,9]


	res <- matrix(NA,sims,14)
	colnames(res) <- c( "Avg.Y1.naive","Avg.Y0.naive","Avg.Y1.lower","Avg.Y1.upper",
                    	"Avg.Y0.lower","Avg.Y0.upper","ATE","ATE lower", "ATE upper",
                    	"Avg.Y1.lessnaive","Avg.Y0.lessnaive", "ATE.lessnaive",
                    	"Lower ATE bound naive","Upper ATE bound naive")

	set.seed(123) # note that this is a different seed from the one used for the published table
	

	m <- 1
	for(m in 1:sims)    {

    	# create schedule of potential outcomes
    	data1 <- matrix(NA,n,6)
    	colnames(data1) <- c("Tr","Y0","Y1","M1","M2","Y")
    	data1 <- as.data.frame(data1)

    	# generate treatment assignments
    	data1$Tr <- rbinom(n,1,.5)

    	# generate potential outcomes
    	data1$Y0 <- rbinom(n,1,y0.prob)
    	data1$Y1 <- rbinom(n,1,y1.prob)

    	# record type of each observation
    	type <- rep(NA,n)
    	type[data1$Y0 == 0 & data1$Y1 == 1] <- "I"      # success only under treatment
    	type[data1$Y0 == 1 & data1$Y1 == 0] <- "II"     # success only under control
    	type[data1$Y0 == 0 & data1$Y1 == 0] <- "III"    # never success
    	type[data1$Y0 == 1 & data1$Y1 == 1] <- "IV"     # always success
		table(type)

    	# generate missingness indicator for first stage as a function of treatment assignment and type
    	data1$M1[data1$Tr == 0 & type == "I"] <-
    	rbinom(n,1,typeI.p.1.Co)[data1$Tr == 0 & type == "I"]
    	data1$M1[data1$Tr == 0 & type == "II"] <-
    	rbinom(n,1,typeII.p.1.Co)[data1$Tr == 0 & type == "II"]
    	data1$M1[data1$Tr == 0 & type == "III"] <-
    	rbinom(n,1,typeIII.p.1.Co)[data1$Tr == 0 & type == "III"]
    	data1$M1[data1$Tr == 0 & type == "IV"] <-
    	rbinom(n,1,typeIV.p.1.Co)[data1$Tr == 0 & type == "IV"]

    	data1$M1[data1$Tr == 1 & type == "I"] <-
    	rbinom(n,1,typeI.p.1.Tr)[data1$Tr == 1 & type == "I"]
    	data1$M1[data1$Tr == 1 & type == "II"] <-
    	rbinom(n,1,typeII.p.1.Tr)[data1$Tr == 1 & type == "II"]
    	data1$M1[data1$Tr == 1 & type == "III"] <-
    	rbinom(n,1,typeIII.p.1.Tr)[data1$Tr == 1 & type == "III"]
    	data1$M1[data1$Tr == 1 & type == "IV"] <-
    	rbinom(n,1,typeIV.p.1.Tr)[data1$Tr == 1 & type == "IV"]


    	# generate missingness indicator for second stage as a function of treatment assignment and type
    	data1$M2[data1$Tr == 0 & type == "I" & data1$M1 == 1] <-
    	rbinom(n,1,typeI.p.2.Co)[data1$Tr == 0 & type == "I" & data1$M1 == 1]
    	data1$M2[data1$Tr == 0 & type == "II" & data1$M1 == 1] <-
    	rbinom(n,1,typeII.p.2.Co)[data1$Tr == 0 & type == "II" & data1$M1 == 1]
    	data1$M2[data1$Tr == 0 & type == "III" & data1$M1 == 1] <-
    	rbinom(n,1,typeIII.p.2.Co)[data1$Tr == 0 & type == "III" & data1$M1 == 1]
    	data1$M2[data1$Tr == 0 & type == "IV" & data1$M1 == 1] <-
    	rbinom(n,1,typeIV.p.2.Co)[data1$Tr == 0 & type == "IV" & data1$M1 == 1]

    	data1$M2[data1$Tr == 1 & type == "I" & data1$M1 == 1] <-
    	rbinom(n,1,typeI.p.2.Tr)[data1$Tr == 1 & type == "I" & data1$M1 == 1]
    	data1$M2[data1$Tr == 1 & type == "II" & data1$M1 == 1] <-
    	rbinom(n,1,typeII.p.2.Tr)[data1$Tr == 1 & type == "II" & data1$M1 == 1]
    	data1$M2[data1$Tr == 1 & type == "III" & data1$M1 == 1] <-
    	rbinom(n,1,typeIII.p.2.Tr)[data1$Tr == 1 & type == "III" & data1$M1 == 1]
    	data1$M2[data1$Tr == 1 & type == "IV" & data1$M1 == 1] <-
    	rbinom(n,1,typeIV.p.2.Tr)[data1$Tr == 1 & type == "IV" & data1$M1 == 1]

    	data1$M2[data1$M1 == 0] <- 0    # obs not missing in first round are not missing in follow-up


    	# determine observed Y
    	data1$Y[data1$Tr == 1 & data1$M1 == 0] <- data1$Y1[data1$Tr == 1 & data1$M1 == 0]
    	data1$Y[data1$Tr == 0 & data1$M1 == 0] <- data1$Y0[data1$Tr == 0 & data1$M1 == 0]

    	# naively calculate average observed outcomes, discarding missing observations
    	res[m,1] <- mean(data1$Y[data1$Tr == 1], na.rm = TRUE)
    	res[m,2] <- mean(data1$Y[data1$Tr == 0], na.rm = TRUE)


    	# randomly choose fraction of obs with initial missingness for detailed follow-up
    	control.M1  <- which(data1$Tr == 0 & data1$M1 == 1)

    	follow.up.Co <- follow.up.Co.no / length(control.M1)
    	if(follow.up.Co > 1) follow.up.Co <- 1

    	control.M2  <- sample(control.M1, floor(length(control.M1)*follow.up.Co), replace = FALSE)
    	control.M2a <- is.element(1:n, control.M2)

    	data1$Y[control.M2a & data1$M2 == 0] <- data1$Y0[control.M2a & data1$M2 == 0]


    	treated.M1  <- which(data1$Tr == 1 & data1$M1 == 1)

    	follow.up.Tr <- follow.up.Tr.no / length(treated.M1)
    	if(follow.up.Tr > 1) follow.up.Tr <- 1

    	treated.M2  <- sample(treated.M1, floor(length(treated.M1)*follow.up.Tr), replace = FALSE)
    	treated.M2a <- is.element(1:n, treated.M2)

    	data1$Y[treated.M2a & data1$M2 == 0] <- data1$Y1[treated.M2a & data1$M2 == 0]

    	# consistency check:
    	stopifnot(sum(is.na(data1[,1:5])) == 0)


    	# calculate Manski bounds
    	temp <- manski2(data1 = data1,
    	control.M1 = control.M1, control.M2 = control.M2, control.M2a = control.M2a,
    	treated.M1 = treated.M1, treated.M2 = treated.M2, treated.M2a = treated.M2a)

    	res[m,3:7] <- temp[1:5]
    	res[m,10:11] <- temp[6:7]

    	# lower and upper bounds on ATE
    	res[m,8] <- res[m,3] - res[m,6]
    	res[m,9] <- res[m,4] - res[m,5]
    	res[m,12] <- res[m,10] - res[m,11]

    	# lower and upper bounds on naive ATE
    	res[m,13] <- temp[10] - temp[9]
    	res[m,14] <- temp[11] - temp[8]

    	cat(m,"\n")

    	}

	# consistency check
	stopifnot(sum(is.na(res)) == 0)

	# naive estimator
	final.tab[k,1 ] <- round(mean(res[,1] - res[,2]),4)
	final.tab[k,2 ] <- round(sd(res[,1] - res[,2]),4)

	# BOUNDS naive estimator
	final.tab[k,3 ] <- round(mean(res[,13]),4)
	final.tab[k,4 ] <- round(mean(res[,14]),4)

	# double sampling estimator, ignoring missingness
	final.tab[k,5 ] <- round(mean(res[,12]),4)
	final.tab[k,6 ] <- round(sd(res[,12]),4)

	# double sampling estimator, bounds
	final.tab[k,7 ] <- round(mean(res[,8]),4)
	final.tab[k,8 ] <- round(mean(res[,9]),4)

	cat("\n","\n","\n",k,"\n","\n","\n")

	}


final.tab
