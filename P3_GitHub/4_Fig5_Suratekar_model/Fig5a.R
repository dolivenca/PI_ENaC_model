#################################################
#								#
#	Surdekar and ENaC ASL model 5 Po 3		#
#								#
#################################################

rm(list=ls())	#clean memory
graphics.off() 	#close plot windows

#install.packages("beepr")
library(beepr)
#beep()

library(deSolve)	#library for solve diferential equations
#citation("deSolve")

# Insert the number of time units that the simulation will run
tmax=3000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

ENaC_op = function(PIP2)
	{.02+.8/(1+181429*exp(-.00112*PIP2))}
ENaC_op(10000)

# Inicial states for variables
# insert, in the parentesis, the inicial value of the dependent variables (ex: X=2) 
state=c(
	PMPI = 0.6633562,
	PI4P = 0.03816765,
	PIP2 = 0.03842087,
	DAG = 0.006625219,
	PMPA = 0.01672982,
	ERPA = 0.1172366,
	CDPDAG = 0.0007857711,
	ERPI = 0.1169003,
	 
 	ENaC = 35,
	ASL = 7
	) 

# insert, in the parentecis, the parameters and independent variables and their values (ex: const=5, rate=2). 
# If no parameters and independent variables just put a zero.
parameters=c(
		dagk_Vmax = 7.545,
		laza_Vmax = 21.196,
		pip5k_Vmax = 1.055,
		
		SPLUNC1 = 7714,

		gamma__ENaC = log(2)/40,
		gamma_ENaC_SPLUNC1 = 9*log(2)/24684800,
		gamma_ENaC_ = log(2)/3200,
		
		gamma_ASL_ENaC = 0.00026
		)   
parameters=c(
		parameters,
		gamma__ASL = (455/2)*as.numeric(parameters[names(parameters)=='gamma_ASL_ENaC']),
		gamma_CFTR_ASL = (735/8)*as.numeric(parameters[names(parameters)=='gamma_ASL_ENaC']),
		gamma_ASL_ = (295/8)*as.numeric(parameters[names(parameters)=='gamma_ASL_ENaC'])	
		)
 

equations=function(t,state,parameters) 	# function containing the diferential equations
	{with(as.list(c(state,parameters)),
		{
		#rate of change (velocities of the variables concentration changes), the actual diferencial equations
		# insert the diferentical equations (ex: dX=const*X^rate or dX=5*X^2)

		V_pip5k = pip5k_Vmax * PI4P/(PI4P + 0.109)
		V_plc = 1 * PIP2/(PIP2 + 0.102)
		V_laza =  laza_Vmax * PMPA/(PMPA + .57)
		V_dagk = dagk_Vmax * DAG/(DAG + 0.076)
		V_sink = 4.27 * DAG /(DAG + 0.097)
		V_patp = .03 * PMPA/(PMPA + .803)
		V_source = 0.273
		V_cds = .521 * ERPA/(ERPA + .106)
		V_pis = 15.943 * CDPDAG/(CDPDAG + 0.045)
		V_pitp = 3.314* ERPI/(ERPI + 1.299)
		V_pi4k = 1.815* PMPI/(PMPI + 3.737)

		dPMPI = V_pitp - V_pi4k
		dPI4P = V_pi4k - V_pip5k
		dPIP2 = V_pip5k - V_plc
		dDAG = V_plc + V_laza - V_dagk - V_sink
		dPMPA = V_dagk - V_laza - V_patp
		dERPA = V_patp + V_source - V_cds
		dCDPDAG = V_cds  - V_pis
		dERPI = V_pis - V_pitp

		################## ENaC ASL model ####################

		V__ENaC = gamma__ENaC
		V_ENaC_SPLUNC1 = gamma_ENaC_SPLUNC1 * ENaC * (SPLUNC1/ASL) 
		V_ENaC_ = gamma_ENaC_ * ENaC

		V__ASL = gamma__ASL 
		V_CFTR_ASL = gamma_CFTR_ASL
		V_ASL_ENaC = gamma_ASL_ENaC * ASL * (ENaC * ENaC_op(PIP2*10000/0.03842087))
		V_ASL_ = gamma_ASL_ * ASL

		dENaC = V__ENaC - V_ENaC_SPLUNC1 - V_ENaC_ 
		dASL = V__ASL + V_CFTR_ASL - V_ASL_ENaC - V_ASL_


		# return the rate of change 
		# insert, in the parentecis, the variables that you want to store the values (ex:dX)
		list(dy=	c(dPMPI,dPI4P,dPIP2,dDAG,dPMPA,dERPA,dCDPDAG,dERPI,dENaC,dASL),
				count=NULL
			)
		})
	}

perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA'))	# initiation of perturbations matrix, if unaltered it will do no perturbations





##################################### Perturb ####################################

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
#ptime=c(0,1000,2000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
#pvar=c('novar','pip5k_Vmax','pip5k_Vmax','novar')
# On the third the new value for the variable.
#pval=c(NA,
#	parameters[names(parameters)=='pip5k_Vmax']*.1,
#	parameters[names(parameters)=='pip5k_Vmax'],
#	NA)
#perturbations=data.frame(time=ptime,var=pvar,val=pval)
#perturbations



# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
#ptime=c(0,1000,2000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
#pvar=c('novar','dagk_Vmax','dagk_Vmax','novar')
# On the third the new value for the variable.
#pval=c(NA,
#
#	parameters[names(parameters)=='dagk_Vmax']*.1,
#	parameters[names(parameters)=='dagk_Vmax'],
#	NA)
#perturbations=data.frame(time=ptime,var=pvar,val=pval)
#perturbations




################################## Cruncher function - ode integrator that deals with perturbations ###############################
Cruncher = function (state,t,equations,parameters,perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA')))
	{

	# running the ode solver
	out_test=NA #cleanning out_MC vector

	# initializating the matrix out1 and out, necessary for the cycle
	out1=rbind(rep(0,1,length(state)+1),c(t=0,state))
	out_test=matrix()

	if (perturbations[1,1]=='NA') # if no perturbations an ordinary ODE is run
		{out_test <- ode(y = state, times = t, func = equations, parms=parameters)}

	if (perturbations[1,1]!='NA') # if there are perturbations, then go to the for cycle
		{
		if (tmax<max(perturbations[1:(nrow(perturbations)-1),1])) {cat("Time error: Last perturbation later than tmax.Please alter tmax or perturbation time.", "\n")}
		for (i in 1:(nrow(perturbations)-1)) # the number of perturbations is nrow(perturbations)-1 because the first line on the perturbations matrix is just to permit the inicial values to be used
			{
			t=seq(max(0,tail(out1[,1],1)+1),perturbations[i+1,1]-1,tstep)  	# set the time intervals #?#
	
			state=out1[nrow(out1),(2:(length(state)+1))]			# retrive the last values of the variables
			for (j in 1:length(state)) 						# check for variables perturbations
				{
				if (names(state)[j]==perturbations[i,2])			# if a perturbation is in a variable this will alter the variable value in the state vector
					{
					state[j]=perturbations[i,3]
					}
				}
			for (k in 1:length(parameters)) 						# check for variables perturbations
				{
				if (names(parameters)[k]==perturbations[i,2])			# if a perturbation is in a variable this will alter the variable value in the state vector
					{
					parameters[k]=perturbations[i,3]
					}
				}

			out1 <- ode(y = state, times = t, func = equations,parms=parameters)	# good and old ode
			if (is.na(out_test[1,1])) {out_test=out1} else {out_test=rbind(out_test,out1)} 			# store the values of the last interval	
			}
		names(out_test)[1]='time'
		}
	out_test
	}

# out_test = Cruncher(state,t,equations,parameters)						# no perturbations
# perturbations1 = Perturb(parameters)
# out_test = Cruncher(state,t,equations,parameters,perturbations)			# with perturbations
# edit(out_test)
  
out = Cruncher(state,t,equations,parameters,perturbations)
# edit(out)

str(out)
varnames=c(
		'Time (s)',
		'PIP2',
		'PMPI',
		'PI4P',
		'DAG',
		'PMPA',
		'ERPA',
		'CDPDAG',
		'ERPI',

		'ENaC',
		'ASL'
		) 

varcolors=c('gray','tan','orange3','orange','orange3','gold3','gold','gold3','yellow','green','cyan')
#  cbind(1:length(varnames),varnames,varcolors)

# naming fluxes in out
#flux_names=NULL
#colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

# cbind(1:length(parameters),parameters)
# cbind(21:(20+length(flux_names)),flux_names)

#par(mfrow=c(3,1))
# time plot
plot(
	out[,1],out[,2],
	type='l',col=1,ylim=c(0,max(out[,2:9])),
	main='PIs Time course',xlab='Time (m)',ylab='Molecule count'	
	)
for (i in c(2:9))
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

windows()
par(mfrow=c(2,1))
plot(
	out[,1],out[,10],
	type='l',col='green',ylim=c(0,max(out[,10])),
	main='ENaC Time course',xlab='Time (m)',ylab='Molecule count'	
	)
plot(
	out[,1],out[,11],
	type='l',col='cyan',ylim=c(0,max(out[,11])),
	main='ASL Time course',xlab='Time (m)',ylab='Molecule count'	
	)
par(mfrow=c(1,1))


# list (primitive)
(head(out,3))[,1:11]
(tail(out,3))[,1:11]

# to get the values in pmoles
state*4.58








###########################
###				###
###	Alma�a's data	###
###				###
###########################

### SPLUNC1 CF / DGK perturbation ###

### WT ###
tmax=30000
tstep=1
t=seq(0,tmax+1,tstep) 	# time

state[names(state)=='ENaC']=35
state[names(state)=='ASL']=7
parameters[names(parameters)=='SPLUNC1']=7714
parameters[names(parameters)=='gamma_CFTR_ASL']= 2.388750e-02

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,10000,20000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','dagk_Vmax','dagk_Vmax','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='dagk_Vmax']*0.01,
	parameters[names(parameters)=='dagk_Vmax'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)
# edit(out)
# View(out)

# time plot
plot(
	out[,1],out[,2],
	type='l',col=1,ylim=c(0,max(out[,2:9])),
	main='PIs Time course',xlab='Time (m)',ylab='Molecule count'	
	)
for (i in c(2:9))
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

windows()
par(mfrow=c(2,1))
plot(
	out[,1],out[,10],
	type='l',col='green',ylim=c(0,max(out[,10])),
	main='ENaC Time course',xlab='Time (m)',ylab='Molecule count'	
	)
plot(
	out[,1],out[,11],
	type='l',col='cyan',ylim=c(0,max(out[,11])),
	main='ASL Time course',xlab='Time (m)',ylab='Molecule count'	
	)
par(mfrow=c(1,1))

time_sample=c(9990,19990)
N_WT_sample = out[time_sample,10]
ASL_WT_sample = out[time_sample,11]
op_WT_sample = ENaC_op(out[time_sample,4]*10000/0.03842087)
activity_WT_sample = N_WT_sample*op_WT_sample

### CF ###
tmax=30000
tstep=1
t=seq(0,tmax+1,tstep) 	# time

state[names(state)=='ENaC'] = 80
state[names(state)=='ASL'] = 4
parameters[names(parameters)=='SPLUNC1'] = 0
parameters[names(parameters)=='gamma_CFTR_ASL'] = 0

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,10000,20000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','dagk_Vmax','dagk_Vmax','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='dagk_Vmax']*0.01,
	parameters[names(parameters)=='dagk_Vmax'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)
# edit(out)
# View(out)

# time plot
plot(
	out[,1],out[,2],
	type='l',col=1,ylim=c(0,max(out[,2:9])),
	main='PIs Time course',xlab='Time (m)',ylab='Molecule count'	
	)
for (i in c(2:9))
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

windows()
par(mfrow=c(2,1))
plot(
	out[,1],out[,10],
	type='l',col='green',ylim=c(0,max(out[,10])),
	main='ENaC Time course',xlab='Time (m)',ylab='Molecule count'	
	)
plot(
	out[,1],out[,11],
	type='l',col='cyan',ylim=c(0,max(out[,11])),
	main='ASL Time course',xlab='Time (m)',ylab='Molecule count'	
	)
par(mfrow=c(1,1))  

time_sample=c(9990,19990)
N_CF_sample = out[time_sample,10]
ASL_CF_sample = out[time_sample,11]
op_CF_sample = ENaC_op(out[time_sample,4]*10000/0.03842087)
activity_CF_sample = N_CF_sample*op_CF_sample

N_sample = c(N_WT_sample,N_CF_sample)
ASL_sample = c(ASL_WT_sample,ASL_CF_sample)
op_sample = c(op_WT_sample,op_CF_sample)
activity_sample = c(activity_WT_sample,activity_CF_sample)
name_perturb = c(	'WT / DGK +',
			'WT / DGK -',
			'CF / DGK +',
			'CF / DGK -'
			)
cbind(name_perturb,N_sample,op_sample,activity_sample,ASL_sample)
joint = rbind(N_sample,op_sample,activity_sample)



graphics.off() 	#close plot windows
windows()
# tiff("P2F2_Fig_suratekar_a.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(2,2))

bp=barplot(joint[1,],
	beside = TRUE,ylim=c(0,max(joint[1,])*1.2),
	ylab="Number of actived ENaC channels (N)",
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[1,], labels = format(round(joint[1,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
#windows()
bp=barplot(joint[2,],
	beside = TRUE,ylim=c(0,1),
	ylab="ENaC open probability (Po)",
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[2,], labels = format(round(joint[2,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
#windows()
bp=barplot(joint[3,],
	beside = TRUE,ylim=c(0,max(joint[3,])*1.2),
	ylab="ENaC activity ( N * Po)",
	col=c('steelblue','steelblue1'),
	cex.main = 2, cex.axis = 1.5, cex.lab = 1.5,
	main = "Suratekar model without PIP5KI regulation"
	)
text(bp,joint[3,], labels = format(round(joint[3,],2), 4),pos = 3, cex = 2)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=1.5)
par(mfrow=c(1,1))

# dev.off()







######################################### Suratekar stuff ###############################################

# Insert the number of time units that the simulation will run
tmax=3000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

##################################### Perturb ####################################

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,500,1000,1500,2000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','laza_Vmax','laza_Vmax','dagk_Vmax','dagk_Vmax','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='laza_Vmax']*.1,
	parameters[names(parameters)=='laza_Vmax'],
	parameters[names(parameters)=='dagk_Vmax']*.1,
	parameters[names(parameters)=='dagk_Vmax'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

  
out = Cruncher(state,t,equations,parameters,perturbations)
# View(out)

str(out)
varnames=c(
		'Time (s)',
		'PMPI',
		'PI4P',
		'PIP2',
		'DAG',
		'PMPA',
		'ERPA',
		'CDPDAG',
		'ERPI',

		'ENaC',
		'ASL'
		) 

varcolors=c('gray','tan','orange3','orange','orange3','gold3','gold','gold3','yellow','green','cyan')
#  cbind(1:length(varnames),varnames,varcolors)

# naming fluxes in out
#flux_names=NULL
#colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

# cbind(1:length(parameters),parameters)
# cbind(21:(20+length(flux_names)),flux_names)

#par(mfrow=c(3,1))
# time plot
#plot(
#	out[,1],out[,2],
#	type='l',col=1,ylim=c(0,max(out[,2:9])),
#	main='PIs Time course',xlab='Time (m)',ylab='Molecule count'	
#	)
#for (i in c(2:9))
#	{
#	points(out[,1],out[,i],type='l',col=varcolors[i])
#	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
#	}

#par(mfrow=c(2,1))
#plot(
#	out[,1],out[,10],
#	type='l',col='green',ylim=c(0,max(out[,10])),
#	main='ENaC Time course',xlab='Time (m)',ylab='Molecule count'	
#	)
#plot(
#	out[,1],out[,11],
#	type='l',col='cyan',ylim=c(0,max(out[,11])),
#	main='ASL Time course',xlab='Time (m)',ylab='Molecule count'	
#	)
#par(mfrow=c(1,1))


# list (primitive)
(head(out,3))[,1:11]
(tail(out,3))[,1:11]


lipid_ratios = c(
			"PI45P2/PItotal",
			"PI4P/PItotal",
			"DAG/PItotal",
			"PAtotal/PItotal",
			"CDPDAG/PItotal"
			)
lipid_values_model = c(
			out[2900,4]/(out[2900,2]+out[2900,9]),				# PI45P2/PItotal
			out[2900,3]/(out[2900,2]+out[2900,9]),				# PI4P/PItotal
			out[2900,5]/(out[2900,2]+out[2900,9]),				# DAG/PItotal
			(out[2900,6]+out[2900,7])/(out[2900,2]+out[2900,9]),		# PAtotal/PItotal
			out[2900,8]/(out[2900,2]+out[2900,9])				# CDPDAG/PItotal
			)
names(lipid_values_model)=NULL
lipid_values = c(.05,.05,.008,.1677,.001)
data.frame(lipid_ratios,lipid_values,lipid_values_model,Dif = abs(lipid_values-lipid_values_model))


windows()
# tiff("P2F2_Fig_suratekar_c.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(2,1))
RdgA3_vals = c(
	((out[1400,6]+out[1400,7])/(out[1400,2]+out[1400,9])) / ((out[1400,6]+out[1400,7])/(out[1400,2]+out[1400,9])),
	((out[1900,6]+out[1900,7])/(out[1900,2]+out[1900,9])) / ((out[1400,6]+out[1400,7])/(out[1400,2]+out[1400,9])),
	(out[1400,5]/(out[1400,2]+out[1400,9])) / (out[1400,5]/(out[1400,2]+out[1400,9])),
	(out[1900,5]/(out[1900,2]+out[1900,9])) / (out[1400,5]/(out[1400,2]+out[1400,9]))	
	)
names(RdgA3_vals) = c('PAtotal / PItotal','PAtotal / PItotal','DAG / PItotal','DAG / PItotal')
barplot(
	RdgA3_vals,
	col=c('darkblue','red','darkblue','red'),
	main = 'RdgA3',cex.main=3,
	cex.axis=2,cex.names=2,
	ylim=c(0,1.2)
		)
text(x=c(.7,1.95,3.1,4.3),y=RdgA3_vals-.1,labels =round(RdgA3_vals,2),col='white',cex=3)


laza22_vals = c(
	((out[400,6]+out[400,7])/(out[400,2]+out[400,9])) / ((out[400,6]+out[400,7])/(out[400,2]+out[400,9])),
	((out[900,6]+out[900,7])/(out[900,2]+out[900,9])) / ((out[400,6]+out[400,7])/(out[400,2]+out[400,9])),
	(out[400,5]/(out[400,2]+out[400,9])) / (out[400,5]/(out[400,2]+out[400,9])),
	(out[900,5]/(out[900,2]+out[900,9])) / (out[400,5]/(out[400,2]+out[400,9]))
	)
names(laza22_vals) = c('PAtotal / PItotal','PAtotal / PItotal','DAG / PItotal','DAG / PItotal')
barplot(
	laza22_vals,
	col=c('darkblue','red','darkblue','red'),
	main = 'Laza22',cex.main=3,cex.names=2,
	cex.axis=2
		)
text(x=c(.7,1.95,3.1,4.3),y=laza22_vals-.20,labels =round(laza22_vals,2),col='white',cex=3)
legend(3,3,legend=c('WT','Mutant'),fill=c('darkblue','red'),bty='n',cex=2)

# dev.off()

beep(sound=1)
