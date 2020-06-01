### Run model in apical configuration first ###

### Test MK18 for perturbations of the PI pathway - P1 phenomena ###



##################################### Perturb function ####################################
Perturb = function(ind)
	{

	perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA'))	# initiation of perturbations matrix, if unaltered it will do no perturbations

	# Mega perturbation
	### PI45P2 sensitive to pi_4K, pi_4K_plus_pip_5KI.
	# Finaly remove the # form the front of all the next code lines
	# On the first column put the time of the perturbation
	ptime=c(0,
		1000,1500,
		2000,2005,2500,2505,
		3000,3005,3500,3505,
		4000,4500,
		5000,5500,6000,6500,
		7000,7500,8000,8500,
		9000,9500,
		max(t));

	# On the second column put the variables to be altered. (ex: 'X')
	pvar=c('novar',
		'gamma_0','gamma_0',
		'pi_4K','pi_4K_pip_5KI','pi_4K','pi_4K_pip_5KI',
		'pip_5KI','pi_4K_pip_5KI','pip_5KI','pi_4K_pip_5KI',
		'MTMR1_6_14','MTMR1_6_14',
		'pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve',
		'PLC','PLC','PLC','PLC',
		'DGK','DGK',
		'novar')
	
	# cbind(1:length(parameters),parameters)
	# On the third the new value for the variable.
	pval=as.numeric(c(NA,
				ind[names(ind)=='gamma_0']*.5,ind[names(ind)=='gamma_0'],
				ind[names(ind)=='pi_4K']*0.1,ind[names(ind)=='pi_4K_pip_5KI']*0.1,ind[names(ind)=='pi_4K'],ind[names(ind)=='pi_4K_pip_5KI'],
				ind[names(ind)=='pip_5KI']*0.5,ind[names(ind)=='pi_4K_pip_5KI']*.5,ind[names(ind)=='pip_5KI'],ind[names(ind)=='pi_4K_pip_5KI'],
				ind[names(ind)=='MTMR1_6_14']*.65,ind[names(ind)=='MTMR1_6_14'],
				ind[names(ind)=='pi_Kfyve']*.1,ind[names(ind)=='pi_Kfyve'],ind[names(ind)=='pi_Kfyve']*.001,ind[names(ind)=='pi_Kfyve'],
				ind[names(ind)=='PLC']*0,ind[names(ind)=='PLC'],ind[names(ind)=='PLC']*5,ind[names(ind)=='PLC'],
				ind[names(ind)=='DGK']*0.33,ind[names(ind)=='DGK'],
				NA))

	perturbations=data.frame(time=ptime,var=pvar,val=pval)
	perturbations
	}

# Perturb(parameters)






############### score of the inicial solution to normalize #################
# creating the base_score to normalize the socres of the diferent solutions. With this normalization each check will be 1. The total will be the number of checks.
base_score=c(1.184110e+02, 1.712800e+00, 9.465740e+00, 1.339280e+02, 2.929200e+03,
		8.705000e+01, 3.365690e+00, 2.840167e+01, 7.791853e-01, 4.999998e-01,
		5.618875e-02, 1.592420e-01, 8.008198e-01, 4.009873e-02, 1.326015e-02,
		6.679521e-04, 5.501719e-02, 7.853117e-03, 1.071201e-01, 3.888513e+00,
		2.695487e-01, 2.163733e-05, 1.000000e-07, 1.487321e-04, 9.488071e-02,
		1.926267e-05, 2.053000e+03
		)





###########################   Score   ##################################

# ind = parameters
# Ss_ind = state
# bscore = base_score

Score = function (ind,Ss_ind,bscore=base_score,Occam=FALSE,only_score=FALSE)
	{
	# Insert the number of time units that the simulation will run
	tmax=10000
	# Insert the step of the simulation
	tstep=1
	t=seq(0,tmax+1,tstep) 	# time
	
	# put names on steady states, if they do not have them.
	names(Ss_ind)=names(state) 

	perturbations2 = Perturb(ind)
	out_score = Cruncher(Ss_ind,t,equations,ind,perturbations2)
	# edit(out_score)

	score = c(
			abs(Ss_ind[3]-Ss_ind[6])/bscore[1],				# similar levels of PI4P and PI45P2
			abs(Ss_ind[2]-Ss_ind[4])/bscore[2],				# similar levels of PI3P and PI5P
			abs(Ss_ind[4]-5*Ss_ind[5])/bscore[3],			# PI5P levels are 5 fold of PI35P2 levels
			abs(Ss_ind[3]-100*Ss_ind[4])/bscore[4],	 		# PI4P or PI45P2 is 100 times more than PI5P or PI3P
			abs(Ss_ind[1]-300000)/bscore[5],				# Ss value of PI close to 300000
			abs(Ss_ind[6]-10000)/bscore[6],				# Ss value of PI45P2 close to 10000
			abs(Ss_ind[2]-100)/bscore[7],					# Ss value of PI3P close to 100
			abs(state[5]-state[7])/bscore[8],				# similar levels of PI35P2 and PI34P2

			abs(out_score[1490,2]/Ss_ind[1]-out_score[1490,7]/Ss_ind[6])/bscore[9],		# PI45P2 will decrease in the same percentage as PI

			abs(0.5-out_score[2490,4]/Ss_ind[3])/bscore[10],	# PI4P drops to .5 after PI4K knockout
			abs(0.5-out_score[2490,7]/Ss_ind[6])/bscore[11],	# PI45P2 drops to .5 after PI4K knockout
			abs(0.5-out_score[3490,7]/Ss_ind[6])/bscore[12],	# PI45P2 drops to .5 after PIP_5KI knockdown

			abs(0.2-out_score[4490,5]/Ss_ind[4])/bscore[13],	# PI5P drops to .2 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)
			abs(1.5-out_score[4490,6]/Ss_ind[5])/bscore[14],	# PI35P2 raizes to 1.5 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)
			abs(0.5-out_score[5490,5]/Ss_ind[4])/bscore[15],	# PI5P should drop to 50% if pi_Kfive is reduced to 10%
			abs(0.5-out_score[5490,6]/Ss_ind[5])/bscore[16],	# PI35P2 should drop to 50% if pi_Kfive is reduced to 10%
			abs(0.15-out_score[6490,5]/Ss_ind[4])/bscore[17],	# PI5P should drop to .15 if pi_Kfyve is knockout
			abs(0.001-out_score[6490,6]/Ss_ind[5])/bscore[18],	# PI35P2 should drop to .001 if pi_Kfyve is knockout
			abs(0.8-out_score[6490,7]/Ss_ind[6])/bscore[19],	# PI45P2 should drop to .8 if pi_Kfyve is knockout 
			abs(5-out_score[6490,3]/Ss_ind[2])/bscore[20],		# PI3P should increase 5 fold if pi_Kfyve is knockout

			abs(1.5-out_score[7490,7]/Ss_ind[6])/bscore[21],		# PI45P2 should increase to 1.5 when PLC is knockout
			abs(1-out_score[7990,7]/Ss_ind[6])/bscore[22],		# PI45P2 should return to normal when PLC is back to basal levels
			ifelse(out_score[8490,7]/state[6]>0.8 & out_score[8490,7]/state[6]<1,0,1)/bscore[23],		# PI45P2 should recover above .8 when PLC is activated, I think, 5 fold
			abs(1-out_score[8990,7]/Ss_ind[6])/bscore[24],		# PI45P2 should return to normal when PLC is back to basal levels
			abs(.428-out_score[9490,14]/Ss_ind[13])/bscore[25],	# PA should drop to .428 when DGK is knockdown
			abs(1-out_score[9990,14]/Ss_ind[13])/bscore[26],		# PA should return to normal when DGK is back to basal levels

			sum(ind[1:18])/base_score[27]					# Occam's razor
			)
	names(score)=NULL
	score=as.numeric(score)

	score_names=c(
	'Similar levels of PI4P and PI45P2',
	'Similar levels of PI3P and PI5P',
	'PI5P levels are 5 fold of PI35P2 levels',
	'PI4P or PI45P2 is 100 times more than PI5P or PI3P',
	'Ss value of PI close to 300000',
	'Ss value of PI45P2 close to 10000',
	'Ss value of PI3P close to 100',
	'Similar levels of PI35P2 and PI34P2',

	'PI45P2 will decrease in the same percentage as PI',

	'PI4P drops to .5 after PI4K knockout',
	'PI45P2 drops to .5 after PI4K knockout',
	'PI45P2 drops to .5 after PIP_5KI knockdown',

	'PI5P drops to .2 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)',
	'PI35P2 raizes to 1.5 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)',
	'PI5P should drop to 50% if pi_Kfive is reduced to 10%',
	'PI35P2 should drop to 50% if pi_Kfive is reduced to 10%',
	'PI5P should drop to .15 if pi_Kfyve is knockout',
	'PI35P2 should drop to .001 if pi_Kfyve is knockout',
	'PI45P2 should drop to .8 if pi_Kfyve is knockout',	
	'PI3P should increase 5 fold if pi_Kfyve is knockout',
	
	'PI45P2 should increase to 1.5 when PLC is knockout',
	'PI45P2 should return to normal when PLC is back to basal levels',
	'PI45P2 should recover above .8 when PLC is activated, I think, 5 fold',
	'PI45P2 should return to normal when PLC is back to basal levels',
	'PA should drop to .428 when DGK is knockdown',
	'PA should return to normal when DGK is back to basal levels',

	"Occam's razor",
	"Score"
	)

	if (Occam==TRUE)												# if we want to calculate score with Occam's razor (less proteins is best)
		{
		score=c(score,sum(score))									# put the total score in the end of the score vector
		if (only_score==TRUE) {tail(score,1)} else {data.frame(score_names,score)}	# if I just want the total score put only_score=TRUE in the function parameters, else we will see the complete score vector
		} else {												# if we want to calculate score without Occam's razor 
		score[length(score)]=sum(score[1:(length(score)-1)])					# put the total score in the place of Occom's 
		score_names = score_names[-(length(score_names)-1)] 
		if (only_score==TRUE) {tail(score,1)} else {data.frame(score_names,score)}	# if I just want the total score put only_score=TRUE in the function parameters, else we will see the complete score vector
		}

	}

# Score(parameters,state,bs=base_score,Occam=FALSE,only_score=FALSE)






#################################   Perturbations   ###################################
perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA'))	# initiation of perturbations matrix, if unaltered it will do no perturbations
# Define perturbations. 
# The first row of the matrix is allways 0,novar,NA for the cycle use the inicial values at the start
# The last row of the matrix is always max(t),novar,NA for the cycle to work.
# The values in the middle rows are the perturbations to the system of diferential equations.






##############
# Figure 2 a #
##############
### PI45P2 sensitive to PI, pi_4K, pi_4K_plus_pip_5KI.
# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,2000,3000,3005,4000,4005,5000,5005,6000,6005,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','gamma_0','gamma_0',
	'pi_4K','pi_4K_pip_5KI','pi_4K','pi_4K_pip_5KI',
	'pip_5KI','pi_4K_pip_5KI','pip_5KI','pi_4K_pip_5KI',
	'novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='gamma_0']*.5,parameters[names(parameters)=='gamma_0'],
	parameters[names(parameters)=='pi_4K']*0.05,parameters[names(parameters)=='pi_4K_pip_5KI']*0.05,parameters[names(parameters)=='pi_4K'],parameters[names(parameters)=='pi_4K_pip_5KI'],
	parameters[names(parameters)=='pip_5KI']*.5,parameters[names(parameters)=='pi_4K_pip_5KI']*.5,parameters[names(parameters)=='pip_5KI'],parameters[names(parameters)=='pi_4K_pip_5KI'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out_test)

str(out)
varnames=c('Time (s)','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3',
		'pi_3KI_a','PTEN_a','DAG','IP3','PA','Ca','PIP5KI_PI4P','pi_4K_pip_5KI_PI4P','PKC','PLC','GCPR') 
varcolors=c('gray','tan','orange3','orange','orange3','gold3','gold','gold3','yellow','red','deepskyblue3','salmon3','yellow2','salmon','cyan','tan3','tan4','magenta','black','red')
#  cbind(1:length(varnames),varnames,varcolors)

# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
head(out_norm)
# edit(out_norm)

windows()
par(mfrow=c(3,3))
for (i in 2:9) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# edit(out_norm)

### bar plots for P1 ###

# PI45P2 and PI4P plots
# for perturbation 3
# Figure 3

# edit(out_norm)
sample1_PI45P2=out_norm[c(500,1500,3500,5500),7]
sample1_PI4P=out_norm[c(500,1500,3500,5500),4]
joint=rbind(sample1_PI45P2,sample1_PI4P)

# tiff("P2F2_Fig_2a.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mar=c(2,6,2,1))
bp=barplot(joint,beside = TRUE,ylim=c(0,max(joint)+.2),
	ylab="Percentage of steady state levels",
	cex.axis=2.5,cex.lab=2.5
	)
bp
text(bp,joint, labels = format(round(joint,2), 4),pos = 3, cex = 2.5)
legend('top', c("PI(4,5)P2","PI(4)P"), cex=2.5, bty="n",fill=c("black","grey"))
points(x=c(4.5,5.5,7.5,8.5,10.5,11.5),y=c(rep(.5,5),1),col='blue', pch = "_", cex = 7)

# dev.off()






##############
# figure 2 b #
##############
### pi_4K, pi_4K_plus_pip_5KI, and PI5P can mantain PI45P2
# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,1500,2000,2500,3000,4000,5000,6000,7000,8000,9000,9500,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','gamma_4','pi_4K_pip_5KI','gamma5_45','pi_4K','pi_4K','pi_4K','pi_4K_pip_5KI','pi_4K_pip_5KI','gamma_4','gamma_4','gamma5_45','gamma5_45','novar')
# On the third the new value for the variable.
pval=c(NA,
	0,0,0,parameters[names(parameters)=='pi_4K']*0.2,
	parameters[names(parameters)=='pi_4K'],parameters[names(parameters)=='pi_4K']*0.2,
	parameters[names(parameters)=='pi_4K_pip_5KI'],0,
	parameters[names(parameters)=='gamma_4'],0,
	parameters[names(parameters)=='gamma5_45'],0,
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out_test)

str(out)
varnames=c('Time (s)','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3',
		'pi_3KI_a','PTEN_a','DAG','IP3','PA','Ca','PIP5KI_PI4P','pi_4K_pip_5KI_PI4P','PKC','PLC','GCPR') 
varcolors=c('gray','tan','orange3','orange','orange3','gold3','gold','gold3','yellow','red','deepskyblue3','salmon3','yellow2','salmon','cyan','tan3','tan4','magenta','black','red')
#  cbind(1:length(varnames),varnames,varcolors)


# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
head(out_norm)
# edit(out_norm)

windows()
par(mfrow=c(3,3))
for (i in 2:9) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# edit(out_norm)

### bar plots for P1 ###

# PI45P2 and PI4P plots 
# for perturbation 1 
# (Figure 2)

# edit(out_norm)
time_sample1=c(0,2500,3000,5000,7000,9000)
sample1_PI45P2=out_norm[time_sample1+490,7]
sample1_PI4P=out_norm[time_sample1+490,4]
joint=rbind(sample1_PI45P2,sample1_PI4P)

# tiff("P2F2_Fig_2b.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mar=c(2,6,2,1))
bp=barplot(joint,beside = TRUE,ylim=c(0,max(joint)+.06),
	ylab="Percentage of steady state levels",
	cex.axis=2.5,cex.lab=2.5	
	)
bp
text(bp,joint, labels = format(round(joint,2), 4),pos = 3, cex = 2)
legend(x=12,y=1, c("PI(4,5)P2","PI(4)P"), cex=2.5, bty="n",fill=c("black","grey"))

# dev.off()






################
# (Figure 2 d) #
################
# edit(out_norm)
time_sample1=c(0,1000,2500,7000)
sample1_PI45P2=out_norm[time_sample1+490,7]
sample1_PI4P=out_norm[time_sample1+490,4]
joint=rbind(sample1_PI45P2,sample1_PI4P)

windows()
# tiff("P2F2_Fig_2d.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mar=c(2,6,2,1))
bp=barplot(joint,beside = TRUE,ylim=c(0,max(joint)+.06),
	ylab="Percentage of steady state levels",
	cex.axis=2.5,cex.lab=2.5	
	)
bp
text(bp,joint, labels = format(round(joint,2), 4),pos = 3, cex = 2.5)
legend("topright", c("PI(4,5)P2","PI(4)P"), cex=2.5, bty="n",fill=c("black","grey"))

# dev.off()







################
# (Figure 2 c) #
################
### knock down of MTMR (myotubularyns) to 70% should make PI5P 20% and PI35P2 150% 
### knock out of SYNJ_TMEM55 should decrease 15%
### Playing with PIKfyve. 50% of PIkfyve should put PI5P and PI35P2 between 35% and 39% of original steady state levels
# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,501,1000,1500,1505,2000,2005,2500,3000,3500,4000,5000,6000,max(t))
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','MTMR1_6_14','MTMR1_6_14','SYNJ','TMEM55','SYNJ','TMEM55','pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='MTMR1_6_14']*.65,parameters[names(parameters)=='MTMR1_6_14'],
	parameters[names(parameters)=='SYNJ']*0,parameters[names(parameters)=='TMEM55']*0,parameters[names(parameters)=='SYNJ'],parameters[names(parameters)=='TMEM55'],
	parameters[names(parameters)=='pi_Kfyve']*.1,parameters[names(parameters)=='pi_Kfyve'],
	parameters[names(parameters)=='pi_Kfyve']*.001,parameters[names(parameters)=='pi_Kfyve'],
	parameters[names(parameters)=='pi_Kfyve']*.5,parameters[names(parameters)=='pi_Kfyve'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out)

# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
head(out_norm)
# edit(out_norm)

windows()

par(mfrow=c(3,3))
for (i in 2:9) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# edit(out_norm)

#       check_pert check_names check_targets check_darts
#1        MTMR 10%        PI5P         0.200 0.300001244
#2        MTMR 10%      PI35P2         1.500 1.744820658
#3            <NA>        <NA>            NA          NA
#4  SYNJ_TMEM55 0%        PI5P         0.850 0.708696492
#5            <NA>        <NA>            NA          NA
#6    pi_Kfyve 10%        PI5P         0.500 0.463738581
#7    pi_Kfyve 10%      PI35P2         0.500 0.377789832
#8            <NA>        <NA>            NA          NA
#9    pi_Kfyve .1%        PI5P         0.150 0.207333684
#10   pi_Kfyve .1%      PI35P2         0.001 0.005400449
#11   pi_Kfyve .1%      PI45P2         0.800 0.901334972
#12   pi_Kfyve .1%        PI3P         5.000 5.467992807

# check for small PIs
check_pert = c('MTMR 10%','MTMR 10%',NA,'SYNJ_TMEM55 0%',NA,'pi_Kfyve 10%','pi_Kfyve 10%',NA,'pi_Kfyve .1%','pi_Kfyve .1%','pi_Kfyve .1%','pi_Kfyve .1%',NA,'pi_Kfyve 50%','pi_Kfyve 50%')   	
check_names = c('PI5P','PI35P2',NA,'PI5P',NA,'PI5P','PI35P2',NA,'PI5P','PI35P2','PI45P2','PI3P',NA,'PI5P','PI35P2')
check_targets = c(.2,1.5,NA,.85,NA,.5,.5,NA,.15,.001,.80,5,NA,35,39)
check_darts = c(out_norm[900,5],out_norm[900,6],NA,out_norm[1900,5],NA,out_norm[2900,5],out_norm[2900,6],NA,out_norm[3900,5],out_norm[3900,6],out_norm[3900,7],out_norm[3900,3],NA,out_norm[4900,5],out_norm[4900,6])
names(check_darts)=NULL
check = data.frame(check_pert,check_names,check_targets,check_darts) 
check

perturbations

# tiff("P2F2_Fig_2c.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

data=c(out_norm[900,5],out_norm[900,6],NA,out_norm[1900,5],NA,out_norm[2900,5],out_norm[2900,6],NA,out_norm[3900,5],out_norm[3900,6],out_norm[3900,7],out_norm[3900,3])	
par(mar=c(8,6,2,1))
bp=barplot(data,
		col=c('gray92','gray','blue','gray92','blue','gray92','gray','blue','gray92','gray','dimgray','black'),
		ylab="Percentage of steady state levels",ylim=c(0,max(data,na.rm = TRUE)+2),
		cex.lab=2.5,cex.axis=2.5,cex.names=2,las=2
		)
bp
text(bp,data ,labels = format(round(data,2), 4),pos = 3, cex = 2.5)
legend("topleft", c("PI(5)P","PI(3,5)P2","PI(4,5)P2","PI(3)"), cex=2.5, bty="n",fill=c('gray92','gray','dimgray','black'))
points(x=bp,
	y=c(.2,1.5,-1,.85,-1,.5,.5,-1,.15,.001,.80,5),
	col='blue', pch = "_", cex = 5)

# dev.off()






