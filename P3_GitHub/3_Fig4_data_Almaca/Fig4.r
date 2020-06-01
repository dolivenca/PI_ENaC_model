### Run model in apical configuration first ###

graphics.off() 	#close plot windows


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
parameters[names(parameters)=='SPLUNC1'] = 7890
parameters[names(parameters)=='gamma_CFTR_ASL'] = 0.0177052536

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,10000,20000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)
# edit(out)
# View(out)

colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

par(mfrow=c(3,1))
# time plot
plot(out[,1],out[,2],type='l',col=1,ylim=c(0,max(out[,2:9])),main='Time course',xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,25000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,3000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))
	
### ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15:16)) 
	{
	plot(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

time_sample=c(9990,19990)
N_WT_sample = out[time_sample,15]
ASL_WT_sample = out[time_sample,16]
op_WT_sample = ENaC_op(out[time_sample,7])
activity_WT_sample = N_WT_sample*op_WT_sample
almaca_plot = activity_WT_sample

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
pvar=c('novar','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)
# edit(out)
# View(out)

colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

par(mfrow=c(3,1))
# time plot
plot(out[,1],out[,2],type='l',col=1,ylim=c(0,max(out[,2:9])),main='Time course',xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,25000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,3000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))
	
### ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15:16)) 
	{
	plot(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],ENaC_op(out[,7]),type='l',col='blue',main='ENaC open probability')   

time_sample=c(9990,19990)
N_CF_sample = out[time_sample,15]
ASL_CF_sample = out[time_sample,16]
op_CF_sample = ENaC_op(out[time_sample,7])
activity_CF_sample = N_CF_sample*op_CF_sample
almaca_plot = c(almaca_plot,activity_CF_sample)

### P2 ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15)) 
	{
	plot(out[1:6000,1],out[1:6000,i],
	type='l',lwd=5,col=varcolors[i],
	xlab='time (min)', 
	ylab='Number of actived ENaC channels'
	)
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[1:6000,1],ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='blue',
	xlab='time (min)',	
	ylab='ENaC open probability')   

plot(out[1:6000,1],out[1:6000,i]*ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='black',
	xlab='time (min)',	
	ylab='ENaC activity')   


N_sample = c(N_WT_sample,N_CF_sample)
ASL_sample = c(ASL_WT_sample,ASL_CF_sample)
op_sample = c(op_WT_sample,op_CF_sample)
activity_sample = c(activity_WT_sample,activity_CF_sample)
name_perturb = c(	'SPLUNC1 + / DGK +',
			'SPLUNC1 + / DGK -',
			'SPLUNC1 - / DGK +',
			'SPLUNC1 - / DGK -'
			)
cbind(name_perturb,N_sample,op_sample,activity_sample,ASL_sample)
joint = rbind(N_sample,op_sample,activity_sample)



# tiff("P2F2_Fig_almaca_a.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(2,2))
#windows()
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
	col=c('steelblue','steelblue1','steelblue','steelblue1'),
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[3,], labels = format(round(joint[3,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
par(mfrow=c(1,1))

# dev.off()






#############################################################################


### PLC inhibition / DGK perturbation ###

### CF ###
tmax=6000
tstep=1
t=seq(0,tmax+1,tstep) 	# time

state[names(state)=='ENaC'] = 80
state[names(state)=='ASL'] = 4
parameters[names(parameters)=='SPLUNC1'] = 0
parameters[names(parameters)=='gamma_CFTR_ASL'] = 0

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,2000,3000,4000,5000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','DGK','DGK','PLC','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	5*0.1,	
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)

colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

par(mfrow=c(3,1))
# time plot
plot(out[,1],out[,2],type='l',col=1,ylim=c(0,max(out[,2:9])),main='Time course',xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,25000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,3000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))
	
### ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15:16)) 
	{
	plot(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],ENaC_op(out[,7]),type='l',col='blue',main='ENaC open probability')   

time_sample=c(990,1990,3990,4990)
N_CF_sample = out[time_sample,15]
ASL_CF_sample = out[time_sample,16]
op_CF_sample = ENaC_op(out[time_sample,7])
activity_CF_sample = N_CF_sample*op_CF_sample
almaca_plot = c(almaca_plot,activity_CF_sample[3:4])

### P2 ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15)) 
	{
	plot(out[1:6000,1],out[1:6000,i],
	type='l',lwd=5,col=varcolors[i],
	xlab='time (min)', 
	ylab='Number of actived ENaC channels'
	)
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[1:6000,1],ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='blue',
	xlab='time (min)',	
	ylab='ENaC open probability')   

plot(out[1:6000,1],out[1:6000,i]*ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='black',
	xlab='time (min)',	
	ylab='ENaC activity')   


N_sample = c(N_CF_sample)
ASL_sample = c(ASL_CF_sample)
op_sample = c(op_CF_sample)
activity_sample = c(activity_CF_sample)

name_perturb = c(	'PLC + / DGK +',
			'PLC + / DGK -',
			'PLC - / DGK +',
			'PLC - / DGK -'
			)
cbind(name_perturb,N_sample,op_sample,activity_sample,ASL_sample)
joint = rbind(N_sample,op_sample,activity_sample)


# tiff("P2F2_Fig_almaca_b.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(2,2))
#windows()
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

# tiff("figS4_c.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)
name_perturb = c(	'SPLUNC1 + \n DGK +',
			'SPLUNC1 + \n DGK -',
			'SPLUNC1 - \n DGK +',
			'SPLUNC1 - \n DGK -')
par(mar=c(5, 6, 4, 2) + 0.1)
bp=barplot(joint[3,],
	beside = TRUE,ylim=c(0,max(joint[3,])*1.2),
	ylab="ENaC activity ( N * Po)",
	col=c('steelblue','steelblue1','steelblue','steelblue1'),
	cex.axis=2.5,cex.lab=3
	)
text(bp,joint[3,], labels = format(round(joint[3,],2), 4),pos = 3, cex = 3)
text(c(.7,1.95,3.15,4.35), par("usr")[3] -0.2, labels = as.character(name_perturb), srt = 0, pos = 1, xpd = TRUE,cex=2)
par(mfrow=c(1,1))

# dev.off()





bp=barplot(joint[3,],
	beside = TRUE,ylim=c(0,max(joint[3,])*1.2),
	ylab="ENaC activity ( N * Po)",
	col=c('steelblue','steelblue1','steelblue','steelblue1'),
	cex.axis=2.5,cex.lab=3
	)
text(bp,joint[3,], labels = format(round(joint[3,],2), 4),pos = 3, cex = 3)
text(c(.7,1.95,3.15,4.35), par("usr")[3] -0.2, labels = as.character(name_perturb), srt = 0, pos = 1, xpd = TRUE,cex=2)
par(mfrow=c(1,1))

# dev.off()


#############################################################################


### PLC activation / DGK perturbation ###

### CF ###
tmax=6000
tstep=1
t=seq(0,tmax+1,tstep) 	# time

state[names(state)=='ENaC'] = 80
state[names(state)=='ASL'] = 4
parameters[names(parameters)=='SPLUNC1'] = 0
parameters[names(parameters)=='gamma_CFTR_ASL'] = 0

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,2000,3000,4000,5000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','DGK','DGK','PLC','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	5*1.2,
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)

colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

par(mfrow=c(3,1))
# time plot
plot(out[,1],out[,2],type='l',col=1,ylim=c(0,max(out[,2:9])),main='Time course',xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,25000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,3000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))
	
### ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15:16)) 
	{
	plot(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],ENaC_op(out[,7]),type='l',col='blue',main='ENaC open probability')   

time_sample=c(990,1990,3990,4990)
N_CF_sample = out[time_sample,15]
ASL_CF_sample = out[time_sample,16]
op_CF_sample = ENaC_op(out[time_sample,7])
activity_CF_sample = N_CF_sample*op_CF_sample
almaca_plot = c(almaca_plot,activity_CF_sample[3:4])

### P2 ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15)) 
	{
	plot(out[1:6000,1],out[1:6000,i],
	type='l',lwd=5,col=varcolors[i],
	xlab='time (min)', 
	ylab='Number of actived ENaC channels'
	)
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[1:6000,1],ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='blue',
	xlab='time (min)',	
	ylab='ENaC open probability')   

plot(out[1:6000,1],out[1:6000,i]*ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='black',
	xlab='time (min)',	
	ylab='ENaC activity')   

N_sample = c(N_CF_sample)
ASL_sample = c(ASL_CF_sample)
op_sample = c(op_CF_sample)
activity_sample = c(activity_CF_sample)
name_perturb = c(	'PLC + / DGK +',
			'PLC + / DGK -',
			'PLC act / DGK +',
			'PLC act / DGK -'
			)
cbind(time_sample,name_perturb,N_sample,op_sample,activity_sample)
joint = rbind(N_sample,op_sample,activity_sample)


# tiff("P2F2_Fig_almaca_c.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(2,2))
#windows()
bp=barplot(joint[1,],
	beside = TRUE,ylim=c(0,max(joint[1,])*1.2),
	ylab="Number of actived ENaC channels (N)",
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[1,], labels = format(round(joint[1,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
#windows()
bp=barplot(joint[2,],
	beside = TRUE,ylim=c(0,1.1),
	ylab="ENaC open probability (Po)",
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[2,], labels = format(round(joint[2,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
#windows()
bp=barplot(joint[3,],
	beside = TRUE,ylim=c(0,max(joint[3,])*1.2),
	ylab="ENaC activity ( N * Po)",
	col=c('steelblue','steelblue1','steelblue','steelblue1'),
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[3,], labels = format(round(joint[3,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
par(mfrow=c(1,1))

# dev.off()





#############################################################################


### PI3KI / DGK perturbation ###

### CF ###
tmax=6000
tstep=1
t=seq(0,tmax+1,tstep) 	# time

state[names(state)=='ENaC'] = 80
state[names(state)=='ASL'] = 4
parameters[names(parameters)=='SPLUNC1'] = 0
parameters[names(parameters)=='gamma_CFTR_ASL'] = 0

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,2000,3000,4000,5000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','DGK','DGK','pi_3KI','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	1000 * .25,	
	parameters[names(parameters)=='DGK']*0.75,
	parameters[names(parameters)=='DGK'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)

colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

par(mfrow=c(3,1))
# time plot
plot(out[,1],out[,2],type='l',col=1,ylim=c(0,max(out[,2:9])),main='Time course',xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,25000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,3000),xlab='Time (m)',ylab='Molecule count')
for (i in c(2:9,12,13,14)) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))
	
### ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15:16)) 
	{
	plot(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],ENaC_op(out[,7]),type='l',col='blue',main='ENaC open probability')   

### P2 ENaC plots ###			
windows()
par(mfrow=c(2,2))
# time plot
for (i in c(15)) 
	{
	plot(out[1:6000,1],out[1:6000,i],
	type='l',lwd=5,col=varcolors[i],
	xlab='time (min)', 
	ylab='Number of actived ENaC channels'
	)
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[1:6000,1],ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='blue',
	xlab='time (min)',	
	ylab='ENaC open probability')   

plot(out[1:6000,1],out[1:6000,i]*ENaC_op(out[1:6000,7]),
	type='l',lwd=5,col='black',
	xlab='time (min)',	
	ylab='ENaC activity')   

time_sample = c(990,1990,3990,4990)
N_sample = out[time_sample,15]
op_sample = ENaC_op(out[time_sample,7])
activity_sample = N_sample*op_sample
almaca_plot = c(almaca_plot,activity_sample[3:4])
name_perturb = c(	'PI3KI + / DGK +',
			'PI3KI + / DGK -',
			'PI3KI - / DGK +',
			'PI3KI - / DGK -'
			)
cbind(time_sample,name_perturb,N_sample,op_sample,activity_sample)
joint = rbind(N_sample,op_sample,activity_sample)


# tiff("P2F2_Fig_almaca_d.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(2,2))
#windows()
bp=barplot(joint[1,],
	beside = TRUE,ylim=c(0,max(joint[1,])*1.2),
	ylab="Number of actived ENaC channels (N)",
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[1,], labels = format(round(joint[1,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
#windows()
bp=barplot(joint[2,],
	beside = TRUE,ylim=c(0,1.1),
	ylab="ENaC open probability (Po)",
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[2,], labels = format(round(joint[2,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
#windows()
bp=barplot(joint[3,],
	beside = TRUE,ylim=c(0,max(joint[3,])*1.2),
	ylab="ENaC activity ( N * Po)",
	col=c('steelblue','steelblue1','steelblue','steelblue1'),
	cex.axis=1,cex.lab=1.5
	)
text(bp,joint[3,], labels = format(round(joint[3,],2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=1, cex.axis=.9)
par(mfrow=c(1,1))

# dev.off()







### Figure 4 ###

graphics.off() 	#close plot windows

# model data for Almaca plot
# if you no not want to run everything # almaca_plot = c(7.943916, 4.652055, 17.884668, 9.822216, 65.230149, 65.139804, 9.731012, 5.509128, 17.865020, 9.798913)
almaca_plot = almaca_plot * 100 / almaca_plot[3]
almaca_plot_redux = almaca_plot[1:8]

# almaca data

wd = getwd()
setwd(file.path(getwd(),'3_Fig4_data_Almaca'))
almaca_data = read.csv("ENaC_act_vs_DGK.csv")
setwd(wd)
almaca_data_redux = almaca_data[1:8,]

name_perturb = c(	'WT *',
			'WT / DGK - *',
			'CF',
			'CF / DGK-',
			'CF / PLC-',
			'CF / PLC- / DGK-',
			'CF / PLC+',
			'CF / PLC+ / DGK-'
			)

# the plot
windows()
par(mar = c(12, 4, 2, 2) + 0.2)
bp=barplot(almaca_data_redux[,3],
	beside = TRUE,ylim=c(0,max(almaca_data[,3])*1.2),
	ylab="ENaC activity %",
	col=c('gray70','gray90','gray70','gray90'),
	cex.axis=1.5,cex.lab=1.5
	)

# intervals for Almaca data
segments(x0=bp[,1]-.2,y0=almaca_data_redux[,4],
	x1=bp[,1]+.2,y1=almaca_data_redux[,4],
	col='black',lwd=1.5) 

segments(x0=bp[,1],y0=almaca_data_redux[,2],
	x1=bp[,1],y1=almaca_data_redux[,4],
	col='black',lwd=1.5) 

segments(x0=bp[,1]-.2,y0=almaca_data_redux[,2],
	x1=bp[,1]+.2,y1=almaca_data_redux[,2],
	col='black',lwd=1.5) 

# model results
segments(x0=bp[,1]-.2,y0=almaca_plot_redux,
	x1=bp[,1]+.2,y1=almaca_plot_redux,
	col='blue',lwd=5)

text(bp,almaca_plot_redux, labels = format(round(almaca_plot_redux,2), 4),pos = 3, cex = 1)
axis(1, at=bp,labels=as.character(name_perturb), las=2, cex.axis=1.5)

windows()
# the plot for the really big values
par(mar = c(12, 4, 2, 2) + 0.2)
bp=barplot(almaca_data_redux[,3],
	beside = TRUE,ylim=c(300,max(almaca_plot_redux)*1.2),
	ylab=" ",
	col=c('gray70','gray90','gray70','gray90'),
	cex.axis=1.5,cex.lab=1.5
	)

# model results
segments(x0=bp[,1]-.2,y0=almaca_plot_redux,
	x1=bp[,1]+.2,y1=almaca_plot_redux,
	col='blue',lwd=5)

text(bp,almaca_plot, labels = format(round(almaca_plot,2), 4),pos = 3, cex = 1)


