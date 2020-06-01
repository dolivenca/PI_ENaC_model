### plot for PI45P2 recovery after PLC activation ###

### Run model in apical configuration first ###



tmax=10000
tstep=1
t=seq(0,tmax+1,tstep) 	# time

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,8000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='PLC']*0.01,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='PLC']*1.5,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='PLC']*2,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='PLC']*3,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='PLC']*4,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='PLC']*5,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='DGK']*0.25,
	parameters[names(parameters)=='DGK'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations




out = Cruncher(state,t,equations,parameters,perturbations)
# edit(out)
# View(out)


# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
#tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
#head(out_norm)
# edit(out_norm)



windows()
par(mar=c(5,5,5,3))
plot(out_norm[,1],out_norm[,7],
	type='l',col=varcolors[7],lwd=5,
	xlab='Time (min)',ylab='Fold change',
	main=expression('PI(4,5)P'[2]*' with PLC activation'),
	cex.lab=2.5,cex.axis=2.5,cex.main=2.5
	)

# tiff("P2F2_Fig_S4.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)
par(mar=c(5,5,5,3))
plot(0:700,out_norm[2900:3600,7],
	type='l',col=varcolors[7],lwd=5,
	ylim = c(.67,1.03),
	xlab='Time (min)',ylab='Fold change',
	main=expression('PI(4,5)P'[2]*' with PLC activation'),
	cex.lab=2.5,cex.axis=2.5,cex.main=2.5
	)
polygon (x=c(100,600,600,100),y=c(0,0,2,2),col='lightpink',border=NA)
points(0:700,out_norm[2900:3600,7],
	type='l',col=varcolors[7],lwd=5,
	ylim = c(.67,1.03),
	xlab='Time (min)',ylab='Fold change',
	main=expression('PI(4,5)P'[2]*' with PLC activation'),
	cex.lab=2.5,cex.axis=2.5,cex.main=2.5
	)
legend(x=100,y=1.05,legend='PLC activation',bty = "n",cex=1.5)
#legend(x=2700,y=.7,legend='PI(4,5)P2 recovery',bty = "n",cex=1.5)
# dev.off()



windows()
# tiff("P2F2_Fig_S4.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)
par(mar=c(5,5,5,3))
plot(0:700,out_norm[2900:3600,7],
	type='l',col=varcolors[7],lwd=5,
	xlim=c(580,620),ylim = c(.67,1.03),
	xlab='Time (min)',ylab='Fold change',
	# main=expression('PI(4,5)P'[2]*' with PLC activation'),
	cex.lab=2.5,cex.axis=2.5,cex.main=2.5
	)
polygon (x=c(100,600,600,100),y=c(0,0,2,2),col='lightpink',border=NA)
points(0:700,out_norm[2900:3600,7],
	type='l',col=varcolors[7],lwd=5,
	xlim=c(580,620),ylim = c(.67,1.03),
	xlab='Time (min)',ylab='Fold change',
	# main=expression('PI(4,5)P'[2]*' with PLC activation'),
	cex.lab=2.5,cex.axis=2.5,cex.main=2.5
	)
legend(x=580,y=1.05,legend='PLC activation',bty = "n",cex=1.5)
#legend(x=2700,y=.7,legend='PI(4,5)P2 recovery',bty = "n",cex=1.5)
# dev.off()



