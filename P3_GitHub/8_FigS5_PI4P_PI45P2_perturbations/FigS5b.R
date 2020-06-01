### Run model in apical configuration first ###


### Vernai_4_45_data ###

wd = getwd()
setwd(file.path(wd,'8_FigS5_PI4P_PI45P2_perturbations//data_Varnai'))
Vernai_data_4_phosphatase = read.csv("Varnai_data_4_phosphatase_PI4P_PI45P2.csv")
setwd(wd)

# tiff("P2F2_Fig_2_Vernai_1.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)
windows()
par(mar=c(5,6,4,5)+.1)
plot(Vernai_data_4_phosphatase$X,Vernai_data_4_phosphatase$PI45P2,
	type='p',col='red',pch=19,
	xlim=c(-5,20),ylim=c(0,120),
	main='4 phosphatase activation',
	xlab='Time (min)',ylab='percentage fold change',
	cex=2, cex.main=3, cex.lab=3, cex.axis=2
	)
points(Vernai_data_4_phosphatase$X,Vernai_data_4_phosphatase$PI4P,type='p',col='blue',pch=19,cex=2)
legend(x=1,y=30,legend=c('PI4P','PI45P2'),pch=c(19,19),col=c('blue','red'),text.col=c('blue','red'),bty = "n",cex=2)
legend(x=8,y=30,legend=c('Model PI4P','Model PI45P2'),lty=1,lwd=5,col=c('blue','red'),text.col=c('blue','red'),bty = "n",cex=2)



# Insert the number of time units that the simulation will run
tmax=2000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,1004,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','SAC1','SAC1','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='SAC1']*8,
	parameters[names(parameters)=='SAC1']*2,
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations


out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out)
# View(out)

# normalized graphs
stst=out[150,]						# getting steady state at time 500
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
# View(out_norm)

points((out_norm[995:1020,1]-1000),out_norm[995:1020,7]*100,type='l',col='red',lwd='5')
points((out_norm[995:1020,1]-1000),out_norm[995:1020,4]*100,type='l',col='blue',lwd='5')

# dev.off()



