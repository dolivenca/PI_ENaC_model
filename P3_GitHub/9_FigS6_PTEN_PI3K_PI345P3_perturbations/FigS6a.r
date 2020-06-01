### Please run the model first in the basolateral configuration (MK20_ENaC_7_basolateral.R) ### 
# or alter the model to simulate basolateral part of the plasma membrane
# To do that, you must divide gammaPTENc_PTENa by 10
# This will decrease the PTEN ability to attach to the plasma membrane (became active)   



### Carsten data - second try ###
# first plot - original values
wd = getwd()
setwd(file.path(wd,'9_FigS6_PTEN_PI3K_PI345P3_perturbations'))
shultz_data_1_PI3KI = read.csv("Shultz_data_1_PI3KI.csv")
setwd(wd)

windows()
plot(shultz_data_1_PI3KI,type='p',col='firebrick',xlim=c(0,50),ylim=c(.6,1.5),main='Carsten data 1')
points(shultz_data_1_PI3KI$x,shultz_data_1_PI3KI$y-0.05,type='p',col='firebrick1',pch='-')
points(shultz_data_1_PI3KI$x,shultz_data_1_PI3KI$y+0.05,type='p',col='firebrick1',pch='-') 
segments(x0=shultz_data_1_PI3KI$x,y0=shultz_data_1_PI3KI$y-0.05,y1=shultz_data_1_PI3KI$y+0.05,col='firebrick1')


wd = getwd()
setwd(file.path(wd,'9_FigS6_PTEN_PI3K_PI345P3_perturbations'))
shultz_data_1_PI345P3 = read.csv("Shultz_data_1_PI345P3.csv")
setwd(wd)

points(shultz_data_1_PI345P3,type='p',col='green')
points(shultz_data_1_PI345P3$x,shultz_data_1_PI345P3$y-0.05,type='p',col='lightgreen',pch='-')
points(shultz_data_1_PI345P3$x,shultz_data_1_PI345P3$y+0.05,type='p',col='lightgreen',pch='-') 
segments(x0=shultz_data_1_PI345P3$x,y0=shultz_data_1_PI345P3$y-0.05,y1=shultz_data_1_PI345P3$y+0.05,col="lightgreen")

legend(x=0,y=1.5,legend=c('PI345P3','PI3K'),lty=1,col=c('green','firebrick'),text.col=c('green','firebrick'),bty = "n")



### Carsten data - second try ###
# first plot - hammered values so PI3KI and PIP3 start at 100%
wd = getwd()
setwd(file.path(wd,'9_FigS6_PTEN_PI3K_PI345P3_perturbations'))
shultz_data_1_PI3KI = read.csv("Shultz_data_1_PI3KI.csv")
setwd(wd)
shultz_data_1_PI3KI$y=shultz_data_1_PI3KI$y+(1-shultz_data_1_PI3KI$y[1])

# windows()
# tiff("P2F2_Fig_3_Carsten1.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)
par(mar=c(5,6,4,5)+.1)
plot(shultz_data_1_PI3KI,type='p',
	col='firebrick',xlim=c(0,50),ylim=c(.6,1.7),
	main='PI3K activation',xlab='Time (min)',ylab='fold change',
	cex.main=3,cex.lab=3,cex.axis=2
	)
points(shultz_data_1_PI3KI$x,shultz_data_1_PI3KI$y-0.05,type='p',col='firebrick1',pch='-')
points(shultz_data_1_PI3KI$x,shultz_data_1_PI3KI$y+0.05,type='p',col='firebrick1',pch='-') 
segments(x0=shultz_data_1_PI3KI$x,y0=shultz_data_1_PI3KI$y-0.05,y1=shultz_data_1_PI3KI$y+0.05,col='firebrick1')

wd = getwd()
setwd(file.path(wd,'9_FigS6_PTEN_PI3K_PI345P3_perturbations'))
shultz_data_1_PI345P3 = read.csv("Shultz_data_1_PI345P3.csv")
setwd(wd)
shultz_data_1_PI345P3$y=shultz_data_1_PI345P3$y+(1-shultz_data_1_PI345P3$y[1])

points(shultz_data_1_PI345P3,type='p',col='green')
points(shultz_data_1_PI345P3$x,shultz_data_1_PI345P3$y-0.05,type='p',col='lightgreen',pch='-')
points(shultz_data_1_PI345P3$x,shultz_data_1_PI345P3$y+0.05,type='p',col='lightgreen',pch='-') 
segments(x0=shultz_data_1_PI345P3$x,y0=shultz_data_1_PI345P3$y-0.05,y1=shultz_data_1_PI345P3$y+0.05,col="lightgreen")

legend(x=0,y=1.6,legend=c('PI345P3','PI3K'),lty=1,col=c('green','firebrick'),text.col=c('green','firebrick'),bty = "n",cex=2)



# Insert the number of time units that the simulation will run
tmax=3000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

### perturbations for Carsten data 1

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,2005,2040,2042,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','gammaPI3KIc_PI3KIa','gammaPI3KIc_PI3KIa','gammaPI3KIc_PI3KIa','pi_3KI_a','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='gammaPI3KIc_PI3KIa'],
	parameters[names(parameters)=='gammaPI3KIc_PI3KIa']*28,	
	parameters[names(parameters)=='gammaPI3KIc_PI3KIa']*1,353.172211855607,
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations
# 353.172211855607
# try 3 for initial model

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

points((out_norm[2000:2060,1]-2000),out_norm[2000:2060,9],type='l',col='green',lwd='5')
points((out_norm[2000:2060,1]-2000),out_norm[2000:2060,10],type='l',col='firebrick',lwd='5')

# dev.off()



