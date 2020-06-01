### Run model in apical configuration first ###

### Figure S5a ###
### Pochynyuk_data_PIP2_ENaC ###

# tiff("P2F2_Fig_2_pochy.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

wd = getwd()
setwd(file.path(wd,'8_FigS5_PI4P_PI45P2_perturbations/data_Pochynyuck'))
pochynyuk_data_1 = read.csv("Pochynyuk_data_PI45P2_PLC.csv")
setwd(wd)

par(mar=c(5,6,4,5)+.1)
plot(pochynyuk_data_1$PLC_active_X,pochynyuk_data_1$PLC_active_Y,
	type='p',col='firebrick',pch=15,cex=2,xlim=c(0,15),ylim=c(0,2),
	main='PLC perturbations',
	xlab='Time (min)',
	ylab='PI45P2 levels (fold change)',
	cex.lab=3,cex.main=3,cex.axis=2
	)
points(pochynyuk_data_1$PLC_inhibited_X,pochynyuk_data_1$PLC_inhibited_Y,
	type='p',col='blue',pch=17,cex=2)
legend(x=0,y=2,legend=c('PLC inhibited','PLC active'),lty=1,col=c('blue','firebrick'),text.col=c('blue','firebrick'),bty = "n",cex=2)



# Insert the number of time units that the simulation will run
tmax=5000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

### perturbations for Pochynyuk PI45P2_PLC data 

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
#ptime=c(0,1000,1500,2000,2500,3000,4000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
#pvar=c('novar','gamma_PLC_i_a','gamma_PLC_i_a','gamma_PLC_i_a','gamma_PLC_i_a','DGK','DGK','novar')
# On the third the new value for the variable.
#pval=c(NA,
#	parameters[names(parameters)=='gamma_PLC_i_a']*0.01,
#	parameters[names(parameters)=='gamma_PLC_i_a'],
#	parameters[names(parameters)=='gamma_PLC_i_a']*20,
#	parameters[names(parameters)=='gamma_PLC_i_a'],
#	parameters[names(parameters)=='DGK']*0.33,
#	parameters[names(parameters)=='DGK'],
#	NA)
#perturbations=data.frame(time=ptime,var=pvar,val=pval)
#perturbations



# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,1500,2000,2500,3000,4000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','PLC','PLC','PLC','PLC','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='PLC']*0.01,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='PLC']*4,
	parameters[names(parameters)=='PLC'],
	parameters[names(parameters)=='DGK']*0.1,
	parameters[names(parameters)=='DGK'],
	NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations


out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out)
# View(out)

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
# View(out_norm)

points((out_norm[1000:1015,1]-1000),out_norm[1000:1015,7],type='l',col='blue',lwd='5')
points((out_norm[2000:2015,1]-2000),out_norm[2000:2015,7],type='l',col='firebrick',lwd='5')

# dev.off()








