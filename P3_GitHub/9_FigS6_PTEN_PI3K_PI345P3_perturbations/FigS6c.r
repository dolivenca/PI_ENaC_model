### Please run the model first in the basolateral configuration (MK20_ENaC_7_basolateral.R) ### 
# or alter the model to simulate basolateral part of the plasma membrane
# To do that, you must divide gammaPTENc_PTENa by 10
# This will decrease the PTEN ability to attach to the plasma membrane (became active)   



### Pochynyuk_data_PIP3_PI3KI ###

# tiff("P2F2_Fig_3_pochy.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

wd = getwd()
setwd(file.path(wd,'9_FigS6_PTEN_PI3K_PI345P3_perturbations'))
pochynyuk_data_2 = read.csv("Pochynyuk_data_PI345P3_PI3KI.csv")
setwd(wd)

par(mar=c(5,6,4,5)+.1)
plot(pochynyuk_data_2$x,pochynyuk_data_2$y,type='b',
	xlim=c(-10,25),ylim=c(0,2),
	main='PI3KI inhibition',xlab='Time (min)',
	ylab='PI345P3 levels (fold change)',
	cex.lab=3,cex.main=3,cex.axis=2	
	)
# legend(x=0,y=2,legend=c('PI345P3 dynamics after PI3KI inhibition'),bty = "n")



# Insert the number of time units that the simulation will run
tmax=3000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

### perturbations for Pochynyuk PI345P3_PI3KI data 

# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,2000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','pi_3KI_a','novar')
# On the third the new value for the variable.
pval=c(NA,35,NA)
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
# edit(out_norm)


points((out_norm[1990:2060,1]-2000),out_norm[1990:2060,9],type='l',col='green',lwd='5')

# dev.off()







