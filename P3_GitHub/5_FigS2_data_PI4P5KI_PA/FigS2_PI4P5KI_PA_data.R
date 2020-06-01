

graphics.off() 	#close plot windows

############################
###				 ###
### Pardo_PA_PI45P2_data ###
###				 ###
############################


wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_data = read.csv("Pardo_PI4P5KI_PA_data.csv")
setwd(wd)

plot(Pardo_data$PA,Pardo_data$PI45P2,type='p',col='black',pch=19,xlim=c(0,.15),ylim=c(0,100),main='Pardo PA PI45P2 data',xlab='PA uM',ylab='PI45P2 nmol/min/mg of protein')

wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_model = read.csv("Pardo_PI4P5KI_PA_model.csv") 
setwd(wd)
points(Pardo_model$PA,Pardo_model$PI45P2,type='l',col='black')


#########################################################################
# V4_45 = gamma4_45 * 
#		PI4P^f4_45_PI4P * 
#		pip_5KI * 
#		65.643 * 1.66 * 10^-18 
#		* HS(PA) * 
#		PI45P2^f4_45_PI45P2
#
# gamma4_45 = 8.493641e+15, 
# f4_45_PI4P = 4.596175e-02, 
# f4_45_PA = 0.2, 
# f4_45_PI45P2 = -0.05,
#
# 65 uM of PI4P = 39143 molecules/um^3 of PI4P and no initial PI45P2
# 8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 
#
#	* 10000^(-.05)
#
##########################################################################


#################################################################
#
# values for activity of PI4P5KI without PA
#
# for 100 uM = 60221 molecules/um^3 of PI4P -> 9 nmol/min/mg = nmol_molecules(9) = .6 molecules/min/per unit of protein 
# for 150 uM = 90332 molecules/um^3 of PI4P -> 16 nmol/min/mg = nmol_molecules(16) = 1.05 molecules/min/per unit of protein 
# for 65 uM = 39143 molecules/um^3 of PI4P -> 3 nmol/min/mg = nmol_molecules(3) = .2 molecules/min/per unit of protein 
#
# it is because of this that all y's on the next plots have a +.2
# because I am considering that PI4P5KI basal activity is of about
# .2 molecules/min/per unit of PI4P5KI protein
#
##################################################################




### Functin to calculate efective surface consentration of PA ###
PA = seq(0.5,260,.5)	# uM
PI4P = 65	# uM
TritonX100 =  1600 # uM | should be, according to Luis 1.6 mM = 1600 uM = 963542.64 molecules/um^3 but there is room for doubt 
PAesc = function(PA) {(PI4P * PA)/(PA + PI4P + TritonX100)}
# plot(PA,PAesc(PA))

### Function to transform uM into molecules / um3 ### 
uM_molecules_um3 = function (uM) {uM * 10^(-6) * 1000 * 6 * 10^23 * 10^(-18)} 
uM_molecules_um3(200)

### Function to transform nmol/min/mg -> molecules/min/unit of PI4P5KI protein ###
#20 nmol/min/mg = 1.3 molecules/min/per unit of PI4P5KI protein 
nmol_molecules = function(nmol) {nmol * 10^(-9) * 6 * 10^23 *65.643 *1.66*10^-18}
nmol_molecules(20)

### Pardo_PA_PI45P2_data - my units ###

wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_data = read.csv("Pardo_PI4P5KI_PA_data.csv")
setwd(wd)

plot(
	uM_molecules_um3(
		Pardo_data$PA*(PI4P + TritonX100) / (PI4P - Pardo_data$PA) 
		),
	(nmol_molecules(Pardo_data$PI45P2)+.2)/.2,
	type='p',col='black',pch=19,
	xlim=c(0,16000),ylim=c(0,30),
	main='Pardo PA PI45P2 data',
	xlab='PA molecules/um^3',
	ylab='PI45P2 nmol/min/mg of protein'
	)


wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_model = read.csv("Pardo_PI4P5KI_PA_model.csv") 
setwd(wd)

points(
	uM_molecules_um3(
		Pardo_model$PA*(PI4P + TritonX100) / (PI4P - Pardo_model$PA) 
		),
	(nmol_molecules(Pardo_model$PI45P2)+.2)/.2,
	type='l',col='black')

### me trying to recreate Pardo's hill function ###
Vmax = 100
Km = .055
points(
	uM_molecules_um3(PA),
	(nmol_molecules((Vmax*PAesc(PA)^2.3)/(Km^2.3+PAesc(PA)^2.3))+.2)/.2,
	type='l',col='cyan')


### My hypersensitivity function ###
HS = function (PA)(1 + 2*PA^2/(PA^2+10000^2))
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')






########################################
### trying to correct for TritonX100 ###
##############################################################################
#
# To get a good fit, I have to assume that:
# V - TritonX100 is decreasing the basal activity of PI4P5KI 
# 	from 1.5 to .2 molecules/min/unit of protein
# V - the maximum activation of PI4P5KI by PA is 2.2 as reported by 
#	Jones et al. and not almost 6 as reported by Pardo. 
#     Therefor I transform the y's by * 2.2/6.
# V - assuming that TritonX100 total consentration is 15000 uM,
#	PA concentrations are OK.
#
##############################################################################

### Pardo_PA_PI45P2_data - my units ###
wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_data = read.csv("Pardo_PI4P5KI_PA_data.csv")
setwd(wd)

plot(
	uM_molecules_um3(
		(Pardo_data$PA*(PI4P + TritonX100) / (PI4P - Pardo_data$PA))
		),
	(nmol_molecules(Pardo_data$PI45P2)/.2)*3/30+1.5,
	type='p',col='black',pch=19,
	xlim=c(0,16000),ylim=c(0,6),
	main='Pardo PA PI45P2 data',
	xlab='PA molecules/um^3',
	ylab='PI45P2 nmol/min/mg of protein'
	)

wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_model = read.csv("Pardo_PI4P5KI_PA_model.csv")
setwd(wd)

points(
	uM_molecules_um3(
		(Pardo_model$PA*(PI4P + TritonX100) / (PI4P - Pardo_model$PA)) 
		),
	(nmol_molecules(Pardo_model$PI45P2)/.2)*3/30+1.5,
	type='l',col='black')

### me trying to recreate Pardo's hill function ###
Vmax = 100
Km = .055
points(
	uM_molecules_um3(PA),
	(nmol_molecules((Vmax*PAesc(PA)^2.3)/(Km^2.3+PAesc(PA)^2.3))/.2)*3/30+1.5,
	type='l',col='cyan')


### My hypersensitivity function ###
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')





#########################
###			    ###
###	Moritz data	    ###
###			    ###
#########################
wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Moritz_data = read.csv("Moritz_PI4P5KI_PA_data.csv")
setwd(wd)

plot(Moritz_data$PA_PIP,Moritz_data$fold_change,
	type='p',col='black',pch=19,
	xlim=c(0,2),ylim=c(0,20),
	main='Pardo PA PI45P2 data',
	xlab='PA/PIP',ylab='PI45P2 nmol/min/mg of protein')
 
### Moritz data in my units ###

plot(
	Moritz_data$PA_PIP*48000,
	Moritz_data$fold_change*1.5,
	type='p',col='black',
	pch=19,xlim=c(0,100000),ylim=c(0,30),
	main='Moritz PA PI45P2 data',
	xlab='PA molecules/um^2',ylab='PI45P2 nmol/min/mg of protein'
	)

### My hypersensitivity function ###
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')


#######################################
### trying to correct of TritonX100 ###
#######################################

plot(
	Moritz_data$PA_PIP*48000,
	Moritz_data$fold_change*3*1.5/30+1.5,
	type='p',col='black',
	pch=17,xlim=c(0,100000),ylim=c(0,5),
	main='Moritz PA PI45P2 data',
	xlab='PA molecules/um^2',ylab='PI45P2 nmol/min/mg of protein'
	)

### My hypersensitivity function ###
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')





#########################
###			    ###
###	Jenkins data    ###
###			    ###
#########################
wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Jenkins_data = read.csv("Jenkins_PI4P5KI_PA_data.csv")
setwd(wd)

plot(Jenkins_data$PA,Jenkins_data$fold_change,
	type='p',col='black',pch=19,
	xlim=c(0,170),ylim=c(0,13),
	main='Jenkins PA PI45P2 data',xlab='PA uM',ylab='Fold Change')
 
### Jenkins data in my units ###

plot(
	uM_molecules_um3(Jenkins_data$PA),
	Jenkins_data$fold_change*1.5,
	type='p',col='black',
	pch=19,xlim=c(0,100000),ylim=c(0,30),
	main='Jenkins PA PI45P2 data',
	xlab='PA molecules/um^3',ylab='PI45P2 nmol/min/mg of protein'
	)

### My hypersensitivity function ###
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')

#######################################
### trying to correct of TritonX100 ###
#######################################

plot(
	uM_molecules_um3(Jenkins_data$PA),
	Jenkins_data$fold_change*3*1.5/13 - 0.1307692 + 1.5,
	type='p',col='black',
	pch=15,xlim=c(0,100000),ylim=c(0,5),
	main='Jenkins PA PI45P2 data',
	xlab='PA molecules/um^2',ylab='PI45P2 nmol/min/mg of protein'
	)

### My hypersensitivity function ###
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')







###########################
### put all in one plot ###
###########################

### quantities for PI4P and TritonX100 ###
PI4P = 65	# uM
TritonX100 =  1600 # uM

### Pardo_PA_PI45P2_data - my units ###
wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_data = read.csv("Pardo_PI4P5KI_PA_data.csv")
setwd(wd)

plot(
	uM_molecules_um3(
		(Pardo_data$PA*(PI4P + TritonX100) / (PI4P - Pardo_data$PA)) 
		),
	(nmol_molecules(Pardo_data$PI45P2)/.2)*3/30+1.5,
	type='p',col='black',pch=19,
	xlim=c(0,100000),ylim=c(0,6),
	main='Pardo PA PI45P2 data',
	xlab='PA molecules/um^3',
	ylab='PI45P2 nmol/min/mg of protein'
	)

wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_model = read.csv("Pardo_PI4P5KI_PA_model.csv")
setwd(wd)

points(
	uM_molecules_um3(
		(Pardo_model$PA*(PI4P + TritonX100) / (PI4P - Pardo_model$PA)) 
		),
	(nmol_molecules(Pardo_model$PI45P2)/.2)*3/30+1.5,
	type='l',col='black')

### me trying to recreate Pardo's hill function ###
Vmax = 100
Km = .055
points(
	uM_molecules_um3(PA),
	(nmol_molecules((Vmax*PAesc(PA)^2.3)/(Km^2.3+PAesc(PA)^2.3))/.2)*3/30+1.5,
	type='l',col='cyan')

### Moritz data ###
points(
	Moritz_data$PA_PIP*48000,
	Moritz_data$fold_change*3*1.5/30+1.5,
	type='p',col='black',
	pch=17,xlim=c(0,20000),ylim=c(0,5),
	main='Moritz PA PI45P2 data',
	xlab='PA molecules/um^2',ylab='PI45P2 molecules/min/unit of PI4P5KI'
	)

### Jerkins data ###
points(
	uM_molecules_um3(Jenkins_data$PA),
	Jenkins_data$fold_change*2.2*1.5/13 - 0.1307692 + 1.5,
	type='p',col='black',
	pch=15,xlim=c(0,100000),ylim=c(0,5),
	main='Jenkins PA PI45P2 data',
	xlab='PA molecules/um^2',ylab='PI45P2 nmol/min/mg of protein'
	)

### My hypersensitivity function ###
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')

legend(
	x=60000, y=1.75,
	legend=c('Moritz`s data','Jarquim-Pardo`s data','Jenkins data'),
	pch=c(17,19,15),
	bty = "n"
	)
legend(
	x=60000, y=1,
	lty=1,
	col=c('black','blue'),
	legend=c('Jarquim-Pardo`s model','Our sigmoid function'),
	text.col=c('black','blue'),
	bty = "n"
	)




###########################################################################
###												###
###	First - adjust Pardo's TritonX100 to better aproach Moritz data	###
###	find the best Hill that fit all the data					###
###												###
###########################################################################
### Pardo's data ###
### quantities for PI4P and TritonX100 ###
PI4P = 65	# uM
TritonX100 =  10*1600 # uM

# tiff("S2.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

### Pardo_PA_PI45P2_data - tritonX100 = 35000 uM ###
wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_data = read.csv("Pardo_PI4P5KI_PA_data.csv")
setwd(wd)

plot(
	uM_molecules_um3(
		(Pardo_data$PA*(PI4P + TritonX100) / (PI4P - Pardo_data$PA)) 
		),
	(nmol_molecules(Pardo_data$PI45P2)/.2)*3/30+1.5,
	type='p',col='black',pch=19,
	xlim=c(0,100000),ylim=c(0,6),
	main='PIP5KI activation by PA',
#	main=expression('PI(4,5)P'[2]*' activation by PA data'),
	xlab=expression('PA molecules/'~ �m^{2}),
	ylab=expression('PI(4,5)P'[2]*' molecules/min/mg of protein')
	)
### Pardo's model ###
wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_model = read.csv("Pardo_PI4P5KI_PA_model.csv") 
setwd(wd)

points(
	uM_molecules_um3(
		(Pardo_model$PA*(PI4P + TritonX100) / (PI4P - Pardo_model$PA)) 
		),
	(nmol_molecules(Pardo_model$PI45P2)/.2)*3/30+1.5,
	type='l',col='black')

### Pardo's data - TritonX100 = 1600 uM ###
wd = getwd()
setwd(file.path(getwd(),'5_FigS2_data_PI4P5KI_PA'))
Pardo_data = read.csv("Pardo_PI4P5KI_PA_data.csv")
setwd(wd)

points(
	uM_molecules_um3(
		(Pardo_data$PA*(PI4P + 1600) / (PI4P - Pardo_data$PA)) 
		),
	(nmol_molecules(Pardo_data$PI45P2)/.2)*3/30+1.5,
	type='p',col='black',pch=1,
	xlim=c(0,16000),ylim=c(0,6),
	main='Pardo PA PI45P2 data',
	xlab='PA molecules/um^3',
	ylab='PI45P2 nmol/min/mg of protein'
	)


### Moritz data ###
points(
	Moritz_data$PA_PIP*48000,
	Moritz_data$fold_change*3*1.5/30+1.5,
	type='p',col='black',
	pch=15,xlim=c(0,20000),ylim=c(0,5),
	main='Moritz PA PI45P2 data',
	xlab='PA molecules/um^2',ylab='PI45P2 molecules/min/unit of PI4P5KI'
	)

### Jerkins data ###
points(
	uM_molecules_um3(Jenkins_data$PA),
	Jenkins_data$fold_change*2.2*1.5/13 - 0.1307692 + 1.5,
	type='p',col='black',
	pch=17,xlim=c(0,100000),ylim=c(0,5),
	main='Jenkins PA PI(4,5)P2 data',
	xlab=expression('PA molecules/'~ �m^{2}),ylab='PI(4,5)P2 nmol/min/mg of protein'
	)

### My hypersensitivity function ###
points(
	uM_molecules_um3(PA),
	8.493641e+15 * 40000^4.596175e-02 * 1 * 65.643 *1.66*10^-18 * 1 * HS(uM_molecules_um3(PA)),
	type='l',col='blue')

legend(
	x=20000, y=2,
	legend=c('Moritz`s data','Jenkins data','Corrected Jarquim-Pardo`s data (TritonX100 = 16000)','Raw Jarquim-Pardo`s data (TritonX100 = 1600)'),
	text.col=c('black','black','black','black'),
	pch=c(17,15,19,1),
	col=c('black','black','black','black'),
	bty = "n"
	)
legend(
	x=14000, y=1,
	lty=1,
	col=c('black','blue'),
	legend=c('Corrected Jarquim-Pardo`s model (TritonX100 = 16000)','Our sigmoid function'),
	text.col=c('black','blue'),
	bty = "n"
	)

# dev.off()









