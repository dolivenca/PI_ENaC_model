#####################################
#						#
#	Local Sensitive analysis	#
#						#
#		For all models		#
#						#
#####################################

### load the apical model steady state (state), parameters and equations  ### 

Ss_old=state	
# for sens analysis alter the state values for ENaC and ASL. 
# They are close but the difference will mess up their sensitivity.
state[14:15] = c(3.553597e+01, 7.192864e+00)
Ss_old[14:15] = c(3.553597e+01, 7.192864e+00)
										# storing steady state to be studied 
Sens=matrix(rep(NA,length(parameters)*length(Ss_old)),ncol=length(Ss_old))	# creating response matrix
colnames(Sens)=names(Ss_old)									# naming the columns after the dependent variables
rownames(Sens)=names(parameters)								# naming the rows after the parameters
Sens

t=seq(0,5000+1,1) 	# time

for (i in 1:length(parameters))				# cycling the parameters
	{
	Pars=parameters						# creating discardable parameter vector to be changed
	Pars[i]=parameters[i]*(1.01)				# changing one parameter in the discardable parameter vector in 1%
	out_new <- ode(y = state, times = t, func = equations,parms=Pars)		#runing the ode solver
	Ss_new=tail(out_new,1)[,2:(length(Ss_old)+1)]	# new steady state
	tsens=vector()						# creating sensitivity vector for this turn 
	for (ii in 1:length(Ss_old))						# cycling all dependent variables
		{
		tsens=c(tsens,(Ss_new[ii]-Ss_old[ii])*100/Ss_old[ii])	#calculating sensitivities
		}
	Sens[i,]=tsens						# puting this parameter sensitivities in sentitivity matrix
	}
Sens									# show all sensitivity matrix
### please note that Sens have gamma6 and gamma4 in different positions relative to the excel file and the supplement table S3. ###



# sum(abs(Sens))							# sum of all sensitivities absolute value
# sum(abs(Sens)>=1)							# amount of significant sensibilities (greater then 1%)
# (Sens_show=abs(Sens)>=1)					# mark significant sensibilities 



# for(iii in 1:length(parameters))				# cycling the rows of the sensitivity matrix
#	{
#	for(iv in 1:length(Ss_old))				# cycling the columns of the sensitivity matrix
#		{
#		if (Sens_show[iii,iv]==1){Sens_show[iii,iv]=Sens[iii,iv]}	# substitute the TRUEs for the significant sensibilities
#		}
#	}
# Sens_show								# show sensitivities matrix with just the significant sensibilities
# sum(abs(Sens_show[Sens_show!=0]))				# sum of the significant sensitivities absolute value

# dim(Sens_show)							# dimentions of sensitivity matrix
# dim(Sens_show)[1]*dim(Sens_show)[2]				# total sensitivities accessed 
# sum(abs(Sens)>=1)/(dim(Sens_show)[1]*dim(Sens_show)[2]) # percentage of high sensitivities



### save and read from disk ###
# as a variable in memory
# save(Sens,Sens_show,file = "Local_sens_analysis_for_ENaC_model.RData")	
# load( "Local_sens_analysis_for_ENaC_model.RData")



