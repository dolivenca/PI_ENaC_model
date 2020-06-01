### Model ###
### This is apical ###

# for apical multiply gammaPTENc_PTENa by 10

### There is a unspected bistability here. Change gamma_5K4P_5K_45P to .075

rm(list=ls())	#clean memory
graphics.off() 	#close plot windows



################## competition and regulation - this is basal ##########################

library(deSolve)	#library for solve diferential equations
#citation("deSolve")

# Insert the number of time units that the simulation will run
tmax=10000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time



HS = function (PA) (1.783743) #(1 + 2*PA^2/(PA^2+10000^2))

ENaC_op = function(PIP2)
	{.02+.8/(1+181429*exp(-.00112*PIP2))}
ENaC_op(10000)

mp2_PI3KI = 0.0006217539	# metaparameter for PI3KI
mp3_PTEN = 0.003			# metaparameter for PTEN




# Inicial states for variables
# insert, in the parentesis, the inicial value of the dependent variables (ex: X=2) 
state=c(
	PI 	   = 	290115.435296034,
	PI3P 	   = 	94.6235111384769,
	PI4P 	   = 	11095.7012602329,
	PI5P     = 	97.176844279141, 
	PI35P2   =	21.3415882859489,
	PI45P2   =	9851.42105324773,
	PI34P2   =	40.3575424621684,
	PI345P3  =	1.12832194419032,
	pi_3KI_a =	5.82680164591602,
	PTEN_a   =	14.5509720205084,
	DAG      =	5787.00059628337,
	IP3      =	499.14011589903,
	PA       =	8025.19773414257,
	ENaC 	   =  35,
	ASL 	   =  7
	) 


# insert, in the parentecis, the parameters and independent variables and their values (ex: const=5, rate=2). 
# If no parameters and independent variables just put a zero.
parameters=c(
		pi_3KI = 1000,			# this and PTEN control the level of PI(345)P3
		PTEN = 50,				# (phosphatase and tensin homologue deleted on chromosome 10) PI_3P family. This and pip_3K control the level of PI(345)P3 

		pi_3KII = 75,	      	# creates PI34P2 from PI4P. In vitro can phosphorilate PI to PI4P so effectly that the original name of this enzyme was PI_4K_type_I 
		pi_3KIII = 5,			# Phosphorilate PI into PI3P.
		pi_4K = 400,		   	# type II. Phosphorilate PI into PI4P. Very similar to PI_3Ks in structure. Type I was after all, PI_3KII 
		pi_Kfyve = 48,			# type III PIP_5K. Really diferent from the other PIP_5K.
		pip_5KI = 330,			# Phosphorilate PI4P into PI45P2. Activated by PA.
		pip_5KII = 845,    		# Phosphorilate PI5P into PI45P2 
		pi_4K_pip_5KI = 230,		# this is here because PI45P2 should not be inflenced by PI4P depletion

		SYNJ = 30, 				# phosphatase that transforms PI45P2 directly into PI
		SAC1 = 100,
		SAC2 = 1,
		SAC3 = 5,
		INPP4 = 50,
		TMEM55 = 20,
		MTMR1_6_14 = 63,
		MTMR78 = 18,

		ORCL1 = 3,
		INPP5BJ = 1,
		INPP5E = 1,
		SKIP = 5,
		SHIP2 = 1,

		PS = 150000,

		PLC = 5,	  		# Phospholipase C. Breaks PI45P2 into DAG and IP3.
		DGK = 50,  			# Diacylglicerol Kinase. Transforms DAG into PA.
		LPP = 100,  		# Phosphatase that transforms PA into DAG, S1P into Shingosine and C1P into ceramide.

		gamma_0 = 15000-787.3, 
		gamma_4 = 150,
		gamma_3 = 10,

		gamma0_3a = 0*5.459875e+14, f0_3a_PI = 1.138439e-01 , 
		gamma0_3b = 8.809591e+14, f0_3b_PI = 0.1761429 , f0_3b_PI4P = -0.07030248,
		gamma0_3c = 5.459875e+14, f0_3c_PI = 1.138439e-01 , 			
		gamma3_0a = 4.557615e+13, f3_0a_PI3P = 0.9995614578, f3_0a_PI4P = -0.04385422, f3_0a_PI5P = -0.0004385422, f3_0a_PI35P2 = -0.0002192711, f3_0a_PI45P2_PI4P = -0.04385422, f3_0a_PI45P2_PI5P = -0.2383381, f3_0a_PI345P3 = -0.001668266, f3_0a_PI45P2_PI = -0.01096355,
		gamma3_0c = 5.096813e+12, f3_0c_PI3P = 0.9993781095, f3_0c_PI4P = -0.06218905, f3_0c_PI5P = -0.0006218905, f3_0c_PI35P2 = -0.0003109453,
		gamma3_0d = 3.059076e+12, f3_0d_PI3P = 9.993364e-01, f3_0d_PI35P2 = -3.317850e-04,
		gamma3_0e = 3.067331e+12, f3_0e_PI3P = 9.993362e-01,


		gamma0_4 = 5.100640e+14, f0_4_PI = 2.864618e-01,  
		gamma4_0a = 4.557615e+13, f4_0a_PI3P = -0.0004385422, f4_0a_PI4P = 0.95614578, f4_0a_PI5P = -0.0004385422, f4_0a_PI35P2 = -0.0002192711, f4_0a_PI45P2_PI4P = -0.04385422, f4_0a_PI45P2_PI5P = -0.2383381, f4_0a_PI345P3 = -0.001668266, f4_0a_PI45P2_PI = -0.01096355,
		gamma4_0c = 5.096813e+12, f4_0c_PI4P = 0.93781095, f4_0c_PI3P = -0.0006218905, f4_0c_PI5P = -0.0006218905, f4_0c_PI35P2 = -0.0003109453,


		gamma0_5 = 4.702825e+11, f0_5_PI = 0.5465332, f0_5_PI3P = -0.000377889,  
		gamma5_0a = 4.557615e+13, f5_0a_PI3P = -0.0004385422, f5_0a_PI4P = -0.04385422, f5_0a_PI5P = 0.9995614578, f5_0a_PI35P2 = -0.0002192711, f5_0a_PI45P2_PI4P = -0.04385422, f5_0a_PI45P2_PI5P = -0.2383381, f5_0a_PI345P3 = -0.001668266, f5_0a_PI45P2_PI = -0.01096355,
		gamma5_0c = 5.096813e+12, f5_0c_PI5P = 0.9993781095, f5_0c_PI3P = -0.0006218905, f5_0c_PI4P = -0.06218905, f5_0c_PI35P2 = -0.0003109453,

		gamma3_35 = 1.618325e+16, f3_35_PI3P = 0.999622111, f3_35_PI = -0.4534668,
		gamma35_3a = 4.557615e+13, f35_3a_PI3P = -0.0004385422, f35_3a_PI4P = -0.04385422, f35_3a_PI5P = -0.0004385422, f35_3a_PI35P2 = 0.9997807289, f35_3a_PI45P2_PI4P = -0.04385422, f35_3a_PI45P2_PI5P = -0.2383381, f35_3a_PI345P3 = -0.001668266, f35_3a_PI45P2_PI = -0.01096355,
		gamma35_3c = 5.096813e+12, f35_3c_PI35P2 = 0.9996890547, f35_3c_PI3P = -0.0006218905, f35_3c_PI4P = -0.06218905, f35_3c_PI5P = -0.0006218905,
		gamma35_3d = 3.060265e+12, f35_3d_PI35P2 = 9.999336e-01,
		gamma35_3e = 5.137691e+12, f35_3e_PI35P2 = 9.996894e-01, f35_3e_PI45P2 = -6.211936e-02, f35_3e_PI345P3 = -2.363093e-03,

		gamma4_45 = 8.493641e+15, f4_45_PI4P = 4.596175e-02, f4_45_PA = 0.2, f4_45_PI45P2 = -0.05,
		gamma45_4a = 4.557615e+13, f45_4a_PI3P = -0.0004385422, f45_4a_PI4P = -0.04385422, f45_4a_PI5P = -0.0004385422, f45_4a_PI35P2 = -0.0002192711, f45_4a_PI45P2_PI4P = 0.95614578, f45_4a_PI45P2_PI5P = -0.2383381, f45_4a_PI345P3 = -0.001668266, f45_4a_PI45P2_PI = -0.01096355,
	 	gamma45_4c = 5.133984e+12, f45_4c_PI45P2 = 9.378613e-01, f45_4c_PI345P3 = -2.363827e-03,
		gamma45_4d = 5.137691e+12, f45_4d_PI45P2 = 9.378806e-01, f45_4d_PI35P2 = -3.105968e-04, f45_4d_PI345P3 = -2.363093e-03,

		gamma5_45 = 2.946989e+13, f5_45_PI5P = 8.784401e-01,
		gamma45_5a = 3.237862e+12, f45_5a_PI3P = -0.0004385422, f45_5a_PI4P = -0.04385422, f45_5a_PI5P = -0.0004385422, f45_5a_PI35P2 = -0.0002192711, f45_5a_PI45P2_PI4P = -0.04385422, f45_5a_PI45P2_PI5P = 0.7616619, f45_5a_PI345P3 = -0.001668266, f45_5a_PI45P2_PI = -0.01096355,
		gamma45_5c = 4.133048e+12, f45_5c_PI45P2 = 6.487218e-01,   	
   	
		gamma45_345 = 1.888500e+14, f45_345_PI45P2 = 3.063327e-01,
		gamma345_45 = 2.904971e+15, f345_45_PI345P3 = 9.555746e-01, f345_45_PI34P2 = -1.101822e-04,
		gamma35_5 = 6.059076e+14, f35_5_PI35P2 = 9.996682e-01, f35_5_PI3P = -6.635700e-04,
	
		gamma34_3 = 1.333061e+12, f34_3_PI34P2 = 9.981983e-01,

		gamma345_34a = 2.476815e+12, f345_34a_PI3P = -0.0004385422, f345_34a_PI4P = -0.04385422, f345_34a_PI5P = -0.0004385422, f345_34a_PI35P2 = -0.0002192711, f345_34a_PI45P2_PI4P = -0.04385422, f345_34a_PI45P2_PI5P = -0.2383381, f345_34a_PI345P3 = 0.998331734, f345_34a_PI45P2_PI = -0.01096355,
		gamma345_34c = 2.790040e+11, f345_34c_PI345P3 = 9.976362e-01, f345_34c_PI45P2 = -6.213866e-02,
		gamma345_34d = 2.792054e+11, f345_34d_PI345P3 = 9.976369e-01, f345_34d_PI35P2 = -3.105968e-04, f345_34d_PI45P2 = -6.211936e-02,
		gamma345_34e = 1.678636e+11, f345_34e_PI345P3 = 9.981984e-01,

		gammai_ = 0.045,
		gamma0_45 = 2.666243e+14, f0_45_PI = 2.864618e-01, f0_45_PA = 0.2, f0_45_PI45P2 = -0.05,
		gamma4_34a = 0*5.229823e+12, f4_34a_PI4P = 5.009150e-01,
		gamma4_34b = 5.638138e+14, f4_34b_PI4P = 0.92969752, f4_34b_PI = -0.8238571 ,
		gamma34_4 = 5.043352e+11, f34_4_PI34P2 = 9.998898e-01, f34_4_PI345P3 = -4.442544e-02,
		gamma45_0 = 1.139404e+12, f45_0_PI3P = 0.9995614578, f45_0_PI4P = -0.04385422, f45_0_PI5P = -0.0004385422, f45_0_PI35P2 = -0.0002192711, f45_0_PI45P2_PI4P = -0.04385422, f45_0_PI45P2_PI5P = -0.2383381, f45_0_PI345P3 = -0.001668266, f45_0_PI45P2_PI = -0.01096355
		)


				
parameters=c(parameters,
		gammaPI3KIc_PI3KIa = mp2_PI3KI * 0.005385984, fPI3KIc_PI3KIa_pi_3KI_c = 1, fPI3KIc_PI3KIa_PI345P3 = .7, 
		gammaPI3KIa_PI3KIc = mp2_PI3KI, fPI3KIa_PI3KIc_pi_3KI_a = 1,
 
		gammaPTENc_PTENa = mp3_PTEN * 4.166667e-05, fPTENc_PTENa_PTEN_c = 1, fPTENc_PTENa_PI45P2 = 1, 
		gammaPTENa_PTENc = mp3_PTEN, fPTENa_PTENc_PTEN_a =1, fPTENa_PTENc_PI345P3 = 1,

		gamma45_DAG = 1.878461e+19, f45_DAG_PI45P2 = 9.770418e-01, f45_DAG_PI4P = -1.260449e-02, f45_DAG_PS = -9.642438e-01,
		gammaPA_DAG = 5.538472e+11, fPA_DAG_PA = 9.652797e-01,						
		gammaDAG_PA = 1.652729e+13, fDAG_PA_DAG = 9.475506e-01,
		gamma_DAG = 0.00001,
		gammaDAG_ = .1, fDAG__DAG = 1,
		gamma_PA = 4,
		gammaPA_ = 0.1, fPA__PA = 1,
		gammaPC_PA = 14.35507, f_PC_PA_PI45P2 = .3,
		gammaIP3_ = 2, fIP3__IP3 = 1,
		gamma4_DAG = 2.546895e+18, f4_DAG_PI4P = 9.873955e-01, f4_DAG_PI45P2 = -2.295818e-02, f4_DAG_PS = -9.642438e-01,

		SPLUNC1 = 7890,

		gamma__ENaC = log(2)/40,
		gamma_ENaC_SPLUNC1 = 3*log(2)/8416000,
		gamma_ENaC_ = log(2)/3200,
		
		gamma_ASL_ENaC = 0.0003922396,
		gamma_CFTR_ASL = 0.0177052536
		)   
parameters=c(
		parameters,
		gamma__ASL = (105)*as.numeric(parameters[names(parameters)=='gamma_ASL_ENaC'])+(4/3)*as.numeric(parameters[names(parameters)=='gamma_CFTR_ASL']),
		gamma_ASL_ = (25/4)*as.numeric(parameters[names(parameters)=='gamma_ASL_ENaC'])+(1/3)*as.numeric(parameters[names(parameters)=='gamma_CFTR_ASL'])	
		)   
#Pars=as.data.frame(t(parameters)) #create a data.frame p with the same information in vector parameters
#attach(Pars) # put all colunms of p as variables in memory
# cbind(1:length(parameters),parameters)

equations=function(t,state,parameters) 	# function containing the diferential equations
	{with(as.list(c(state,parameters)),
		{
		#rate of change (velocities of the variables concentration changes), the actual diferencial equations
		# insert the diferentical equations (ex: dX=const*X^rate or dX=5*X^2)

		V_0 = gamma_0
		V_4 = gamma_4
		V_3 = gamma_3

		V0_3a = gamma0_3a*PI^f0_3a_PI * pi_3KI_a * 123.245 * 1.66*10^-18
		V0_3b = gamma0_3b * PI^f0_3b_PI * PI4P^f0_3b_PI4P * pi_3KII * 180.983 * 1.66*10^-18
		V0_3c = gamma0_3c*PI^f0_3c_PI*pi_3KIII  * 101.549 * 1.66*10^-18
		V3_0a = gamma3_0a * PI3P^f3_0a_PI3P * PI4P^f3_0a_PI4P * PI5P^f3_0a_PI5P * PI35P2^f3_0a_PI35P2 * PI45P2^f3_0a_PI45P2_PI4P * PI45P2^f3_0a_PI45P2_PI5P * PI345P3^f3_0a_PI345P3 * PI45P2^f3_0a_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18
		V3_0c = gamma3_0c * PI3P^f3_0c_PI3P * PI4P^f3_0c_PI4P * PI5P^f3_0c_PI5P * PI35P2^f3_0c_PI35P2 * SAC1 * 66.967 * 1.66*10^-18
		V3_0d = gamma3_0d * PI3P^f3_0d_PI3P * PI35P2^f3_0d_PI35P2 * MTMR1_6_14 * 69.932 * 1.66*10^-18
		V3_0e = gamma3_0e*PI3P^f3_0e_PI3P * MTMR78 * 75.833 * 1.66*10^-18

		V0_4 = gamma0_4*PI^f0_4_PI* pi_4K * 109.711 * 1.66*10^-18
		V4_0a = gamma4_0a * PI3P^f4_0a_PI3P * PI4P^f4_0a_PI4P * PI5P^f4_0a_PI5P * PI35P2^f4_0a_PI35P2 * PI45P2^f4_0a_PI45P2_PI4P * PI45P2^f4_0a_PI45P2_PI5P * PI345P3^f4_0a_PI345P3 * PI45P2^f4_0a_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18
		V4_0c = gamma4_0c * PI4P^f4_0c_PI4P * PI3P^f4_0c_PI3P * PI5P^f4_0c_PI5P * PI35P2^f4_0c_PI35P2 * SAC1 * 66.967 * 1.66*10^-18

		V0_5 = gamma0_5 * PI^f0_5_PI * PI3P^f0_5_PI3P * pi_Kfyve * 237.136 * 1.66*10^-18
		V5_0a = gamma5_0a * PI3P^f5_0a_PI3P * PI4P^f5_0a_PI4P * PI5P^f5_0a_PI5P * PI35P2^f5_0a_PI35P2 * PI45P2^f5_0a_PI45P2_PI4P * PI45P2^f5_0a_PI45P2_PI5P * PI345P3^f5_0a_PI345P3 * PI45P2^f5_0a_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18
		V5_0c = gamma5_0c * PI5P^f5_0c_PI5P * PI3P^f5_0c_PI3P * PI4P^f5_0c_PI4P * PI35P2^f5_0c_PI35P2 * SAC1 * 66.967 * 1.66*10^-18

		V3_35 = gamma3_35 * PI3P^f3_35_PI3P * PI^f3_35_PI * pi_Kfyve * 237.136 * 1.66*10^-18
		V35_3a = gamma35_3a * PI3P^f35_3a_PI3P * PI4P^f35_3a_PI4P * PI5P^f35_3a_PI5P * PI35P2^f35_3a_PI35P2 * PI45P2^f35_3a_PI45P2_PI4P * PI45P2^f35_3a_PI45P2_PI5P * PI345P3^f35_3a_PI345P3 * PI45P2^f35_3a_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18
		V35_3c = gamma35_3c * PI35P2^f35_3c_PI35P2 * PI3P^f35_3c_PI3P * PI4P^f35_3c_PI4P * PI5P^f35_3c_PI5P * SAC1 * 66.967 * 1.66*10^-18
		V35_3d = gamma35_3d*PI35P2^f35_3d_PI35P2* SAC3 * 103.635 * 1.66*10^-18
		V35_3e = gamma35_3e * PI35P2^f35_3e_PI35P2 * PI45P2^f35_3e_PI45P2 * PI345P3^f35_3e_PI345P3 * INPP5E * 70.205 * 1.66*10^-18

		V4_45 = gamma4_45 * PI4P^f4_45_PI4P * pip_5KI * 65.643 *1.66*10^-18 * HS(PA) * PI45P2^f4_45_PI45P2
		V45_4a = gamma45_4a * PI3P^f45_4a_PI3P * PI4P^f45_4a_PI4P * PI5P^f45_4a_PI5P * PI35P2^f45_4a_PI35P2 * PI45P2^f45_4a_PI45P2_PI4P * PI45P2^f45_4a_PI45P2_PI5P * PI345P3^f45_4a_PI345P3 * PI45P2^f45_4a_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18
		V45_4c = gamma45_4c * PI45P2^f45_4c_PI45P2 * PI345P3^f45_4c_PI345P3 * (ORCL1 * 104.205 * 1.66*10^-18 + INPP5BJ * 96.755 * 1.66*10^-18 + SKIP * 51.090 * 1.66*10^-18 + SAC2 * 128.407* 1.66*10^-18)
 		V45_4d = gamma45_4d * PI45P2^f45_4d_PI45P2 * PI35P2^f45_4d_PI35P2 * PI345P3^f45_4d_PI345P3 * INPP5E * 70.205 * 1.66*10^-18

		V5_45 = gamma5_45*PI5P^f5_45_PI5P*pip_5KII * 46.968 *1.66*10^-18
		V45_5a = gamma45_5a * PI3P^f45_5a_PI3P * PI4P^f45_5a_PI4P * PI5P^f45_5a_PI5P * PI35P2^f45_5a_PI35P2 * PI45P2^f45_5a_PI45P2_PI4P * PI45P2^f45_5a_PI45P2_PI5P * PI345P3^f45_5a_PI345P3 * PI45P2^f45_5a_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18
		V45_5c = gamma45_5c*PI45P2^f45_5c_PI45P2 * TMEM55 * 28.470 * 1.66*10^-18

		V45_345 = gamma45_345*PI45P2^f45_345_PI45P2 * pi_3KI_a * 123.245 * 1.66*10^-18
		V345_45 = gamma345_45 * PI345P3^f345_45_PI345P3 * PI34P2^f345_45_PI34P2 * PTEN_a * 47.166 * 1.66*10^-18

		V35_5 = gamma35_5 * PI35P2^f35_5_PI35P2 * PI3P^f35_5_PI3P * MTMR1_6_14 * 69.932 * 1.66*10^-18

		V34_3 = gamma34_3*PI34P2^f34_3_PI34P2 * INPP4 * 107.347 * 1.66*10^-18

		V345_34a = gamma345_34a * PI3P^f345_34a_PI3P * PI4P^f345_34a_PI4P * PI5P^f345_34a_PI5P * PI35P2^f345_34a_PI35P2 * PI45P2^f345_34a_PI45P2_PI4P * PI45P2^f345_34a_PI45P2_PI5P * PI345P3^f345_34a_PI345P3 * PI45P2^f345_34a_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18
		V345_34c = gamma345_34c * PI345P3^f345_34c_PI345P3 * PI45P2^f345_34c_PI45P2 * (ORCL1 * 104.205 * 1.66*10^-18 + INPP5BJ * 96.755 * 1.66*10^-18 + SKIP  * 51.090 * 1.66*10^-18 + SAC2 * 128.407* 1.66*10^-18)
		V345_34d = gamma345_34d * PI345P3^f345_34d_PI345P3 * PI35P2^f345_34d_PI35P2 * PI45P2^f345_34d_PI45P2 * INPP5E * 70.205 * 1.66*10^-18
		V345_34e = gamma345_34e*PI345P3^f345_34e_PI345P3* SHIP2 * 138.599 * 1.66*10^-18

		V45_  = gammai_*PI45P2
		V0_   = gammai_*PI
		V4_   = gammai_*PI4P
		V345_ = gammai_*PI345P3
		V3_   = gammai_*PI3P
		V35_  = gammai_*PI35P2
		V5_   = gammai_*PI5P
		V34_  = gammai_*PI34P2

		V0_45 = gamma0_45 * PI^f0_45_PI * pi_4K_pip_5KI * 248.608 * 1.66*10^-18  * HS(PA) * PI45P2^f0_45_PI45P2

		V4_34a = gamma4_34a*PI4P^f4_34a_PI4P* pi_3KI_a * 123.245 * 1.66*10^-18
		V4_34b = gamma4_34b * PI4P^f4_34b_PI4P * PI^f4_34b_PI * pi_3KII * 180.983 * 1.66*10^-18 
		V34_4 = gamma34_4 * PI34P2^f34_4_PI34P2 * PI345P3^f34_4_PI345P3 * PTEN_a * 47.166 * 1.66*10^-18

		V45_0 = gamma45_0 * PI3P^f45_0_PI3P * PI4P^f45_0_PI4P * PI5P^f45_0_PI5P * PI35P2^f45_0_PI35P2 * PI45P2^f45_0_PI45P2_PI4P * PI45P2^f45_0_PI45P2_PI5P * PI345P3^f45_0_PI345P3 * PI45P2^f45_0_PI45P2_PI * SYNJ * 185.272 *1.66*10^-18

		pi_3KI_c = pi_3KI - pi_3KI_a
		VPI3KIc_PI3KIa = gammaPI3KIc_PI3KIa * (pi_3KI - pi_3KI_a) ^ fPI3KIc_PI3KIa_pi_3KI_c * PI345P3 ^ fPI3KIc_PI3KIa_PI345P3
		VPI3KIa_PI3KIc = gammaPI3KIa_PI3KIc * pi_3KI_a

		PTEN_c = PTEN - PTEN_a
		VPTENc_PTENa = gammaPTENc_PTENa * (PTEN - PTEN_a) ^ fPTENc_PTENa_PTEN_c * PI45P2 ^ fPTENc_PTENa_PI45P2
		VPTENa_PTENc = gammaPTENa_PTENc * PTEN_a ^ fPTENa_PTENc_PTEN_a 

		V45_DAG = gamma45_DAG * PI45P2^(f45_DAG_PI45P2) * PI4P^(f45_DAG_PI4P) * PS^(f45_DAG_PS) * PLC * 88.4215 *1.66*10^-18
		VPA_DAG = gammaPA_DAG * PA^fPA_DAG_PA * LPP * 32.574 * 1.66*10^-18
		VDAG_PA = gammaDAG_PA * DAG^fDAG_PA_DAG * DGK * 116.997 *1.66*10^-18
		V_DAG = gamma_DAG 
		VDAG_ = gammaDAG_ * DAG^fDAG__DAG 
		V_PA = gamma_PA 
		VPA_ = gammaPA_ * PA^fPA__PA
		VPC_PA = gammaPC_PA * PI45P2^f_PC_PA_PI45P2 
		VIP3_ = gammaIP3_ * IP3^fIP3__IP3

		V_PI4P_DAG = gamma4_DAG * PI4P^(f4_DAG_PI4P) * PI45P2^(f4_DAG_PI45P2) * PS^(f4_DAG_PS) * PLC * 88.4215 *1.66*10^-18

		dPI = V3_0a + V3_0c + V3_0d + V3_0e + V4_0a + V4_0c + V5_0a  + V5_0c + V45_0 - V0_3a - V0_3b - V0_3c - V0_4 - V0_5 - V0_ - V0_45 + 19.05246 * VPA_
		dPI3P = V_3 + V0_3a + V0_3b + V0_3c + V35_3a + V35_3c + V35_3d + V35_3e + V34_3 - V3_0a - V3_0c - V3_0d - V3_0e - V3_35 - V3_
		dPI4P = V_4 + V0_4 + V45_4a + V45_4c + V45_4d + V34_4 - V4_0a  - V4_0c - V4_45 - V4_ - V4_34a - V4_34b - V_PI4P_DAG
		dPI5P = V0_5 + V35_5 + V45_5a + V45_5c - V5_0a - V5_0c - V5_45 - V5_
		dPI35P2 = V3_35 - V35_3a - V35_3c - V35_3d - V35_3e - V35_5 - V35_
		dPI45P2 = V4_45 + V5_45 + V345_45 + V0_45 - V45_4a - V45_4c - V45_4d - V45_5a  - V45_5c - V45_345 - V45_ - V45_0 - V45_DAG
		dPI34P2 = V345_34a + V345_34c + V345_34d + V345_34e + V4_34a + V4_34b - V34_3 - V34_ - V34_4
		dPI345P3 = V45_345 - V345_45 - V345_34a - V345_34c - V345_34d - V345_34e - V345_		
	
		dpi_3KI_a = (VPI3KIc_PI3KIa - VPI3KIa_PI3KIc) 
		dPTEN_a = VPTENc_PTENa - VPTENa_PTENc 

		dDAG = V_DAG + VPA_DAG  + V_PI4P_DAG + V45_DAG - VDAG_PA - VDAG_
		dIP3 = V45_DAG - VIP3_
		dPA = V_PA + VDAG_PA + VPC_PA - VPA_DAG - VPA_ 

		######################## ENaC SPLUNC1 model #######################

		V__ENaC = gamma__ENaC
		V_ENaC_SPLUNC1 = gamma_ENaC_SPLUNC1 * ENaC * (SPLUNC1/ASL) 
		V_ENaC_ = gamma_ENaC_ * ENaC

		V__ASL = gamma__ASL 
		V_CFTR_ASL = gamma_CFTR_ASL 
		V_ASL_ENaC = gamma_ASL_ENaC * ASL * (ENaC * ENaC_op(PI45P2))
		V_ASL_ = gamma_ASL_ * ASL

		dENaC = V__ENaC - V_ENaC_SPLUNC1 - V_ENaC_ 
		dASL = V__ASL + V_CFTR_ASL - V_ASL_ENaC - V_ASL_



		# return the rate of change 
		# insert, in the parentecis, the variables that you want to store the values (ex:dX)
		list(dy=c(
				dPI, dPI3P, dPI4P, dPI5P,dPI35P2,dPI45P2,dPI34P2,dPI345P3,
				dpi_3KI_a,dPTEN_a,
				dDAG,dIP3,dPA,
				dENaC,dASL
				),
				count=c(
					V_0,V0_3a,V0_3b,V0_3c,V3_0a,V3_0c,V3_0d,V3_0e,V0_4,V4_0a,V4_0c,V0_5,V5_0a,V5_0c,V3_35,V35_3a,V35_3c,V35_3d,V35_3e,V4_45,V45_4a,V45_4c,V45_4d,V5_45,V45_5a,V45_5c,V45_345,V345_45,V_4,V35_5,V_3,V34_3,V345_34a,V345_34c,V345_34d,V345_34e,
					V45_,V0_,V4_,V345_,V3_,V35_,V5_,V34_,
					V0_45,V4_34a,V4_34b,V34_4,V45_0,
					VPI3KIc_PI3KIa,VPI3KIa_PI3KIc,VPTENc_PTENa,VPTENa_PTENc,
					V45_DAG,VPA_DAG,VDAG_PA,V_DAG,VDAG_,V_PA,VPA_,VPC_PA,VIP3_,
					V_PI4P_DAG,
					V__ENaC,V_ENaC_SPLUNC1,V_ENaC_,
					V__ASL,V_CFTR_ASL,V_ASL_ENaC,V_ASL_	
					))
		})
	}

perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA'))	# initiation of perturbations matrix, if unaltered it will do no perturbations






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
		if (max(t)<max(perturbations[1:(nrow(perturbations)-1),1])) {cat("Time error: Last perturbation later than tmax.Please alter tmax or perturbation time.", "\n")}
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




# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
#ptime=c(0,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,8000,max(t));
# On the second column put the variables to be altered. (ex: 'X')
#pvar=c('novar','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','PLC','DGK','DGK','novar')
# On the third the new value for the variable.
#pval=c(NA,
#	parameters[names(parameters)=='PLC']*0.01,
#	parameters[names(parameters)=='PLC'],
#	parameters[names(parameters)=='PLC']*1.5,
#	parameters[names(parameters)=='PLC'],
#	parameters[names(parameters)=='PLC']*2,
#	parameters[names(parameters)=='PLC'],
#	parameters[names(parameters)=='PLC']*3,
#	parameters[names(parameters)=='PLC'],
#	parameters[names(parameters)=='PLC']*4,
#	parameters[names(parameters)=='PLC'],
#	parameters[names(parameters)=='PLC']*5,
#	parameters[names(parameters)=='PLC'],
#	parameters[names(parameters)=='DGK']*0.25,
#	parameters[names(parameters)=='DGK'],
#	NA)
#perturbations=data.frame(time=ptime,var=pvar,val=pval)
#perturbations





out = Cruncher(state,t,equations,parameters,perturbations)
# edit(out)
# View(out)

str(out)
varnames=c('Time (s)','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3',
		'pi_3KI_a','PTEN_a','DAG','IP3','PA',
		'ENaC','ASL'
		) 
varcolors=c('gray','tan','orange3','orange','orange3','gold3','gold','gold3','yellow','red','deepskyblue3','salmon3','yellow2','salmon',
		'green','cyan'
		)
#  cbind(1:length(varnames),varnames,varcolors)

# naming fluxes in out
flux_names=c(
		'V_0','V0_3a','V0_3b','V0_3c','V3_0a','V3_0c','V3_0d','V3_0e','V0_4','V4_0a','V4_0c','V0_5','V5_0a','V5_0c',
		'V3_35','V35_3a','V35_3c','V35_3d','V35_3e','V4_45','V45_4a','V45_4c','V45_4d','V5_45','V45_5a','V45_5c',
		'V45_345','V345_45','V_4','V35_5','V_3','V34_3','V345_34a','V345_34c','V345_34d','V345_34e','V45_','V0_','V4_',
		'V345_','V3_','V35_','V5_','V34_','V0_45','V4_34a','V4_34b','V34_4','V45_0',
		'VPI3KIc_PI3KIa','VPI3KIa_PI3KIc',
		'VPTENc_PTENa','VPTENa_PTENc',
		'V45_DAG','VPA_DAG','VDAG_PA','V_DAG','VDAG_','V_PA','VPA_','VPC_PA','VIP3_',
		'V_PI4P_DAG',
		'V__ENaC','V_ENaC_SPLUNC1','V_ENaC_',
		'V__ASL','V_CFTR_ASL','V_ASL_ENaC','V_ASL_'
		)
colnames(out)[(length(state)+2):dim(out)[2]]=flux_names

(tail(out,3))
(tail(out,3))[,1:16]



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

windows()
par(mfrow=c(2,1))
plot(out[,1],out[,10],type='l',col=1,ylim=c(0,max(out[,10:11])),main='PIK3 PTEN',xlab='Time (m)',ylab='Molecule count')
points(out[,1],out[,10],type='l',col=varcolors[10])
text(tmax,out[tmax-1,10],varnames[10],col=varcolors[10])
points(out[,1],out[,11],type='l',col=varcolors[11])
text(tmax,out[tmax-1,11],varnames[11],col=varcolors[11])

plot(out[,1],out[,12],type='l',col='blue',ylim=c(0,max(out[,12:14])),main='Time course',xlab='Time (m)',ylab='Molecule count')
points(out[,1],out[,12],type='l',col=varcolors[12])
text(tmax,out[tmax-1,12],varnames[12],col=varcolors[12])
points(out[,1],out[,13],type='l',col=varcolors[13])
text(tmax,out[tmax-1,13],varnames[13],col=varcolors[13])
points(out[,1],out[,14],type='l',col=varcolors[14])
text(tmax,out[tmax-1,14],varnames[14],col=varcolors[14])
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
         
# list (primitive)
(head(out,3))[,1:16]
(tail(out,3))[,1:16]

sum(out[nrow(out),2:9])

# edit(out)



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
par(mfrow=c(3,3))
for (i in c(2:9)) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

windows()
par(mfrow=c(3,3))
for (i in c(12,13,14,15,16)) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# par(mfrow=c(2,1))
# i=7
# plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i],lwd=5)
# text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
# i=14
# plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i],lwd=5)
# text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
# par(mfrow=c(1,1))
# out_norm[7900,7]

(tail(out,3))[,1:16]









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
pvar=c('novar','DGK','DGK','novar')
# On the third the new value for the variable.
pval=c(NA,
	parameters[names(parameters)=='DGK']*0.25,
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
	parameters[names(parameters)=='DGK']*0.25,
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

graphics.off() 	#close plot windows

# tiff("P2F2_Fig_11c.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

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


