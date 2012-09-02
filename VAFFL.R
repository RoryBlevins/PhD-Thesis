# function



gillespied=function (N, T=100, dt=1, ...)
{
        tt=0
        n=T%/%dt
        x=N$M
        S=t(N$Post-N$Pre)
        u=nrow(S)
        v=ncol(S)
        xmat=matrix(0,ncol=u,nrow=n)
	i=1
	target=0
	repeat {
            h=N$h(x, th, ...)
		h0=sum(h)
		if (h0<1e-10)
			tt=1e99
		else
	                tt=tt+rexp(1,h0)
		while (tt>=target) {
			xmat[i,]=x
			i=i+1
			target=target+dt
			if (i>n)
			        return(ts(xmat,start=0,deltat=dt))
		}
                j=sample(v,1,prob=h)
                x=x+S[,j]
        }
}





# test of miRNA-containing FFL

FFL=list()
#starting values
FFL$M=c(
100, 	# TF mRNA
100, 	# TF protein
100,	# miRNA
100,	# Output mRNA
100,	# Output protein
0	# miRNA-protein complex
)

#Reaction matrices
FFL$Pre=matrix(c(
0,0,0,0,0,0,	# TF transcription
1,0,0,0,0,0,	# TF mRNA degradation
1,0,0,0,0,0,	# TF translation
0,1,0,0,0,0,	# TF protein degradation
0,0,0,0,0,0,	# miRNA transcription
0,0,1,0,0,0,	# miRNA degradation
0,0,1,1,0,0,	# miRNA-output mRNA binding
0,0,0,0,0,1,	# miRNA-output mRNA dissociation
0,0,0,0,0,1,	# miRNA-output mRNA degradation of mRNA
0,0,0,0,0,0,	# Output transcription
0,0,0,1,0,0,	# Output mRNA degradation
0,0,0,1,0,0,	# Output translation
0,0,0,0,1,0	  	# Output protein degradation
),ncol=6,byrow=TRUE)


FFL$Post=matrix(c(
1,0,0,0,0,0,	# TF transcription
0,0,0,0,0,0,	# TF mRNA degradation
1,1,0,0,0,0,	# TF translation
0,0,0,0,0,0,	# TF protein degradation
0,0,1,0,0,0,	# miRNA transcription
0,0,0,0,0,0,	# miRNA degradation
0,0,0,0,0,1,	# miRNA-output mRNA binding
0,0,1,1,0,0,	# miRNA-output mRNA dissociation
0,0,1,0,0,0,	# miRNA-output mRNA degradation of mRNA
0,0,0,1,0,0,	# Output transcription
0,0,0,0,0,0,	# Output mRNA degradation
0,0,0,1,1,0,	# Output translation
0,0,0,0,0,0	  	# Output protein degradation
),ncol=6,byrow=TRUE)





#rate constant calculation
FFL$h=function(x,th)
{
 return(c(
th[1], 								#TF transcription
th[2]*x[1],								#TF mRNA degradation
th[3]*x[1],								#TF translation
th[4]*x[2],								#TF protein degradation								
(th[5]*x[2]^th[7])/(th[6]^th[7]+x[2]^th[7]),		#miRNA transcription
th[8]*x[3],								#miRNA degradation
th[9]*x[3]*x[4],							# miRNA-output mRNA binding
th[10]*x[6]	,							# miRNA-output mRNA dissociation
th[11]*x[6]	,							# miRNA-output mRNA degradation of mRNA
(th[12]*x[2]^th[7])/(th[6]^th[7]+x[2]^th[7]),	#output transcription
th[2]*x[4],							#output mRNA degradation
th[13]*x[4],							#output translation
th[14]*x[5]								#ouput degradation
))}



# test of miRNA-independent circuit

TF=list()
#starting values
TF$M=c(
100, 	# TF mRNA
100, 	# TF protein
100,	# miRNA
100,	# Output mRNA
100,	# Output protein
0	# miRNA-protein complex
)

#Reaction matrices
TF$Pre=matrix(c(
0,0,0,0,0,0,	# TF transcription
1,0,0,0,0,0,	# TF mRNA degradation
1,0,0,0,0,0,	# TF translation
0,1,0,0,0,0,	# TF protein degradation
0,0,0,0,0,0,	# Output transcription
0,0,0,1,0,0,	# Output mRNA degradation
0,0,0,1,0,0,	# Output translation
0,0,0,0,1,0	  	# Output protein degradation
),ncol=6,byrow=TRUE)


TF$Post=matrix(c(
1,0,0,0,0,0,	# TF transcription
0,0,0,0,0,0,	# TF mRNA degradation
1,1,0,0,0,0,	# TF translation
0,0,0,0,0,0,	# TF protein degradation
0,0,0,1,0,0,	# Output transcription
0,0,0,0,0,0,	# Output mRNA degradation
0,0,0,1,1,0,	# Output translation
0,0,0,0,0,0	  	# Output protein degradation
),ncol=6,byrow=TRUE)




#rate constant calculation
TF$h=function(x,th)
{
 return(c(
th[1], 								#TF transcription
th[2]*x[1],								#TF mRNA degradation
th[3]*x[1],								#TF translation
th[4]*x[2],								#TF protein degradation								
(th[12]*x[2]^th[7])/(th[6]^th[7]+x[2]^th[7]),	#output transcription
th[2]*x[4],							#output mRNA degradation
th[13]*x[4],							#output translation
th[14]*x[5]								#ouput degradation
))}


#Timecourse data

T<- 5000
dt<- 10

measuretime <- 4000/dt

#number of simulations
n <- 10000


	


	#rate constants
	th=c(
	0.06, 	# transcription rate
	0.006, 	# TF mRNA degradation rate
	0.04, 	# TF translation rate
	0.002,	# TF protein degradation rate
	0.5, 		# base miRNA transcription rate
	200, 		# miRNA dissociation coefficient
	2,		# hill coefficient for miRNA 
	0.006,	# miRNA degradation rate
	0.0001,	# miRNA complex formation rate
	0.0001,	# miRNA-mRNA dissociation rate
	0.02, 	# miRNA-mRNA degradation rate
	0.8,		# base output transcription rate
	#200, 		# mRNA dissociation coefficient
	#2,		# hill coefficient for transcription
	#0.006, 	# output mRNA degradation rate
	0.04, 	# base translation rate
	0.002		# output protein degradation
	)	

	#Run simulations


		FFLresults <- 0
		TFresults <- 0


		for (ii in (1:n)){
			FFLOutput <-compgillespied(FFL,T,dt)
			TFOutput <-compgillespied(TF,T,dt)

			if (ii%%10==0){
				print(ii)
				}

			FFLresults <- rbind(FFLresults,FFLOutput[measuretime,])
			TFresults <- rbind(TFresults, TFOutput[measuretime,])
			}

		TFresults <- TFresults[-1,]

		FFLresults <- FFLresults[-1,]

		TFMean <- mean(TFresults[,5])

		FFLMean <- mean(FFLresults[,5])

		TFSD <- sd(TFresults[,5])

		FFLSD <- sd(FFLresults[,5])

		RelativeMean <- TFMean/FFLMean

		FFLCV <- FFLSD/FFLMean

		TFCV <- TFSD/TFMean

		RelativeCV <- TFCV/FFLCV

		RelativeMean

		RelativeCV



write.csv(FFLresults, "VAFFLresults.csv",)
write.csv(TFresults, "VATFresults.csv",)

# eof

