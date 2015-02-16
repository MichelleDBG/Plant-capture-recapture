# Plant-capture-recapture
R and MARK methods of capture-recapture


################################################
# Rcapture Shefferson et al. 2001
# http://cran.r-project.org/web/packages/Rcapture/Rcapture.pdf

library(Rcapture)

library(plyr)
library(plotrix)

#Open population model with Cormack-Jolly-Seber


library(plyr)
library(plotrix)
library(ggplot2)
library(reshape2)



#	"Q:/Research/All_Projects_by_Species/Astragalus SPECIES/
#	Astragalus_microcymbus/
#	AsMi_metadata/Scouting/20140421_Data/April 21 2014_May 13_.csv"
#	"Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/R_Analysis/R_code_multiyearvisits/"

setwd(path.expand("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/R_Analysis/R_code_multiyearvisits/"))

## Delete non-data rows in the .csv
asmi14 <- #read.csv("April 21 2014_May 13_June 2.csv", 
	    #read.csv("April 21 2014_May 13_June 2_June 25_July 15.csv",
	    #read.csv("April 21 2014_May 13_June 2_June 25_July 15_August5.csv",
	read.csv("2014_AsMi_3Weekly.csv",
	header = F, as.is = TRUE,
	blank.lines.skip = T ,	
	skip = 1,
	colClasses = c('numeric', 'numeric', 'numeric', 'factor', 'character',
				'numeric','character',						#2013 Ln and comments 1-7
				'numeric','factor','numeric','factor','character',	# 8-12
				'character', 							# april seedlings 13
				'numeric','factor','numeric','factor','character',	#May 14-18
				'numeric','factor','numeric','factor','character',	#June 19-23
				'numeric','factor','numeric','factor','character',	#June2 24-28
				'numeric','factor','numeric','factor','character',	#Jul 29-33
				'numeric','factor','numeric','factor','character',	#Aug 34-38
				'numeric','factor','numeric','factor','character',	#Sept	39-43
				'numeric','factor','numeric','factor','character',	#Oct	44-48
				'numeric','factor','numeric','factor','character',	#Oct2	49-53
				'numeric'))	#tag number 54

head(asmi14)
str(asmi14)
colnames(asmi14)<- c('Tag','Plot','Dist','Dir','Location','Ln13','Comments13',	# 1-7
				'LnApr14','FlApr14','FrApr14','BrApr14','CommentsApr14',	# 8-12
				'AprSeedlings',								# 13
				'LnMay14','FlMay14','FrMay14','BrMay14','CommentsMay14',	# 14-18
				'LnJun14','FlJun14','FrJun14','BrJun14','CommentsJun14',
				'LnJun214','FlJun214','FrJun214','BrJun214','CommentsJun214',
				'LnJul14','FlJul14','FrJul14','BrJul14','CommentsJul14',
				'LnAug14','FlAug14','FrAug14','BrAug14','CommentsAug14',
				'LnSep14','FlSep14','FrSep14','BrSep14','CommentsSep14',
				'LnOct14','FlOct14','FrOct14','BrOct14','CommentsOct14',
				'LnOct214','FlOct214','FrOct214','BrOct214','CommentsOct214',
				'Tag2')		

str(asmi14)
head(asmi14)
foo <- asmi14

head(foo)
str(foo)


foo$Site[foo$Plot <= 610 & foo$Plot >= 606] <- "5"
foo$Site[foo$Plot == 8 | foo$Plot == 578 | foo$Plot == 581 | foo$Plot == 799] <- "15"
foo$Site[foo$Plot == 300 | foo$Plot <= 515 & foo$Plot >= 512] <- "19"
foo$Site[foo$Plot == 238 | foo$Plot == 480 | 
		foo$Plot == 598 | foo$Plot == 611 | foo$Plot == 614] <- "26"
foo$Site[foo$Plot == 22 | foo$Plot == 1 | foo$Plot == 32 | foo$Plot == 42
		| foo$Plot == 52 ] <- "Cebolla Mid"

foo$Site[foo$Plot == 16 | foo$Plot == 96 | foo$Plot == 58
		| foo$Plot == 89 | foo$Plot == 98 | foo$Plot == 90 ] <- "Cebolla North"

# foo <- subset(foo, !is.na(foo$Site))

names(foo)
str(foo)
hist(foo$Ln13)
hist(foo$LnApr14[foo$LnApr14 != 0])
hist(foo$LnMay14[foo$LnMay14 != 0])
hist(foo$LnAug14[foo$Site == 'Cebolla North'])

head(foo)
mcapt <- foo[,c("Tag","Plot","Site",
				"LnMay14","LnJun14","LnJun214",
				"LnJul14","LnAug14","LnSep14","LnOct14","LnOct214")]
head(mcapt)
names(mcapt)

# http://www.afsc.noaa.gov/Publications/ProcRpt/PR2013-04.pdf
# http://rpubs.com/brouwern/package_marked
# Laake et al 2013 

## Phi is survival probability. 

HMMLikelihood=function(x,first,m,T,dmat,gamma,delta)
{
# Arguments:
# x: observed sequence (capture (encounter) history)
# first: occasion to start sequence
# m: number of states
# T: number of occasions; sequence length
# dmat: array of occasion specific observation probabilty matrices
# gamma: array of occasion specific transition matrices
# delta: initial state distribution
# Other variables:
# lnl: log likelihood value
# phi: alpha/sum(alpha) sequence as defined in Zucchini/MacDonald
# v: temp variable to hold phi calculations
# u: sum(v)
# Assign prob state vector for initial observation: delta*p(x_first)
v=delta%*%diag(dmat[first,x[first],])
# Compute log-likelihood contribution for first observation; for
# models that condition on first observation u=1,lnl=0
u=sum(v)
phi=v/u
lnl=log(u)
# Loop over occasions for this encounter history (x)
for(t in (first+1):T)
{
# Compute likelihood contribution for this occasion
v=phi%*%gamma[t-1,,]%*%diag(dmat[t,x[t],])
u=sum(v)
lnl=lnl+log(u)
# Compute updated state vector
phi=v/u
}
return(lnl)
}


library(marked)


data(dipper)
head(dipper)
str(dipper)

process.data(dipper, model = "cjs", begin.time = 1981)

#Convert intrayear AsMi to look like dipper data:

vals <- mcapt[,4:11]
vals[is.na(vals)] <- 0
vals[vals>0]<-1
head(vals)

CR <- data.frame(mcapt[,1:3],vals)
head(CR)
str(CR)



# Vector to string
CR <- data.frame(CR[,1:3],ch = apply(CR[,-(1:3)], 1, paste, collapse = ""))
head(CR)
str(CR)

# Add 

# Week transitions
# May 13, June 2, June 25, July 15, Aug 5, Sept 7, Oct 8, Oct 28
wks <- c(20, 23, 26, 29, 32, 36, 41, 44)
wks.int <- diff(wks)

# Add fencing field


setwd(path.expand("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/R_Analysis/R_tables"))
PS<-read.csv("2014_PlotSummary_asmi.csv", header = T, as.is = T)
	

head(PS)
names(PS)

fen <- unique(PS[,c("Plot","Fence")])
names(fen)
names(CR)

CR <- merge(fen, CR, by = "Plot")


# Fruit production
Fr.cols <- grep("Fr", names(foo))
Fruit <- foo[,Fr.cols]
head(Fruit)
Fruit[is.na(Fruit)] <- 0
names(Fruit) <- c("Fr21","Fr22","Fr23","Fr24","Fr25","Fr26","Fr27","Fr28","Fr29")

CR <- cbind(CR,Fruit)
head(CR)

# Length at each interval

unique(CR$ch)

CR <- transform(CR, ch = as.character(ch))

?crm

# remove tags never above ground? Can't find 'min' when all 0
CR1 <- subset(CR, ch != "00000000")
CR1 <- transform(CR1, Fence = as.factor(Fence))

# 1
CR1.pro <- process.data(CR1, model = "HMMSCJS", begin.time = 20, groups = "Fence") 



names(CR1.pro)
CR1.pro$nstrata # need 4, seedling, vegetative, reproductive, dormant???

####2
# How do I set things like fencing or fruiting or important to my system?
dd.asmi <- make.design.data(CR1.pro) #, parameters = list(time.varying = names(Fruit))) 

names(dd.asmi)

crm(CR1.pro, model = "hmmCJS", 
	model.parameters = list(Phi = list(formula = ~1),
			p=list(formula = ~time)))
## p is capture probability over time (weeks over the summer)



design.Phi <- list(static = c("Site"), time.varying = c("Fr"))

model.parameters <- list(Phi=list(formula=~Site + Fr),
				p=list(formula = ~time+Fr))

# Design data with static (fencing?) and time-varying covariates, fruit

design.Phi <- list(static = c('Site'))

design.p <- list(time.varying = c('Fr'))

make.design.data(CR1.pro, list(Phi = design.Phi, p = design.p))

## define the model
Phi.fr <- list(formula = ~ Site)

## define the p (observation) model
p.fr <- list(formula = ~Site+Fr)

## Run model
cjs.initial(list(design.Phi,design.p))

model_mcmc <- crm(CR1.pro, 
			model = "probitCJS",	# won't accumulate records with similar CR histories
							# a MCMC implementation
			begin.time = 20,
			design.parameters = list(Phi = design.Phi, p = design.p),
			model.parameters = list(Phi = Phi.fr, p = p.fr),
			brunin = 1000,
			iter = 10000)



# 2
# time.intervals for unequal time intervals
crm(CR1.pro)
print(crm(CR1.pro))



# annual data

setwd(path.expand("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/R_Analysis/R_tables"))
 asmi <-read.csv("2014_RawData_AsMi.csv", header = T, as.is = T, na.strings = 'na')

names(asmi)
asmi <- subset(asmi, AsMi_site_id != 1 & AsMi_site_id != 2 &
		AsMi_plot_id != 598)	# Should keep Beaver and Cebolla seperate, 
						# remove plot 598
ASMIF4 <- asmi[order(asmi$AsMi_site_id, asmi$AsMi_tag_id, asmi$year),] 
str(asmi)
str(ASMIF4)

table(ASMIF4$status)		# one blank untill get rid of last row of asmi
table(ASMIF4$year, ASMIF4$AsMi_site_id)


#### Set the order of the classes
stages <- c("seedling", "vegetative", "reproductive", "dormant","dead") 
	#Orders the status by the stage levels
  ASMIF4$status <- ordered(ASMIF4$status, levels = stages)	

## check the data for errors
table(ASMIF4$status)
table(ASMIF4$year,ASMIF4$status)
table(ASMIF4$year, ASMIF4$AsMi_site_id)
table(ASMIF4$length)
ASMIF4$status[is.na(ASMIF4$length)]
	subset(ASMIF4, status == 'vegetative' & is.na(length))
ASMIF4$status[is.na(ASMIF4$length) & ASMIF4$status == 'vegetative'] 

str(ASMIF4)

# One length that was recorded crazy
ASMIF4$length[ASMIF4$length == 921] <- NA	
		
# Change insect browsing, mammal browsing, 
#	both to True or False any browsing 
ASMIF4$browsed[ASMIF4$browsing=="None"]<-FALSE		
ASMIF4$browsed[ASMIF4$browsing!="None"]<-TRUE
table(ASMIF4$browsed)
table(ASMIF4$browsing)
mrasmi <- subset(merge(ASMIF4[,1:8,10:11], ASMIF4[,1:8,10:11],
	 by = "AsMi_tag_id", sort = FALSE), year.x == year.y - 1) 	
names(mrasmi)
head(mrasmi)

# annual data long to wide
names(ASMIF4[,c(1:8,10:11)])

fooasmi <- ASMIF4
foo2 <- transform(fooasmi, status = as.character(status))
str(foo2)
foo2$status

foo2$status[foo2$status == "seedling"] <- "S"
foo2$status[foo2$status == "vegetative"] <- "V"
foo2$status[foo2$status == "reproductive"] <- "R"
foo2$status[foo2$status == "dormant"] <- "D"
foo2$status[foo2$status == "dead"] <- "0"
foo2$status[is.na(foo2$status)] <- "0"


longasmi <- dcast(foo2, AsMi_tag_id ~ year, value.var = "status", fill = "0")
str(longasmi)
head(longasmi)

Mark1 <- data.frame(Pl = longasmi[,1], 
	cr = apply(longasmi[,-1], 1, paste, collapse = ""))
head(Mark1)
Mark2 <- ddply(Mark1, .(cr), summarize,
	Freq = length(cr))


head(Mark2)

###
# Frequency of transitions




write.table(Mark2, "clipboard", sep="\t", row.names=FALSE)

longlength <- dcast(ASMIF4, AsMi_tag_id ~ year, value.var = "length", fill = 0)
longfruit <- dcast(ASMIF4, AsMi_tag_id ~ year, value.var = "fruit", fill = 0)

head(longlength)
head(longfruit)


# http://www.doc.govt.nz/Documents/science-and-technical/inventory-monitoring/im-toolbox-herpetofauna-population-estimates.pdf
### www.phidot.org	# help for mark recapture data things
#RMark
library(RMark)




#EasyMARK?

################################################
# Unequal time intervals (multiple times a year??) from Dennis and Ponciano 2014



################################################
# Chen et al. 2009



###################################################
# BaSTA

 


#######################################################
## Format data for MARK

