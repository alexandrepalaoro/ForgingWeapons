rm(list=ls())
##
###PACKAGES NEEDED
##
require(betareg)
require(lmtest)

##LOAD AND REARRANGE DATA

full.table<-read.table("dados-morfo1.csv",h=T,sep=',')
names(full.table)
full.table$species<-factor(full.table$species,levels=c("longirostri","denticulata","abtao"))
full.table<-full.table[order(full.table$species),]

#
## I am reordering the table to make more visible plots, otherwise denticulata's data
## would be hidden behind abtao's data.  

## 
## This is for legend plotting only
##

full.table$sp.plot<-factor(full.table$species,levels=c("Aegla longirostri","Aegla denticulata","Aegla abtao"))

######################################################################################
## Scaling and centering the continuous variables that will be used in the analyses ##
######################################################################################

full.table$lncs.scale<-as.vector(scale(full.table$lncs,center=T,scale=T))
full.table$cc.scale<-as.vector(scale(log(full.table$cc),center=T,scale=T))
full.table$proc.scale<-as.vector(scale(full.table$dist.proc,center=T,scale=T))

###################################
## INVESTMENT IN CLAW EFFICIENCY ##
###################################

beta1<-betareg(maT~lncs.scale*species*sex | species*sex,data=full.table)
summary(beta1)
plot(beta1,which=1)

## Since anova() does not work with betareg() we will build three extra nested models
## and compare them using a likelihood ratio test with 'lrtest' (in lmtest package)

# Building models

beta2<-betareg(maT~lncs.scale*species+species*sex+lncs.scale*sex | species*sex,data=full.table)
beta3<-betareg(maT~lncs.scale*species+species*sex+sex | species*sex,data=full.table)
beta4<-betareg(maT~lncs.scale*species+species+lncs.scale*sex | species*sex,data=full.table)

# Comparing between them

lrtest(beta1,beta2,beta3,beta4)

###########################
## PARAMETERS ESTIMATION ##
###########################
#
###1. Model with no intercept to increase standard deviation precision 
#-------------(Schielzeth, Methods Ecol. Evol., 2010)-----------------#
###2. Beta regression uses a logit link to bound the response variable
###   between 0 and 1. This means that, to get estimates within the 
###   boundaries of the response variable. The inverse of the logit
###   is: exponential of the linear predictor divided by one plus the 
###   exponential of the linear predictor. In S+: 
###              exp(coefficient)/(1+exp(coefficient))

#Note: Do not take into account P-values calculated in summary(betaEST) below,  
#they only say whether the parameters differ from zero. To get comparable 
#P-values look summary(glm5)

betaEST<-betareg(maT~lncs.scale*species*sex-1 | species*sex,data=full.table)

beta.coef<-coefficients(betaEST)
summary(betaEST) #getting their placements in the summary table

# Since R used females of AEGLA LONGIROSTRI as their base, begin calculating by them
###
#FEMALE AEGLA LONGIROSTRI
##

exp(beta.coef[2])/(1+exp(beta.coef[2])) #INTERCEPT

exp(beta.coef[1])/(1+exp(beta.coef[1])) #SLOPE

##
#MALE AEGLA LONGIROSTRI
##
# Male's intercept is the female's plus the value in the "sexmale" factor

exp(beta.coef[2]+beta.coef[5])/(1+exp(beta.coef[2]+beta.coef[5])) #INTERCEPT

# Male's slope is the female's slope plus the value in the interaction 
# cc:sexmale

exp(beta.coef[1]+beta.coef[8])/(1+exp(beta.coef[1]+beta.coef[8])) #SLOPE

##
#FEMALE AEGLA ABTAO
##

## Since we removed the intercept, we get the real value here, no need to add

exp(beta.coef[4])/(1+exp(beta.coef[4])) #INTERCEPT

exp(beta.coef[1]+beta.coef[7])/(1+exp(beta.coef[1]+beta.coef[7])) #SLOPE

##
#MALE AEGLA ABTAO
##

## Things get really tricky here. Simplifying, the intercept for the male is the
## intercept for the female plus the male factor and the interaction species:factor

exp(beta.coef[4]+beta.coef[5]+beta.coef[10])/(1+exp(beta.coef[4]+beta.coef[5]+beta.coef[10])) #INTERCEPT

## Meanwhile the slope follows the same logic: female slope plus every interaction
## between size and male (i.e. lncs:species, lncs:sex, lncs:spp:sex)

exp(beta.coef[1]+beta.coef[7]+beta.coef[8]+beta.coef[12])/(1+exp(beta.coef[1]+beta.coef[7]+beta.coef[8]+beta.coef[12])) #SLOPE

##
#FEMALE AEGLA DENTICULATA
##

exp(beta.coef[3])/(1+exp(beta.coef[3])) #INTERCEPT

exp(beta.coef[1]+beta.coef[6])/(1+exp(beta.coef[1]+beta.coef[6])) #SLOPE

##
#MALE AEGLA DENTICULATA
##

exp(beta.coef[3]+beta.coef[5]+beta.coef[9])/(1+exp(beta.coef[3]+beta.coef[5]+beta.coef[9])) #INTERCEPT

exp(beta.coef[1]+beta.coef[6]+beta.coef[8]+beta.coef[11])/(1+exp(beta.coef[1]+beta.coef[6]+beta.coef[8]+beta.coef[11])) #SLOPE


#####################################
## CONFIDENCE INTERVALS ESTIMATION ##
#####################################

#Looking at their positions

confint(betaEST)

##
#FEMALE AEGLA LONGIROSTRI
##

exp(confint(betaEST))[2,1]/(1+exp(confint(betaEST))[2,1]) #INTERCEPT 2.5%
exp(confint(betaEST))[2,2]/(1+exp(confint(betaEST))[2,2]) #INTERCEPT 97.5%

exp(confint(betaEST))[1,1]/(1+exp(confint(betaEST))[1,1]) #SLOPE 2.5%
exp(confint(betaEST))[1,2]/(1+exp(confint(betaEST))[1,2]) #SLOPE 97.5%

##
#MALE AEGLA LONGIROSTRI
##

# Here it is a bit trickier. The confidence intervals for males are:
# female parameter plus the value of a given confidence interval for
# the male.

exp(beta.coef[2]+confint(betaEST)[5,1])/(1+exp(beta.coef[2]+confint(betaEST)[5,1])) #INTERCEPT 2.5%
exp(beta.coef[2]+confint(betaEST)[5,2])/(1+exp(beta.coef[2]+confint(betaEST)[5,2])) #INTERCEPT 97.5%

exp(beta.coef[2]+confint(betaEST)[8,1])/(1+exp(beta.coef[2]+confint(betaEST)[8,1])) #SLOPE 2.5%
exp(beta.coef[2]+confint(betaEST)[8,2])/(1+exp(beta.coef[2]+confint(betaEST)[8,2])) #SLOPE 97.5%

##
#FEMALE AEGLA ABTAO
##

exp(confint(betaEST)[4,1])/(1+exp(confint(betaEST)[4,1])) #INTERCEPT 2.5%
exp(confint(betaEST)[4,2])/(1+exp(confint(betaEST)[4,2])) #INTERCEPT 97.5%

exp(beta.coef[1]+confint(betaEST)[7,1])/(1+exp(beta.coef[1]+confint(betaEST)[7,1])) #SLOPE 2.5%
exp(beta.coef[1]+confint(betaEST)[7,2])/(1+exp(beta.coef[1]+confint(betaEST)[7,2])) #SLOPE 97.5%

##
#MALE AEGLA ABTAO
##

exp(beta.coef[4]+beta.coef[5]+confint(betaEST)[10,1])/(1+exp(beta.coef[4]+beta.coef[5]+confint(betaEST)[10,1])) #INTERCEPT 2.5%
exp(beta.coef[4]+beta.coef[5]+confint(betaEST)[10,2])/(1+exp(beta.coef[4]+beta.coef[5]+confint(betaEST)[10,2])) #INTERCEPT 97.5%


exp(beta.coef[1]+beta.coef[7]+beta.coef[8]+confint(betaEST)[12,1])/(1+exp(beta.coef[1]+beta.coef[7]+beta.coef[8]+confint(betaEST)[12,1])) #SLOPE 2.5%
exp(beta.coef[1]+beta.coef[7]+beta.coef[8]+confint(betaEST)[12,2])/(1+exp(beta.coef[1]+beta.coef[7]+beta.coef[8]+confint(betaEST)[12,2])) #SLOPE 97.5%

##
#FEMALE AEGLA DENTICULATA
##

exp(confint(betaEST)[3,1])/(1+exp(confint(betaEST)[3,1])) #INTERCEPT 2.5%
exp(confint(betaEST)[3,2])/(1+exp(confint(betaEST)[3,2])) #INTERCEPT 97.5%

exp(beta.coef[1]+confint(betaEST)[6,1])/(1+exp(beta.coef[1]+confint(betaEST)[6,1])) #SLOPE 2.5%
exp(beta.coef[1]+confint(betaEST)[6,2])/(1+exp(beta.coef[1]+confint(betaEST)[6,2])) #SLOPE 97.5%

##
#MALE AEGLA DENTICULATA
##

exp(beta.coef[3]+beta.coef[5]+confint(betaEST)[9,1])/(1+exp(beta.coef[3]+beta.coef[5]+confint(betaEST)[9,1])) #INTERCEPT 2.5%
exp(beta.coef[3]+beta.coef[5]+confint(betaEST)[9,2])/(1+exp(beta.coef[3]+beta.coef[5]+confint(betaEST)[9,2])) #INTERCEPT 97.5%

exp(beta.coef[1]+beta.coef[6]+beta.coef[8]+confint(betaEST)[11,1])/(1+exp(beta.coef[1]+beta.coef[6]+beta.coef[8]+confint(betaEST)[11,1])) #SLOPE 2.5%
exp(beta.coef[1]+beta.coef[6]+beta.coef[8]+confint(betaEST)[11,2])/(1+exp(beta.coef[1]+beta.coef[6]+beta.coef[8]+confint(betaEST)[11,2])) #SLOPE 97.5%

#######################################
## PLOTTING DATA AND PREDICTED LINES ##
#######################################

## For plotting we use the unscaled and uncentered variable to
## know how many units of centroid size grows when one unit of 
## LnCS grows

betaPlot<-betareg(maT~lncs*species*sex | species*sex,data=full.table)

######
## Creating new x-axis values for each species
#####
## This is needed to use the predict() function

lncs.dentM<-seq(min(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="male"]),
                max(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="male"]),
                length.out=100)

lncs.abtaoM<-seq(min(full.table$lncs[full.table$species=="abtao"&full.table$sex=="male"]),
                 max(full.table$lncs[full.table$species=="abtao"&full.table$sex=="male"]),
                 length.out=100)

lncs.longM<-seq(min(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="male"]),
                max(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="male"]),
                length.out=100)

lncs.dentF<-seq(min(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="female"]),
                max(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="female"]),
                length.out=100)

lncs.abtaoF<-seq(min(full.table$lncs[full.table$species=="abtao"&full.table$sex=="female"]),
                 max(full.table$lncs[full.table$species=="abtao"&full.table$sex=="female"]),
                 length.out=100)

lncs.longF<-seq(min(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="female"]),
                max(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="female"]),
                length.out=100)

## Now the factors

sp.dent<-factor(rep("denticulata",length(lncs.dentM)))
sp.abtao<-factor(rep("abtao",length(lncs.abtaoM)))
sp.long<-factor(rep("longirostri",length(lncs.longM)))
malesNew<-factor(rep("male",length(lncs.longM)))
femalesNew<-factor(rep("female",length(lncs.longF)))

## Now we can start plotting. If you want to save the graph in tiff
## format in our computer, please remove the # in the tiff() below
## and in dev.off() at the bottom

#tiff(file="MAxLnCS.tiff",units="mm",width=170,height=150,res=600,
#     compression="lzw")
plot(maT~lncs,pch=c(24,21)[as.numeric(sex)],
     bg=adjustcolor(c("red","gold","gray")[as.numeric(species)],0.7),
     cex=1.5,data=full.table,bty='l',col=NA,las=1,ylab="Mechanical advantage",
     xlab="Centroid size (log)")

## Rather awkward, but I never got around to learn how to plot the triangles and circles
## side-by-side without adding another whole column. Therefore, I plot everything here,
## open it in Photoshop and move legends around there.
## If anyone knows how to work around my legends issue, please let me know. It would be of 
## tremendous help!!

legend("bottomright", legend=c(expression(italic("Aegla abtao"),
                                          italic("Aegla denticulata"),
                                          italic("Aegla longirostri")),
                               "Male","Female","M"),
       pch = c(24,24,24,21,21,21),col=NA,
       pt.bg=adjustcolor(c("red","gold","gray","red","gold","gray"),0.7),bty='n',
       cex=1.2)
legend("topright",legend="(c)",cex=1.2,bty='n')

## Plotting the lines 

## First we build a data.frame with the new data and then plot the lines. 
## Over the years I found that data.frames give less error messages than 
## using list(). For instance, the 'betareg' package does not work really
## well if you use list(). Thus, data.frames it is. 

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(lncs=lncs.dentM,species=sp.dent,sex=malesNew)
lines(dentM$lncs,predict(betaPlot,type='response',newdata=dentM),lwd=3,col="yellow4")

##
#AEGLA DENTICULATA - FEMALES
##

dentF<-data.frame(lncs=lncs.dentF,species=sp.dent,sex=femalesNew)
lines(dentF$lncs,predict(betaPlot,type='response',newdata=dentF),lwd=3,lty=2,col="yellow4")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs=lncs.abtaoM,species=sp.abtao,sex=malesNew)
lines(abtaoM$lncs,predict(betaPlot,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA ABTAO - FEMALES
##

abtaoF<-data.frame(lncs=lncs.abtaoF,species=sp.abtao,sex=femalesNew)
lines(abtaoF$lncs,predict(betaPlot,type='response',newdata=abtaoF),col="black",lwd=3,lty=2)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs=lncs.longM,species=sp.long,sex=malesNew)
lines(longM$lncs,predict(betaPlot,type='response',newdata=longM),col="darkred",lwd=3)

##
#AEGLA LONGIROSTRI - FEMALES
##

longF<-data.frame(lncs=lncs.longF,species=sp.long,sex=femalesNew)
lines(longF$lncs,predict(betaPlot,type='response',newdata=longF),col="darkred",lwd=3,lty=2)

#dev.off()

#DONE :D
