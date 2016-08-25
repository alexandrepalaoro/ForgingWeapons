rm(list=ls())
##
###PACKAGES NEEDED
##
require(bbmle)
require(betareg)

##LOAD AND REARRANGE DATA
full.table<-read.table("dados-morfo.csv",h=T,sep=',')
names(full.table)
full.table$species<-factor(full.table$species,levels=c("longirostri","denticulata","abtao"))
full.table<-full.table[order(full.table$species),]

## I am reordering the table to make more visible plots, otherwise denticulata's data
## would be hidden behind abtao's data.  

## 
## This is for plotting the legend later on.
##
full.table$sp.plot<-factor(full.table$species,levels=c("Aegla longirostri","Aegla denticulata","Aegla abtao"))

######################################################################################
## Scaling and centering the continuous variables that will be used in the analyses ##
######################################################################################
full.table$lncs.scale<-as.vector(scale(full.table$lncs,center=T,scale=T))
full.table$cc.scale<-as.vector(scale(log(full.table$cc),center=T,scale=T))
full.table$proc.scale<-as.vector(scale(full.table$dist.proc,center=T,scale=T))


#############################
## INVESTMENT IN CLAW SIZE ##
#############################


lm1<-lm(lncs~cc.scale*species*sex,data=full.table)
plot(lm1,which=1)

glm2<-glm(lncs~cc.scale*species*sex,data=full.table,
          family=gaussian(link="log"))
plot(glm2,which=1)

glm3<-glm(lncs~cc.scale*species*sex,data=full.table,
          family=Gamma(link="log"))
plot(glm3,which=1)

## The gamma distribution seems to homonogenize the residuals pretty well.

anova(glm3,test="Chisq")
summary(glm3)

###########################
## PARAMETERS ESTIMATION ##
###########################
#
###1. Model with no intercept to increase standard deviation precision 
#-------------(Schielzeth, Methods Ecol. Evol., 2010)-----------------#
#Note: Do not take into account P-values calculated here, they only 
#say whether the parameters differ from zero. To get comparable P-values
#look in summary(glm3)

glmEST<-glm(lncs~cc.scale*species*sex-1,data=full.table,
            family=Gamma(link="log"))

lncs.coef<-coefficients(glmEST)
summary(glmEST) #getting their placements in the summary table

# Since R used females of AEGLA LONGIROSTRI as their base, begin calculating by them
###
#FEMALE AEGLA LONGIROSTRI
##

exp(lncs.coef[2]) #INTERCEPT

exp(lncs.coef[1]) #SLOPE

##
#MALE AEGLA LONGIROSTRI
##
# Male's intercept is the female's plus the value in the "sexmale" factor

exp(lncs.coef[2]+lncs.coef[5]) #INTERCEPT; P < 0.00001

# Male's slope is the female's slope plus the value in the interaction 
# cc:sexmale

exp(lncs.coef[1]+lncs.coef[8]) #SLOPE; P < 0.0001

##
#FEMALE AEGLA ABTAO
##

## Since we removed the intercept, we get the real value here

exp(lncs.coef[4]) #INTERCEPT

exp(lncs.coef[1]+lncs.coef[7]) #SLOPE

##
#MALE AEGLA ABTAO
##

## Things get really tricky here. Simplifying, the intercept for the male is the
## intercept for the female plus the male factor and the interaction species:factor

exp(lncs.coef[4]+lncs.coef[5]+lncs.coef[10]) #INTERCEPT

## Meanwhile the slope follows the same logic: female slope plus every interaction
## between size and male (i.e. lncs:species, lncs:sex, lncs:spp:sex)

exp(lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+lncs.coef[12]) #SLOPE

##
#FEMALE AEGLA DENTICULATA
##

exp(lncs.coef[3]) #INTERCEPT

exp(lncs.coef[1]+lncs.coef[6]) #SLOPE

##
#MALE AEGLA DENTICULATA
##

exp(lncs.coef[3]+lncs.coef[5]+lncs.coef[9]) #INTERCEPT

exp(lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+lncs.coef[11]) #SLOPE

#####################################
## CONFIDENCE INTERVALS ESTIMATION ##
#####################################

#Looking at the positions in the table

confint(glmEST)

##
#FEMALE AEGLA LONGIROSTRI
##

exp(confint(glmEST))[2,1] #INTERCEPT 2.5%
exp(confint(glmEST))[2,2] #INTERCEPT 97.5%

exp(confint(glmEST))[1,1] #SLOPE 2.5%
exp(confint(glmEST))[1,2] #SLOPE 97.5%

##
#MALE AEGLA LONGIROSTRI
##

# Here it is a bit trickier. The confidence intervals for males are:
# female parameter plus the value of a given confidence interval for
# the male. The same reasoning is applied for the males

exp(lncs.coef[2]+confint(glmEST)[5,1]) #INTERCEPT 2.5%
exp(lncs.coef[2]+confint(glmEST)[5,2]) #INTERCEPT 97.5%

exp(lncs.coef[1]+confint(glmEST)[8,1]) #SLOPE 2.5%
exp(lncs.coef[1]+confint(glmEST)[8,2]) #SLOPE 97.5%

##
#FEMALE AEGLA ABTAO
##

exp(confint(glmEST)[4,1]) #INTERCEPT 2.5%
exp(confint(glmEST)[4,2]) #INTERCEPT 97.5%

exp(lncs.coef[1]+confint(glmEST)[7,1]) #SLOPE 2.5%
exp(lncs.coef[1]+confint(glmEST)[7,2]) #SLOPE 97.5%

##
#MALE AEGLA ABTAO
##

exp(lncs.coef[4]+lncs.coef[5]+confint(glmEST)[10,1]) #INTERCEPT 2.5%
exp(lncs.coef[4]+lncs.coef[5]+confint(glmEST)[10,2]) #INTERCEPT 97.5%

exp(lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+confint(glmEST)[12,1]) #SLOPE 2.5%
exp(lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+confint(glmEST)[12,2]) #SLOPE 97.5%

##
#FEMALE AEGLA DENTICULATA
##

exp(confint(glmEST)[3,1]) #INTERCEPT 2.5%
exp(confint(glmEST)[3,2]) #INTERCEPT 97.5%

exp(lncs.coef[1]+confint(glmEST)[6,1]) #SLOPE 2.5%
exp(lncs.coef[1]+confint(glmEST)[6,2]) #SLOPE 97.5%

##
#MALE AEGLA DENTICULATA
##

exp(lncs.coef[3]+lncs.coef[5]+confint(glmEST)[9,1]) #INTERCEPT 2.5%
exp(lncs.coef[3]+lncs.coef[5]+confint(glmEST)[9,2]) #INTERCEPT 97.5%

exp(lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+confint(glmEST)[11,1]) #SLOPE 2.5%
exp(lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+confint(glmEST)[11,2]) #SLOPE 97.5%

#######################################
## PLOTTING DATA AND PREDICTED LINES ##
#######################################

## We scaled and centered the log of the cephalothorax length.
## Thus, we need to plot the log of the CL.

full.table$cc.log<-log(full.table$cc)

## For plotting we use the unscaled and uncentered variable to
## know how many units of centroid size grows when one unit of 
## CL grows

glmPlot<-glm(lncs~cc.log*species*sex,data=full.table,
             family=Gamma(link="log"))

######
## Creating new x-axis values for each species
#####
## This is needed to use the predict() function

cc.dentM<-seq(min(full.table$cc.log[full.table$species=="denticulata"&full.table$sex=="male"]),
              max(full.table$cc.log[full.table$species=="denticulata"&full.table$sex=="male"]),
              length.out=100)

cc.abtaoM<-seq(min(full.table$cc.log[full.table$species=="abtao"&full.table$sex=="male"]),
               max(full.table$cc.log[full.table$species=="abtao"&full.table$sex=="male"]),
               length.out=100)

cc.longM<-seq(min(full.table$cc.log[full.table$species=="longirostri"&full.table$sex=="male"]),
              max(full.table$cc.log[full.table$species=="longirostri"&full.table$sex=="male"]),
              length.out=100)

cc.dentF<-seq(min(full.table$cc.log[full.table$species=="denticulata"&full.table$sex=="female"]),
              max(full.table$cc.log[full.table$species=="denticulata"&full.table$sex=="female"]),
              length.out=100)

cc.abtaoF<-seq(min(full.table$cc.log[full.table$species=="abtao"&full.table$sex=="female"]),
               max(full.table$cc.log[full.table$species=="abtao"&full.table$sex=="female"]),
               length.out=100)

cc.longF<-seq(min(full.table$cc.log[full.table$species=="longirostri"&full.table$sex=="female"]),
              max(full.table$cc.log[full.table$species=="longirostri"&full.table$sex=="female"]),
              length.out=100)

## Now the factors

sp.dent<-factor(rep("denticulata",length(cc.longM)))
sp.abtao<-factor(rep("abtao",length(cc.abtaoM)))
sp.long<-factor(rep("longirostri",length(cc.dentM)))
malesNew<-factor(rep("male",length(cc.longM)))
femalesNew<-factor(rep("female",length(cc.longM)))

## Now we can start plotting. If you want to save the graph in tiff
## in your computer, just remove the # in the tiff() below
## and in dev.off() at the bottom

#tiff(file="LnCSxCL.tiff",units="mm",width=170,height=150,res=600,
#     compression="lzw")
plot(lncs~cc.log,pch=c(24,21)[as.numeric(sex)],
     bg=adjustcolor(c("red","gold","gray")[as.numeric(species)],0.7),
     cex=1.5,data=full.table,bty='l',col=NA,las=1,ylab="Centroid size (log)",
     xlab="Cephalothorax length (log)")

## Rather awkward, but I prefer to build figure pannels in photoshop
## So I plot everything here and move the legends around in PS

legend("topleft", legend=c(expression(italic("Aegla abtao"),
                                      italic("Aegla denticulata"),
                                      italic("Aegla longirostri")),
                           "Male","Female","M"),
       pch = c(24,24,24,21,21,21),col=NA,
       pt.bg=adjustcolor(c("red","gold","gray","red","gold","gray"),0.7),bty='n',
       cex=1.2,ncol=2)
legend(3.45,6.65,legend="(a)",cex=1.2,bty='n')

## Plotting the lines 

## First we build a data.frame with the new data and then plot the lines. 
## Over the years I found that data.frames give less error messages than 
## using list(). For instance, the 'betareg' package does not work really
##well if you use list(). Thus, data.frames it is. 

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(cc.log=cc.dentM,species=sp.dent,sex=malesNew)
lines(dentM$cc.log,predict(glmPlot,type='response',newdata=dentM),lwd=3,col="yellow4")

##
#AEGLA DENTICULATA - FEMALES
##
dentF<-data.frame(cc.log=cc.dentF,species=sp.dent,sex=femalesNew)
lines(dentF$cc.log,predict(glmPlot,type='response',newdata=dentF),lwd=3,lty=2,col="yellow4")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(cc.log=cc.abtaoM,species=sp.abtao,sex=malesNew)
lines(abtaoM$cc.log,predict(glmPlot,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA ABTAO - FEMALES
##

abtaoF<-data.frame(cc.log=cc.abtaoF,species=sp.abtao,sex=femalesNew)
lines(abtaoF$cc.log,predict(glmPlot,type='response',newdata=abtaoF),col="black",lwd=3,lty=2)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(cc.log=cc.longM,species=sp.long,sex=malesNew)
lines(longM$cc.log,predict(glmPlot,type='response',newdata=longM),col="darkred",lwd=3)

##
#AEGLA LONGIROSTRI - FEMALES
##

longF<-data.frame(cc.log=cc.longF,species=sp.long,sex=femalesNew)
lines(longF$cc.log,predict(glmPlot,type='response',newdata=longF),col="darkred",lwd=3,lty=2)
#dev.off()

#DONE! :D
