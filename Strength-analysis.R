rm(list=ls())
#
###
### PACKAGED NEEDED
###
#

library(nlme)

##LOAD AND REARRANGE DATA

full.table<-read.table("morpho-data.csv",h=T,sep=',')
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

########################################################################
## Scaling the continuous variables that will be used in the analyses ##
########################################################################

full.table$lncs.scale<-as.vector(scale(full.table$lncs,center=T,scale=T))
full.table$cc.scale<-as.vector(scale(log(full.table$cc),center=T,scale=T))
full.table$proc.scale<-as.vector(scale(full.table$dist.proc,center=T,scale=T))

#################################
## INVESTMENT IN CLAW STRENGTH ##
#################################

## This is a hard response variable to analyze. 
## First, let's look a the plot:

plot(apodeme~lncs.scale,data=full.table,cex=1.5,las=1,ylab="Apodeme area (square root)",
     xlab="Centroid size (log)")

## It does not seem that species differ in their investment. Let's color this graph up

plot(apodeme~lncs.scale,pch=c(24,21)[as.numeric(sex)],
     bg=adjustcolor(c("red","gold","gray")[as.numeric(species)],0.7),
     cex=1.5,data=full.table,bty='l',col=NA,las=1,ylab="Apodeme area (square root)",
     xlab="Centroid size (log)")

## By looking at it, I can only see differences in sheer claw size, not in strength
## But maybe there are small differences in intercepts and slopes that we cannot see.
## Since there are values really close to one, a Gamma distribution with a log link 
## should work.

glm1<-glm(apodeme~lncs.scale*species*sex,data=full.table,
          family=Gamma(link="log"))

summary(glm1)

## Hmm, the coefficientes for denticulata seems way too big for that region of the graph
## Let's check it out.

exp(coefficients(glm1)[2]+coefficients(glm1)[6]) # Seems way too high. Let's compare.
exp(coefficients(glm1)[2]+coefficients(glm1)[8]) # Higher than a fighting species, that's just wrong

## Plotting the predicted line just to get an idea of what is going on

pred1<-predict(glm1,type="response")[full.table$species=="denticulata"&full.table$sex=="female"]
lines(full.table$lncs.scale[full.table$species=="denticulata"&full.table$sex=="female"],pred1,lwd=2)

## If you take a careful look, you will notice that by the end of the line it start going up
## (as an exponential should do). Apparently, this property of the exponential is making parameter
## estimation unrealiable. To put that unrealiability to the test, I will do one thing here:

plot(log(apodeme)~lncs,data=full.table)

## Please note the first data points. They all belong to denticulata females and they show a large
## increase. That's because a log() of a number near zero is a large value, and even very close numbers
## will get large values. If you want further tests, check the confidence intervals 
## (code further down for the real analysis, just add exp()) - you will notice that all intervals overlap
## This means that there should be no difference among species or sexes and the model is unstable.
## I have tried several other distributions and data transformation that I am not going to show here,
## but you can try changing for a gaussian distribution, an identity link, a lm() and you will get the  
## same results.

## To get around this issue I square-rooted the apodeme value. Why? Beucase I measured the apodeme area, 
## which is a quadratic measure. By square-rooting it, I transformed it from mm² to mm. Now I can use
## a linear model and get a stable model with reasonable parameters. Let the coding begin. 

# Note: I used a gls() because allometric data usually show auto-correlation and it seemed reasonable to
# control for that. The residuals also seem pleased with my choice.

gls1<-gls(sqrt(apodeme)~lncs.scale*species*sex,data=full.table,
          weights=varComb(varExp(form=~lncs.scale),varIdent(form=~1|species)))

## Plotting residuals against fitted values
RESID1<-resid(gls1,type="normalized")
FIT1<-fitted(gls1)
plot(FIT1,RESID1, xlab="Fitted",ylab="Residuals")

## Residuals seem well dispersed.

summary(gls1)
anova(gls1)

###########################
## PARAMETERS ESTIMATION ##
###########################
#
###1. Model with no intercept to increase standard deviation precision 
#-------------(Schielzeth, Methods Ecol. Evol., 2010)-----------------#
#Note: Do not take into account P-values calculated here, they only 
#say whether the parameters differ from zero. To get comparable P-values
#look in summary(glm5)

glsEST<-gls(sqrt(apodeme)~lncs.scale*species*sex-1,data=full.table,
            weights=varComb(varExp(form=~lncs.scale),varIdent(form=~1|species)))

lncs.coef<-coefficients(glsEST)

summary(glsEST) #getting their placements in the summary table

#Since R used females of AEGLA LONGIROSTRI as their base, begin calculating by them

##
#FEMALE AEGLA LONGIROSTRI
##

lncs.coef[2] #INTERCEPT

lncs.coef[1] #SLOPE

##
#MALE AEGLA LONGIROSTRI
##

# Male's intercept is the female's plus the value in the "sexmale" factor

lncs.coef[2]+lncs.coef[5] #INTERCEPT

# Male's slope is the female's slope plus the value in the interaction 
# cc:sexmale

lncs.coef[1]+lncs.coef[8] #SLOPE

##
#FEMALE AEGLA ABTAO
##

## Since we removed the intercept, we get the real value here

lncs.coef[4] #INTERCEPT

lncs.coef[1]+lncs.coef[7] #SLOPE

##
#MALE AEGLA ABTAO
##

lncs.coef[4]+lncs.coef[5]+lncs.coef[10] #INTERCEPT

lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+lncs.coef[12] #SLOPE

##
#FEMALE AEGLA DENTICULATA
##

lncs.coef[3] #INTERCEPT

lncs.coef[1]+lncs.coef[6] #SLOPE

##
#MALE AEGLA DENTICULATA
##

lncs.coef[3]+lncs.coef[9] #INTERCEPT

lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+lncs.coef[11] #SLOPE

#####################################
## CONFIDENCE INTERVALS ESTIMATION ##
#####################################

#Looking at the positions in the table

confint(glmEST)

##
#FEMALE AEGLA LONGIROSTRI
##

confint(glsEST)[2,1] #INTERCEPT 2.5%
confint(glsEST)[2,2] #INTERCEPT 97.5%

confint(glsEST)[1,1] #SLOPE 2.5%
confint(glsEST)[1,2] #SLOPE 97.5%

##
#MALE AEGLA LONGIROSTRI
##

# Here it is a bit trickier. The confidence intervals for males are:
# female parameter plus the value of a given confidence interval for
# the male. The same reasoning is applied for the males

lncs.coef[2]+confint(glsEST)[5,1] #INTERCEPT 2.5%
lncs.coef[2]+confint(glsEST)[5,2] #INTERCEPT 97.5%

lncs.coef[1]+confint(glsEST)[8,1] #SLOPE 2.5%
lncs.coef[1]+confint(glsEST)[8,2] #SLOPE 97.5%

##
#FEMALE AEGLA ABTAO
##

confint(glsEST)[4,1] #INTERCEPT 2.5%
confint(glsEST)[4,2] #INTERCEPT 97.5%

lncs.coef[1]+confint(glsEST)[7,1] #SLOPE 2.5%
lncs.coef[1]+confint(glsEST)[7,2] #SLOPE 97.5%

##
#MALE AEGLA ABTAO
##

lncs.coef[4]+lncs.coef[5]+confint(glsEST)[10,1] #INTERCEPT 2.5%
lncs.coef[4]+lncs.coef[5]+confint(glsEST)[10,2] #INTERCEPT 97.5%

lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+confint(glsEST)[12,1] #SLOPE 2.5%
lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+confint(glsEST)[12,2] #SLOPE 97.5%

##
#FEMALE AEGLA DENTICULATA
##

confint(glsEST)[3,1] #INTERCEPT 2.5%
confint(glsEST)[3,2] #INTERCEPT 97.5%

lncs.coef[1]+confint(glsEST)[6,1] #SLOPE 2.5%
lncs.coef[1]+confint(glsEST)[6,2] #SLOPE 97.5%

##
#MALE AEGLA DENTICULATA
##

lncs.coef[3]+confint(glsEST)[9,1] #INTERCEPT 2.5%
lncs.coef[3]+confint(glsEST)[9,2] #INTERCEPT 97.5%

lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+confint(glsEST)[11,1] #SLOPE 2.5%
lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+confint(glsEST)[11,2] #SLOPE 97.5%

#######################################
## PLOTTING DATA AND PREDICTED LINES ##
#######################################

## For plotting we use the unscaled and uncentered variable to
## know how many units of centroid size grows when one unit of 
## LnCS grows

glsPlot<-gls(sqrt(apodeme)~lncs*species*sex,data=full.table,
             weights=varComb(varExp(form=~lncs),varIdent(form=~1|species)))

######
## Creating new x-axis values for each species
#####
## This is needed to use the predict() function

lncs.dentM<-seq(min(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="male"]),
                max(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="male"]),length.out=100)

lncs.abtaoM<-seq(min(full.table$lncs[full.table$species=="abtao"&full.table$sex=="male"]),
                 max(full.table$lncs[full.table$species=="abtao"&full.table$sex=="male"]),length.out=100)

lncs.longM<-seq(min(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="male"]),
                max(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="male"]),length.out=100)

lncs.dentF<-seq(min(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="female"]),
                max(full.table$lncs[full.table$species=="denticulata"&full.table$sex=="female"]),length.out=100)

lncs.abtaoF<-seq(min(full.table$lncs[full.table$species=="abtao"&full.table$sex=="female"]),
                 max(full.table$lncs[full.table$species=="abtao"&full.table$sex=="female"]),length.out=100)

lncs.longF<-seq(min(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="female"]),
                max(full.table$lncs[full.table$species=="longirostri"&full.table$sex=="female"]),length.out=100)

## Now the factors

sp.dent<-factor(rep("denticulata",length(lncs.dentM)))
sp.abtao<-factor(rep("abtao",length(lncs.abtaoM)))
sp.long<-factor(rep("longirostri",length(lncs.longM)))
malesNew<-factor(rep("male",length(lncs.longM)))
femalesNew<-factor(rep("female",length(lncs.longM)))


## Now we can start plotting. If you want to save the graph in tiff
## format in our computer, please remove the # in the tiff() below
## and in dev.off() at the bottom

#tiff(file="apodxLnCS.tiff",units="mm",width=170,height=150,res=600,
#     compression="lzw")
plot(sqrt(apodeme)~lncs,pch=c(24,21)[as.numeric(sex)],
     bg=adjustcolor(c("red","gold","gray")[as.numeric(species)],0.7),
     cex=1.5,data=full.table,bty='l',col=NA,las=1,ylab="Apodeme area (square root)",
     xlab="Centroid size (log)")

## Rather awkward, but I never got around to learn how to plot the triangles and circles
## side-by-side without adding another whole column. Therefore, I plot everything here,
## open it in Photoshop and move legends around there.
## If anyone knows how to work around my legends issue, please let me know. It would be of 
## tremendous help!!

legend("topleft", legend=c(expression(italic("Aegla abtao"),
                                      italic("Aegla denticulata"),
                                      italic("Aegla longirostri")),
                           "Male","Female","M"),
       pch = c(24,24,24,21,21,21),col=NA,
       pt.bg=adjustcolor(c("red","gold","gray","red","gold","gray"),0.7),bty='n',
       cex=1.2,ncol=2)
legend(6.5,1.3,legend="(b)",cex=1.2,bty='n')

## Plotting the lines 

## First we build a data.frame with the new data and then plot the lines. 
## Over the years I found that data.frames give less error messages than 
## using list(). For instance, the 'betareg' package does not work really
## well if you use list(). Thus, data.frames it is. 

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(lncs=lncs.dentM,species=sp.dent,sex=malesNew)
lines(dentM$lncs,predict(glsPlot,type='response',newdata=dentM),lwd=3,col="yellow4")

##
#AEGLA DENTICULATA - FEMALES
##

dentF<-data.frame(lncs=lncs.dentF,species=sp.dent,sex=femalesNew)
lines(dentF$lncs,predict(glsPlot,type='response',newdata=dentF),lwd=3,lty=2,col="yellow4")

lines(dentM$lncs,predict(glsPlot,type='response',newdata=dentM),lwd=3,col="yellow4")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs=lncs.abtaoM,species=sp.abtao,sex=malesNew)
lines(abtaoM$lncs,predict(glsPlot,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA ABTAO - FEMALES
##

abtaoF<-data.frame(lncs=lncs.abtaoF,species=sp.abtao,sex=femalesNew)
lines(abtaoF$lncs,predict(glsPlot,type='response',newdata=abtaoF),col="black",lwd=3,lty=2)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs=lncs.longM,species=sp.long,sex=malesNew)
lines(longM$lncs,predict(glsPlot,type='response',newdata=longM),col="darkred",lwd=3)

##
#AEGLA LONGIROSTRI - FEMALES
##

longF<-data.frame(lncs=lncs.longF,species=sp.long,sex=femalesNew)
lines(longF$lncs,predict(glsPlot,type='response',newdata=longF),col="darkred",lwd=3,lty=2)

#dev.off()

#DONE :D
