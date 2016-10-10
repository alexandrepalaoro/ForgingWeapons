rm(list=ls())
#
###
### PACKAGES NEEDED
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

#######################################
## INTERPLAY CLAW STRENGTH AND SHAPE ##
######################################

gls2<-gls(sqrt(apodeme)~proc.scale*species*sex,weights=varIdent(form=~1|species*sex),
            data=full.table)

## Plotting residuals against fitted values

RESID1<-resid(gls2,type="normalized")
FIT1<-fitted(gls2)
plot(FIT1,RESID1, xlab="Fitted",ylab="Residuals")

## Residuals seem well dispersed.

summary(gls2)
anova(gls2)

###########################
## PARAMETERS ESTIMATION ##
###########################
#
###1. Model with no intercept to increase standard deviation precision 
#-------------(Schielzeth, Methods Ecol. Evol., 2010)-----------------#
#Note: Do not take into account P-values calculated here, they only 
#say whether the parameters differ from zero. To get comparable P-values
#look in summary(gls2)

gls2EST<-gls(sqrt(apodeme)~proc.scale*species*sex-1,data=full.table,
            weights=varIdent(form=~1|species*sex))

lncs.coef<-coefficients(gls2EST)

summary(gls2EST) #getting their placements in the summary table

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

confint(gls2EST)

##
#FEMALE AEGLA LONGIROSTRI
##

confint(gls2EST)[2,1] #INTERCEPT 2.5%
confint(gls2EST)[2,2] #INTERCEPT 97.5%

confint(gls2EST)[1,1] #SLOPE 2.5%
confint(gls2EST)[1,2] #SLOPE 97.5%

##
#MALE AEGLA LONGIROSTRI
##

# Here it is a bit trickier. The confidence intervals for males are:
# female parameter plus the value of a given confidence interval for
# the male. The same reasoning is applied for the males

lncs.coef[2]+confint(gls2EST)[5,1] #INTERCEPT 2.5%
lncs.coef[2]+confint(gls2EST)[5,2] #INTERCEPT 97.5%

lncs.coef[1]+confint(gls2EST)[8,1] #SLOPE 2.5%
lncs.coef[1]+confint(gls2EST)[8,2] #SLOPE 97.5%

##
#FEMALE AEGLA ABTAO
##

confint(gls2EST)[4,1] #INTERCEPT 2.5%
confint(gls2EST)[4,2] #INTERCEPT 97.5%

lncs.coef[1]+confint(gls2EST)[7,1] #SLOPE 2.5%
lncs.coef[1]+confint(gls2EST)[7,2] #SLOPE 97.5%

##
#MALE AEGLA ABTAO
##

lncs.coef[4]+lncs.coef[5]+confint(gls2EST)[10,1] #INTERCEPT 2.5%
lncs.coef[4]+lncs.coef[5]+confint(gls2EST)[10,2] #INTERCEPT 97.5%

lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+confint(gls2EST)[12,1] #SLOPE 2.5%
lncs.coef[1]+lncs.coef[7]+lncs.coef[8]+confint(gls2EST)[12,2] #SLOPE 97.5%
##
#FEMALE AEGLA DENTICULATA
##

confint(gls2EST)[3,1] #INTERCEPT 2.5%
confint(gls2EST)[3,2] #INTERCEPT 97.5%

lncs.coef[1]+confint(gls2EST)[6,1] #SLOPE 2.5%
lncs.coef[1]+confint(gls2EST)[6,2] #SLOPE 97.5%

##
#MALE AEGLA DENTICULATA
##

lncs.coef[3]+confint(gls2EST)[9,1] #INTERCEPT 2.5%
lncs.coef[3]+confint(gls2EST)[9,2] #INTERCEPT 97.5%

lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+confint(gls2EST)[11,1] #SLOPE 2.5%
lncs.coef[1]+lncs.coef[6]+lncs.coef[8]+confint(gls2EST)[11,2] #SLOPE 97.5%

#######################################
## PLOTTING DATA AND PREDICTED LINES ##
#######################################


## For plotting we use the unscaled and uncentered variable to
## know how many units of centroid size grows when one unit of 
## LnCS grows

gls2Plot<-gls(sqrt(apodeme)~dist.proc*species*sex,data=full.table,
             weights=varIdent(form=~1|species))

######
## Creating new x-axis values for each species
#####
## This is needed to use the predict() function

proc.dentM<-seq(min(full.table$dist.proc[full.table$species=="denticulata"&full.table$sex=="male"]),
                max(full.table$dist.proc[full.table$species=="denticulata"&full.table$sex=="male"]),length.out=100)

proc.abtaoM<-seq(min(full.table$dist.proc[full.table$species=="abtao"&full.table$sex=="male"]),
                 max(full.table$dist.proc[full.table$species=="abtao"&full.table$sex=="male"]),length.out=100)

proc.longM<-seq(min(full.table$dist.proc[full.table$species=="longirostri"&full.table$sex=="male"]),
                max(full.table$dist.proc[full.table$species=="longirostri"&full.table$sex=="male"]),length.out=100)

proc.dentF<-seq(min(full.table$dist.proc[full.table$species=="denticulata"&full.table$sex=="female"]),
                max(full.table$dist.proc[full.table$species=="denticulata"&full.table$sex=="female"]),length.out=100)

proc.abtaoF<-seq(min(full.table$dist.proc[full.table$species=="abtao"&full.table$sex=="female"]),
                 max(full.table$dist.proc[full.table$species=="abtao"&full.table$sex=="female"]),length.out=100)

proc.longF<-seq(min(full.table$dist.proc[full.table$species=="longirostri"&full.table$sex=="female"]),
                max(full.table$dist.proc[full.table$species=="longirostri"&full.table$sex=="female"]),length.out=100)

## Now the factors

sp.dent<-factor(rep("denticulata",length(proc.dentM)))
sp.abtao<-factor(rep("abtao",length(proc.abtaoM)))
sp.long<-factor(rep("longirostri",length(proc.longM)))
malesNew<-factor(rep("male",length(proc.longM)))
femalesNew<-factor(rep("female",length(proc.longM)))


## Now we can start plotting. If you want to save the graph in tiff
## format in our computer, please remove the # in the tiff() below
## and in dev.off() at the bottom

#tiff(file="apodxProcDist.tiff",units="mm",width=170,height=150,res=600,
#     compression="lzw")
plot(sqrt(apodeme)~dist.proc,pch=c(24,21)[as.numeric(sex)],
     bg=adjustcolor(c("red","gold","gray")[as.numeric(species)],0.7),
     cex=1.5,data=full.table,bty='l',col=NA,las=1,ylab="Apodeme area (square root)",
     xlab="Procrustes' distance")

## Rather awkward, but I prefer to build the pannels in photoshop
## So I plot everything here and move legends around in that

legend("topleft", legend=c(expression(italic("Aegla abtao"),
                                      italic("Aegla denticulata"),
                                      italic("Aegla longirostri")),
                           "Male","Female","M"),
       pch = c(24,24,24,21,21,21),col=NA,
       pt.bg=adjustcolor(c("gray","gold","red","gray","gold","red"),0.7),bty='n',
       cex=1.2,ncol=2)
legend(0.15,3.5,legend="(d)",cex=1.2,bty='n')

## Plotting the lines 

## First we build a data.frame with the new data and then plot the lines. 
## Over the years I found that data.frames give less error messages than 
## using list(). For instance, the 'betareg' package does not work really
## well if you use list(). Thus, data.frames it is. 

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(dist.proc=proc.dentM,species=sp.dent,sex=malesNew)
lines(dentM$dist.proc,predict(gls2Plot,type='response',newdata=dentM),lwd=3,col="yellow4")

##
#AEGLA DENTICULATA - FEMALES
##

dentF<-data.frame(dist.proc=proc.dentF,species=sp.dent,sex=femalesNew)
lines(dentF$dist.proc,predict(gls2Plot,type='response',newdata=dentF),lwd=3,lty=2,col="yellow4")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(dist.proc=proc.abtaoM,species=sp.abtao,sex=malesNew)
lines(abtaoM$dist.proc,predict(gls2Plot,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA ABTAO - FEMALES
##

abtaoF<-data.frame(dist.proc=proc.abtaoF,species=sp.abtao,sex=femalesNew)
lines(abtaoF$dist.proc,predict(gls2Plot,type='response',newdata=abtaoF),col="black",lwd=3,lty=2)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(dist.proc=proc.longM,species=sp.long,sex=malesNew)
lines(longM$dist.proc,predict(gls2Plot,type='response',newdata=longM),col="darkred",lwd=3)

##
#AEGLA LONGIROSTRI - FEMALES
##

longF<-data.frame(dist.proc=proc.longF,species=sp.long,sex=femalesNew)
lines(longF$dist.proc,predict(gls2Plot,type='response',newdata=longF),col="darkred",lwd=3,lty=2)

#dev.off()

#DONE :D
