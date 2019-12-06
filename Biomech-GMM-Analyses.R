rm(list=ls())
#
###
### PACKAGES NEEDED
###
#

library(geomorph)
library(nlme)
library(phia)
library(betareg)
library(lmtest)
library(scales)

##LOAD DATA

full.table<-read.table("males-morpho-data.csv",h=T,sep=',')

## This is the main factor I am interested in, so I will making it easier to use while coding 
## I usually don't do this, but some functions behave weirdly if they are not attached to R
## (at least on my pc).

species<-full.table$species

# First, I am doing the biomechanical analyses and then I will proceed with the shape analyses. 
# However, before I do that, I need to calculate my proxy for claw size - the centroid size. 
# To do so, I need to perform the GPA analysis.
# Thus, this first section is a preamble to my biomechanical analyses.

###############################
## PERFORMING GPA ANALYSIS   ##
###############################

## Reading the TPS file with the landmarks already digitized

aeglas<-readland.tps("machos_aegla.TPS",specID="ID")

## Making a sliders file to show where the semilandmarks are located

# slid.aeglas<-define.sliders(aeglas[,,1],nsliders=11) 
# the information regarding which are semi-landmarks are in the paper and it will be uploaded
# alongside the TPS and xls file

# Perform this once and run the write.csv command below. By saving
# the sliders file, you can just load the csv file when reanalyzing

# write.csv(slid.aeglas,"sliders.csv",row.names=F) 

slid.aeglas<-read.table("sliders.csv",sep=",",h=T)

###### GPA analysis ########
 
gpa.aegla<-gpagen(aeglas,curves=slid.aeglas,ProcD=T)

plot(gpa.aegla) #just checking fi everything is OK

full.table$csize<-log(gpa.aegla$Csize) #saving csize on the main data.frame 

########################################################################
## Scaling the proxy for claw size that will be used in the analyses  ##
########################################################################

head(full.table)

abtao<-full.table[full.table$species=='abtao',]
dent<-full.table[full.table$species=='denticulata',]
long<-full.table[full.table$species=='longirostri',]

abtao$lncs.center<-as.vector(scale(abtao$csize,center=T,scale=F))
dent$lncs.center<-as.vector(scale(dent$csize,center=T,scale=F))
long$lncs.center<-as.vector(scale(long$csize,center=T,scale=F))

abtao$ma.center<-as.vector(scale(abtao$ma2,center=T,scale=F))
dent$ma.center<-as.vector(scale(dent$ma2,center=T,scale=F))
long$ma.center<-as.vector(scale(long$ma2,center=T,scale=F))

abtao$mat.center<-as.vector(scale(abtao$maT,center=T,scale=F))
dent$mat.center<-as.vector(scale(dent$maT,center=T,scale=F))
long$mat.center<-as.vector(scale(long$maT,center=T,scale=F))

abtao$apod.center<-as.vector(scale(sqrt(abtao$apodeme),center=T,scale=F))
dent$apod.center<-as.vector(scale(sqrt(dent$apodeme),center=T,scale=F))
long$apod.center<-as.vector(scale(sqrt(long$apodeme),center=T,scale=F))

abtao$icf.center<-as.vector(scale(abtao$icf,center=T,scale=F))
dent$icf.center<-as.vector(scale(dent$icf,center=T,scale=F))
long$icf.center<-as.vector(scale(long$icf,center=T,scale=F))


head(abtao)

full.table$lncs.center<-c(abtao$lncs.center,dent$lncs.center,long$lncs.center)
full.table$ma.center<-c(abtao$ma.center,dent$ma.center,long$ma.center)
full.table$mat.center<-c(abtao$mat.center,dent$mat.center,long$mat.center)  
full.table$apod.center<-c(abtao$apod.center,dent$apod.center,long$apod.center) 
full.table$icf.center<-c(abtao$icf.center,dent$icf.center,long$icf.center) 

head(full.table)


# Now we can move to the biomechanical analyses

############################
###### CLAW STRENGTH #######
############################

## I will square-root the apodeme because it is an area measurement, while
## centroid size is a linear measurement. Thus, any significant correlation
## I find between these two variables could be entirely spurious due to
## the differences in how they increase. It has nothing to do with why we
## usually transform our data. For more information on this, please see:

#-- Pélabon C, et al. (2014). Evolution of morphological allometry. --#
#-- Annals of the New York Academy of Sciences, 1320(1), 58-75.     --#

hist(sqrt(full.table$apodeme))

## Well, this seems a bit bimodal. I will probably have to do something
## about the error structure of this analysis

## Just so we have an idea with what we are dealing with here

plot(sqrt(apodeme)~lncs.center,col=species,data=full.table)

## Okay, one species is clearly different (the red one), the other two
## (the weapons) we will have to test to be sure.
## I won't enter in details about how I got to the gls() function.
## Suffice to say that the fitted values were a bit off to use a conventional 
## linear regression. Therefore, I used an analysis that could fix that 
## issue by adding an exponential error structure to my analysis.
## I will run two models: one without the error structure (gls1), and one 
## with the error structure (gls2). Then, I will plot the fitted residuals 
## so you can look whether I did the right thing or not.

gls1<-gls(sqrt(apodeme)~lncs.center*species,data=full.table)

gls2<-gls(sqrt(apodeme)~lncs.center*species,data=full.table,
          weights=varExp(form=~lncs.center))

## Plotting residuals against fitted values

par(mfrow=c(1,2))

RESID1<-resid(gls1,type="normalized")
FIT1<-fitted(gls1)
plot(FIT1,RESID1, xlab="Fitted",ylab="Residuals",
     main="With no weights")

RESID2<-resid(gls2,type="normalized")
FIT2<-fitted(gls2)
plot(FIT2,RESID2, xlab="Fitted",ylab="Residuals",
     main="With weights")

par(mfrow=c(1,1))

## Residuals seem well dispersed when using the error structure.
## Thus, I will proceed with that analysis.

summary(gls2)
anova(gls2)

## Cool, so there is a significant interaction between size and species.
## To provide a pairwise test, I will use the function testInteractions in
## the 'phia' package. This function is able to calculate the pairwise 
## differences while accounting for the interaction term. However, 
## before continuing run the following lines.
##  These introduce an extension for 'gls' objects.

model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}

# The slope argument gives you the pairwise differences of the interaction.
# If you remove that argument, it will provide pairwise differences between
# the main factor 'species'.

testInteractions(gls2,pairwise='species',slope='lncs.center')
testInteractions(gls2,pairwise='species')

testInteractions(teste1,pairwise='species',slope='log(ceph.length)')

# Now, the plot.

plot(log(apodeme)~log(ceph.length),pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,bty='l',las=1,ylab="Apodeme area (cm)",
     xlab="Centroid size (centred)")

legend("topleft", legend=c("Aegla abtao",
                           "Aegla denticulata",
                           "Aegla longirostri"),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),"white"),bty='n',
       cex=1.2)

# To plot predicted lines, follow the script below.

######
## Creating new x-axis values for each species
#####
## This is needed to use the predict() function

lncs.dentM<-seq(min(full.table$lncs.center[full.table$species=="denticulata"&full.table$sex=="male"]),
                max(full.table$lncs.center[full.table$species=="denticulata"&full.table$sex=="male"]),length.out=100)

lncs.abtaoM<-seq(min(full.table$lncs.center[full.table$species=="abtao"&full.table$sex=="male"]),
                 max(full.table$lncs.center[full.table$species=="abtao"&full.table$sex=="male"]),length.out=100)

lncs.longM<-seq(min(full.table$lncs.center[full.table$species=="longirostri"&full.table$sex=="male"]),
                max(full.table$lncs.center[full.table$species=="longirostri"&full.table$sex=="male"]),length.out=100)


## Now the factors

sp.dent<-factor(rep("denticulata",length(lncs.dentM)))
sp.abtao<-factor(rep("abtao",length(lncs.abtaoM)))
sp.long<-factor(rep("longirostri",length(lncs.longM)))

#Plotting the lines

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(lncs.center=lncs.dentM,species=sp.dent)
lines(dentM$lncs.center,predict(gls2,type='response',newdata=dentM),lwd=3,col="darkgray")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs.center=lncs.abtaoM,species=sp.abtao)
lines(abtaoM$lncs.center,predict(gls2,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs.center=lncs.longM,species=sp.long)
lines(longM$lncs.center,predict(gls2,type='response',newdata=longM),col="black",lwd=3,lty=2)

## Cool. Now moving on to the mechanical advantage analysis.


#####################
## CLAW EFFICIENCY ##
#####################

## Here I used the mechanical advantage (MA) of the claw as a proxy for claw mechanical efficiency.
## MA goes from 0 to 1, so it is a censored data. Luckily, the data range goes from ~0.3 to ~0.7
## which means we can use the gls() function as well - the data never reaches the censored part of
## its distribution. Here, I also changed the error structure for the exponential of the co-variable.

hist(full.table$maT)

b1<-betareg(maT~lncs.center*species,data=full.table)
b2<-betareg(maT~lncs.center+species,data=full.table)
b3<-betareg(maT~lncs.center,data=full.table)
b4<-betareg(maT~species,data=full.table)
b5<-betareg(maT~1,data=full.table)

lrtest(b1,b2)
lrtest(b4,b5)
lrtest(b3,b5)


summary(b1)
summary(b2)

## Now we can plot.

plot(maT~lncs.center,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,bty='l',las=1,ylab="Mechanical advantage",
     xlab="Centroid size (centred)")

lines(dentM$lncs.center,predict(b1,type='response',newdata=dentM),lwd=3,col="darkgray")
lines(abtaoM$lncs.center,predict(b1,type='response',newdata=abtaoM),col="black",lwd=3)
lines(longM$lncs.center,predict(b1,type='response',newdata=longM),col="black",lwd=3,lty=2)


plot(resid.ma~lncs.center,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,bty='l',las=1,ylab="Residuals of Mechanical advantage on apodeme",
     xlab="Centroid size (centred)")

plot(resid(m1)~lncs.center,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,bty='l',las=1,ylab="Mechanical advantage",
     xlab="Centroid size (centred)")

legend("topleft", legend=c("Aegla abtao",
                           "Aegla denticulata",
                           "Aegla longirostri"),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),"white"),bty='n',
       cex=1.2)

#Now the lines

dentM<-data.frame(lncs.center=lncs.dentM,species=sp.dent)
lines(dentM$lncs.center,predict(gls3,type='response',newdata=dentM),lwd=3,col="darkgray")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs.center=lncs.abtaoM,species=sp.abtao)
lines(abtaoM$lncs.center,predict(gls3,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs.center=lncs.longM,species=sp.long)
lines(longM$lncs.center,predict(gls3,type='response',newdata=longM),col="black",lwd=3,lty=2)



############################################################
## Now, to the plot that will be provided in the paper :D ##
############################################################
pdf("FIGURE1.pdf",width=10,height=5)
#png(file="FIGURE1.png",units="mm",width=250,height=130,res=600)
#tiff(file="FIGURE1.tiff",units="mm",width=450,height=130,res=600,
#     compression="lzw")
par(mfrow=c(1,2),bty='l',mar=c(5,5,4,2))
plot(sqrt(apodeme)~lncs.center,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=2,data=full.table,las=1,ylab=expression(sqrt("Apodeme area (cm)")),
     xlab="Centroid size (centred)")

legend("topleft", legend=c(expression(italic("Aegla abtao")),
                           expression(italic("Aegla denticulata")),
                           expression(italic("Aegla longirostri"))),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),"white"),bty='n')
legend(0.5,1.2,legend="(a)",bty='n')

## Plotting the lines 

## First we build a data.frame with the new data and then plot the lines. 
## Over the years I found that data.frames give less error messages than 
## using list(). For instance, the 'betareg' package does not work really
## well if you use list(). Thus, data.frames it is. 

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(lncs.center=lncs.dentM,species=sp.dent)
lines(dentM$lncs.center,predict(gls2,type='response',newdata=dentM),lwd=3,col="darkgray")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs.center=lncs.abtaoM,species=sp.abtao)
lines(abtaoM$lncs.center,predict(gls2,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs.center=lncs.longM,species=sp.long)
lines(longM$lncs.center,predict(gls2,type='response',newdata=longM),col="black",lwd=3,lty=2)

####


plot(maT~lncs.center,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,bty='l',las=1,ylab="Mechanical advantage",
     xlab="Centroid size (centred)")

lines(dentM$lncs.center,predict(b1,type='response',newdata=dentM),lwd=3,col="darkgray")
lines(abtaoM$lncs.center,predict(b1,type='response',newdata=abtaoM),col="black",lwd=3)
lines(longM$lncs.center,predict(b1,type='response',newdata=longM),col="black",lwd=3,lty=2)

legend(0.5,1,legend="(b)",bty='n')

dev.off()

#plot(maT~lncs.center,pch=21,cex.lab=1.2,cex.axis=1.2,
#     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
#     cex=3,data=full.table,las=1,ylab="Mechanical advantage",
#    xlab="Centroid size (centred)")

#legend(0.6,1,legend="(b)",cex=1.5,bty='n')

##
#AEGLA DENTICULATA - MALES
##

#dentM<-data.frame(lncs.center=lncs.dentM,species=sp.dent)
#lines(dentM$lncs.center,predict(gls3,type='response',newdata=dentM),lwd=5,col="grey")

##
#AEGLA ABTAO - MALES
##

#abtaoM<-data.frame(lncs.center=lncs.abtaoM,species=sp.abtao)
#lines(abtaoM$lncs.center,predict(gls3,type='response',newdata=abtaoM),col="black",lwd=5)

##
#AEGLA LONGIROSTRI - MALES
##

#longM<-data.frame(lncs.center=lncs.longM,species=sp.long)
#lines(longM$lncs.center,predict(gls3,type='response',newdata=longM),col="black",lwd=5,lty=2)

## apodeme vs. mechanical advantage ##


##### Figure to show the correlation between MA and muscle size

png("FigS1.png",res=600,width=15,height=12,units='cm')
plot(maT~apod.center,data=full.table,pch=21,bty='l',las=1,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=2,las=1,ylab="Mechanical advantage",
     xlab="Apodeme area (centred)")

legend("bottomright", legend=c(expression(italic("Aegla abtao")),
                           expression(italic("Aegla denticulata")),
                           expression(italic("Aegla longirostri"))),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),"white"),bty='n')
dev.off()


## Now, off to the shape analyses.
# Here, it is good practice to put all data you are going to use in
# a geomorph data frame. You can skip this step if you want, but 
# this data frame shortens the code, and makes analysis run smoother.

agdt<-geomorph.data.frame(gpa.aegla,species=species,
                          icf=full.table$icf.center)


#####################
# SUMMARIZING SHAPE #
#####################

## PCA analysis on Procrustes coordinates

pca.aegla<-plotTangentSpace(agdt$coords,groups=agdt$species)

# write.csv(pca.aegla$pc.shapes,"aegla-pc_shapes.csv") 
# This will save the consensus and other shapes in a .csv file. Useful when 
# plotting deformation grids.

## To know how much each PC axis explains.


pca.aegla$pc.summary

tiff(file="FIGURE2.tiff",units="mm",width=160,height=140,res=600,
     compression="lzw")

plot(pca.aegla$pc.scores[,1],pca.aegla$pc.scores[,2],xlab="PC1 (56.88%)",
     ylab="PC2 (17.54%)",pch=21,ylim=c(-0.07,0.07),xlim=c(-0.11,0.11),
     bg=c(alpha("grey30",0.75),alpha("black",0.75),
          "white")[as.numeric(species)],
     cex=2,data=full.table,las=1,bty='l')

legend("topright", 
       legend=c(expression(italic("Aegla abtao")),
                expression(italic("Aegla denticulata")),
                expression(italic("Aegla longirostri"))),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),
               "white"),bty='n',cex=1.2)
text(-0.1,0.06,"(a)")

dev.off()

par(mfrow=c(1,1))

####################################
## ANALYSIS OF SHAPE AND FUNCTION ##
####################################

## First, let's use the procD.lm() function to get eh full ANOVA table with
## F-values for all our variables. 

model1<-procD.lm(coords~icf*species,data=agdt,RRPP=T,iter=9999)
summary(model1)

## Now, for pairwise comparisons it is a bit more complicated. This analysis tests one model (f1)
## against another more complex model (f2), and does the pairwise comparisons based on that analysis.
## Thus, you have to be very careful on what you add in both models to make sure you are interpreting
## correctly. For more info on this: 
# https://github.com/geomorphR/geomorph/wiki/advanced.procD.lm-for-pairwise-tests-and-model-comparisons #

## In my case, I want to test if different functions generate different amounts of shape changes.  
## This means that I am interested in the interaction between shape and claw efficiency (i.e. size * species and). 
## My first model is thus a model without slopes (i.e. size is not interacting with the factor)
## and my second model is the full model. Morever, I chose to analyse the angles in degrees
## simply because they seem more intuitive to me.

## Two sets of p-values will be provided. The first set will provide the magnitude of the shape change 
## (i.e. how much shape changes by each unit of claw efficiency). The second set provides the direction
## of shape change.

model2<-pairwise(fit=model1,groups=agdt$species,covariate =agdt$icf)

######################################
## TESTING FUNCTIONAL AMPLIFICATION ##
######################################

## Here, we will use morphological disparity to test for functional amplification. 
## Functional amplification occurs when several morphological characters correlate
## and/or coevolve to enhance a desired function. This means that, in my case, shape
## and function should have a high correlation in the fighting weapon, but not in the others.
## The function morphol.disparity() calculates Procrustes variance, which is a correlate of 
## calculating the residuals of a linear regression. This function receives a formula, or an 
## object from procD.lm() (or PGLS) and calculates the residuals from that fit. 
## Then, it performs pairwise comparisons between the residuals using the levels of the factor
## you want as the groups (argument 'groups'). In my case, my linear fit the full model. 
## I want to test whether the residuals of the fighting weapon are larger than the other functions
## Thus, I am accounting for the different allometric slopes and then performing the
## pairwise comparisons between the residuals of the species. 

mdisp<-morphol.disparity(f1 = model1,groups = ~species,data=agdt,iter=9999)
mdisp

## This is it for the main analyses.

###------
## Now, for the standard log-log allometry analyses that can be seen
## in several papers I will just do a bunch of log-log regressions.
## The data are freely available, so you can play around if you want.

### Allometry calculations----
head(full.table)

alom1.abt<-lm(log(apodeme)~log(ceph.length),data=abtao)
plot(alom1.abt,which=1)
summary(alom1.abt)
confint(alom1.abt)

alom1.dent<-lm(log(apodeme)~log(ceph.length),data=dent)
plot(alom1.dent,which=1)
summary(alom1.dent)
confint(alom1.dent)

alom1.long<-lm(log(apodeme)~log(ceph.length),data=long)
plot(alom1.long,which=1)
summary(alom1.long)
confint(alom1.long)

alom2.abt<-lm(log(maT)~log(ceph.length),data=abtao)
plot(alom2.abt,which=1)
summary(alom2.abt)
confint(alom2.abt)

alom2.dent<-lm(log(maT)~log(ceph.length),data=dent)
plot(alom2.dent,which=1)
summary(alom2.dent)
confint(alom2.dent)

alom2.long<-lm(log(maT)~log(ceph.length),data=long)
plot(alom2.long,which=1)
summary(alom2.long)
confint(alom2.long)

alom3.abt<-lm(log(claw.length)~log(ceph.length),data=abtao)
plot(alom3.abt,which=1)
summary(alom3.abt)
confint(alom3.abt)

alom3.dent<-lm(log(claw.length)~log(ceph.length),data=dent)
plot(alom3.dent,which=1)
summary(alom3.dent)
confint(alom3.dent)

alom3.long<-lm(log(claw.length)~log(ceph.length),data=long)
plot(alom3.long,which=1)
summary(alom3.long)
confint(alom3.long)
  