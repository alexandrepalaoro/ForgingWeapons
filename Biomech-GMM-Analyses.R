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

full.table$lncs.scale<-as.vector(scale(full.table$csize,center=T,scale=T))

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

plot(sqrt(apodeme)~lncs.scale,col=species,data=full.table)

## Okay, one species is clearly different (the red one), the other two
## (the weapons) we will have to test to be sure.
## I won't enter in details about how I got to the gls() function.
## Suffice to say that the fitted were a bit off to use a conventional 
## linear regression. Therefore, I used an analysis that could fix that 
## issue by adding an exponential error structure to my analysis.
## I will run two models: one without the error structure (gls1), and one 
## with the error structure (gls2). Then, I will plot the fitted residuals 
## so you can look whether I did the right thing or not.

gls1<-gls(sqrt(apodeme)~lncs.scale*species,data=full.table)

gls2<-gls(sqrt(apodeme)~lncs.scale*species,data=full.table,
          weights=varExp(form=~lncs.scale))

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
## before continuing 4un the following lines.
##  These introduce methods for 'gls' objects.
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

testInteractions(gls2,pairwise='species',slope='lncs.scale')
testInteractions(gls2,pairwise='species')

# Now, the plot.

plot(sqrt(apodeme)~lncs.scale,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,bty='l',las=1,ylab="Apodeme area (square root)",
     xlab="Centroid size (centered and scaled)")

legend("topleft", legend=c("Display+intense fighting",
                           "No display+no intense fighting",
                           "No display+intense fighting"),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),"white"),bty='n',
       cex=1.2)
# To plot predicted lines, follow the script below.

######
## Creating new x-axis values for each species
#####
## This is needed to use the predict() function

lncs.dentM<-seq(min(full.table$lncs.scale[full.table$species=="denticulata"&full.table$sex=="male"]),
                max(full.table$lncs.scale[full.table$species=="denticulata"&full.table$sex=="male"]),length.out=100)

lncs.abtaoM<-seq(min(full.table$lncs.scale[full.table$species=="abtao"&full.table$sex=="male"]),
                 max(full.table$lncs.scale[full.table$species=="abtao"&full.table$sex=="male"]),length.out=100)

lncs.longM<-seq(min(full.table$lncs.scale[full.table$species=="longirostri"&full.table$sex=="male"]),
                max(full.table$lncs.scale[full.table$species=="longirostri"&full.table$sex=="male"]),length.out=100)


## Now the factors

sp.dent<-factor(rep("denticulata",length(lncs.dentM)))
sp.abtao<-factor(rep("abtao",length(lncs.abtaoM)))
sp.long<-factor(rep("longirostri",length(lncs.longM)))

#Plotting the lines

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(lncs.scale=lncs.dentM,species=sp.dent)
lines(dentM$lncs.scale,predict(gls2,type='response',newdata=dentM),lwd=3,col="darkgray")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs.scale=lncs.abtaoM,species=sp.abtao)
lines(abtaoM$lncs.scale,predict(gls2,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs.scale=lncs.longM,species=sp.long)
lines(longM$lncs.scale,predict(gls2,type='response',newdata=longM),col="black",lwd=3,lty=2)

## Cool. Now moving on to the mechanical advantage analysis.


#####################
## CLAW EFFICIENCY ##
#####################

## Here I used the mechanical advantage (MA) of the claw as a proxy for claw mechanical efficiency.
## MA goes from 0 to 1, so it is a censored data. Luckily, the data range goes from ~0.3 to ~0.7
## which means we can use the gls() function as well - the data never reaches the censored part of
## its distribution. Here, I also changed the error structure for the exponential of the co-variable.

hist(full.table$ma2)

gls3<-gls(ma2~lncs.scale*species,data=full.table,
            weights = varExp(form=~lncs.scale))

RESID3<-resid(gls3,type="normalized")
FIT3<-fitted(gls3)
plot(FIT3,RESID3, xlab="Fitted",ylab="Residuals")

summary(gls3)
anova(gls3)

testInteractions(gls3,pairwise = 'species',slope='lncs.scale')
testInteractions(gls3,pairwise = 'species')

## Now we can make the plot.

plot(ma2~lncs.scale,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,bty='l',las=1,ylab="Mechanical advantage",
     xlab="Centroid size (log)")

legend("topleft", legend=c("Display+intense fighting",
                           "No display+no intense fighting",
                           "No display+intense fighting"),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),"white"),bty='n',
       cex=1.2)

#Now the lines

dentM<-data.frame(lncs.scale=lncs.dentM,species=sp.dent)
lines(dentM$lncs.scale,predict(gls3,type='response',newdata=dentM),lwd=3,col="darkgray")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs.scale=lncs.abtaoM,species=sp.abtao)
lines(abtaoM$lncs.scale,predict(gls3,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs.scale=lncs.longM,species=sp.long)
lines(longM$lncs.scale,predict(gls3,type='response',newdata=longM),col="black",lwd=3,lty=2)


############################################################
## Now, to the plot that will be provided in the paper :D ##
############################################################

tiff(file="Figure3-PANEL.tiff",units="mm",width=300,height=150,res=600,
     compression="lzw")
plot(sqrt(apodeme)~lncs.scale,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,las=1,ylab="Apodeme area (square root)",
     xlab="Centroid size (centred and scaled)")

legend("topleft", legend=c("Display+fierce fighting",
                           "No display+no fierce fighting",
                           "No display+fierce fighting"),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),"white"),bty='n',
       cex=1.2)
legend(1.7,1.1,legend="(a)",cex=1.2,bty='n')

## Plotting the lines 

## First we build a data.frame with the new data and then plot the lines. 
## Over the years I found that data.frames give less error messages than 
## using list(). For instance, the 'betareg' package does not work really
## well if you use list(). Thus, data.frames it is. 

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(lncs.scale=lncs.dentM,species=sp.dent)
lines(dentM$lncs.scale,predict(gls2,type='response',newdata=dentM),lwd=3,col="darkgray")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs.scale=lncs.abtaoM,species=sp.abtao)
lines(abtaoM$lncs.scale,predict(gls2,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs.scale=lncs.longM,species=sp.long)
lines(longM$lncs.scale,predict(gls2,type='response',newdata=longM),col="black",lwd=3,lty=2)

####

plot(ma2~lncs.scale,pch=21,
     bg=c(alpha("grey30",0.75),alpha("black",0.75),"white")[as.numeric(species)],
     cex=1.5,data=full.table,las=1,ylab="Mechanical advantage",
     xlab="Centroid size (centred and scaled)")

legend(1.7,0.56,legend="(b)",cex=1.2,bty='n')

##
#AEGLA DENTICULATA - MALES
##

dentM<-data.frame(lncs.scale=lncs.dentM,species=sp.dent)
lines(dentM$lncs.scale,predict(gls3,type='response',newdata=dentM),lwd=3,col="darkgray")

##
#AEGLA ABTAO - MALES
##

abtaoM<-data.frame(lncs.scale=lncs.abtaoM,species=sp.abtao)
lines(abtaoM$lncs.scale,predict(gls3,type='response',newdata=abtaoM),col="black",lwd=3)

##
#AEGLA LONGIROSTRI - MALES
##

longM<-data.frame(lncs.scale=lncs.longM,species=sp.long)
lines(longM$lncs.scale,predict(gls3,type='response',newdata=longM),col="black",lwd=3,lty=2)


dev.off()

## Now, off to the shape analyses.
# Here, it is good practice to put all data you are going to use in
# a geomorph data frame. You can skip this step if you want, but 
# this data frame shortens the code, and makes analysis run smoother.

agdt<-geomorph.data.frame(gpa.aegla,species=species,
                          icf=full.table$icf)

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

####################################
## ANALYSIS OF SHAPE AND FUNCTION ##
####################################

icf.scale<-as.vector(scale(log(agdt$icf)))

## First, let's use the procD.lm() function to get eh full ANOVA table with
## F-values for all our variables. 

model1<-procD.lm(coords~icf.scale*species,data=agdt,RRPP=T,iter=9999)
summary(model1)

plot(model1,predictor=icf.scale,type=c("regression"),reg.type=c("PredLine"),
     cex=1.3,pch=21,bg=c(alpha("grey30",0.75),alpha("black",0.75),
                         "white")[as.numeric(agdt$species)],
     xlab="Index of closing force (scaled and centered)",bty='l',las=1)

tiff(file="Figure4.tiff",units="mm",width=180,height=280,res=600,
     compression="lzw")

par(mfrow=c(2,1),bty='l')

plot(pca.aegla$pc.scores[,1],pca.aegla$pc.scores[,2],xlab="PC1 (56.88%)",
     ylab="PC2 (17.54%)",pch=21,ylim=c(-0.07,0.07),xlim=c(-0.11,0.11),
     bg=c(alpha("grey30",0.75),alpha("black",0.75),
          "white")[as.numeric(species)],
     cex=1.5,data=full.table,las=1)

legend("topleft",legend="(a)",cex=1.2,bty='n')

plot(model1,predictor=icf.scale,type=c("regression"),reg.type=c("RegScore"),
     cex=1.5,pch=21,bg=c(alpha("grey30",0.75),alpha("black",0.75),
                         "white")[as.numeric(agdt$species)],
     xlab="Index of closing force (centred and scaled)",las=1)

legend("topright",legend="(c)",cex=1.2,bty='n')

legend("bottomleft", 
       legend=c("Display+fierce fighting",
                "No display+no fierce fighting",
                "No display+fierce fighting"),
       pch = c(21,21,21),
       pt.bg=c(alpha("grey30",0.75),alpha("black",0.75),
               "white"),bty='n',cex=1.2)


dev.off()

par(mfrow=c(1,1))


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

model2<-advanced.procD.lm(f1=coords~icf.scale+species,
                          f2=~icf.scale*species,
                          data=agdt,slope=~log(icf),groups=~species,
                          RRPP=T,iter=9999,angle.type = "deg")
summary(model2)

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
# DONE :D
