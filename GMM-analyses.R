rm(list=ls())
#
###
### PACKAGES NEEDED
###
#

library(geomorph)

##LOAD AND REARRANGE DATA

full.table<-read.table("dados-morfo1.csv",h=T,sep=',')

## Making it easier for later on

species<-full.table$species
sex<-full.table$sex

## Reading the TPS file with landmark information in it

aeglas<-readland.tps("aeglas1.TPS",specID="ID")

## Making a sliders file to show what are my semilandmarks

#slid.aeglas<-define.sliders(aeglas[,,1],nsliders=12) 

# Perform this once and run the write.csv command below. By saving
# the sliders file, you can just load the csv file when reanalyzing

#write.csv(slid.aeglas,"sliders.csv",row.names=F) 

slid.aeglas<-read.table("sliders.csv",sep=",",h=T)

#GPA analysis

gpa.aegla<-gpagen(aeglas,curves=slid.aeglas,ProcD=T)

## Now let's put all variables within a geomorph data frame.
## You can skip this step if you want, but this data frame shortens
## the code.

agdt<-geomorph.data.frame(gpa.aegla,species=species,sex=sex)


#####################
# SUMMARIZING SHAPE #
#####################

## PCA analysis on Procrustes coordinates

pca.aegla<-plotTangentSpace(agdt$coords,groups=agdt$species)
write.csv(pca.aegla$pc.shapes,"aegla-pc_shapes.csv") #This will save the consensus
# and other shapes in a .csv file. Useful when plotting deformation grids.

## To know how much each PC axis explains.

pca.aegla$pc.summary

## Now, plotting the data. Note that if you run the plotTangentSpace without
## assigning it to an object, you get a plot right away.

#tiff(file="PCA-with-size.tiff",units="mm",width=170,height=150,res=600,
#     compression="lzw")
plot(pca.aegla$pc.scores[,1],pca.aegla$pc.scores[,2],xlab="PC1 (51.11%)",
     ylab="PC2 (18.4%)",pch=c(24,21)[as.numeric(sex)],
     bg=adjustcolor(c("gray","gold","red")[as.numeric(species)],0.7),
     cex=1.5,data=full.table,bty='l',col=NA,las=1)
#dev.off()

## Since size can influence shape (Klingenberg, 2016), let's account for 
## a common allometric slope. 

#== http://link.springer.com/article/10.1007/s00427-016-0539-2 ==#

anova.aegla<-procD.lm(coords~Csize,data=agdt,
                      RRPP=T,iter=999,logsz=T)
summary(anova.aegla)

## Adding the consensus and adjusting the data frame.

shape.resid<-arrayspecs(anova.aegla$residuals,p=dim(gpa.aegla$coords)[1],
                        k=dim(gpa.aegla$coords)[2])
adj.shape<-shape.resid+array(gpa.aegla$consensus,dim(shape.resid))


## Running the PCA analysis again, but this time with the size-corrected shape.

pca.nosize<-plotTangentSpace(adj.shape,groups=agdt$species,warpgrids=F)
write.csv(pca.nosize$pc.shapes,"aegla-nosize-pc_shapes.csv") #This will save the consensus
# and other shapes in a .csv file. Useful when plotting deformation grids.

## To know how much each PC axis explains.

pca.nosize$pc.summary

## Now, plotting the data. Note that if you run the plotTangentSpace without
## assigning it to an object, you get a plot right away.

#tiff(file="PCA-no-size.tiff",units="mm",width=170,height=150,res=600,
#     compression="lzw")
plot(pca.nosize$pc.scores[,1],pca.nosize$pc.scores[,2],xlab="PC1 (34%)",ylab="PC2 (24.15%)",
     pch=c(24,21)[as.numeric(sex)],
     bg=adjustcolor(c("gray","gold","red")[as.numeric(species)],0.7),
     cex=1.5,data=full.table,bty='l',col=NA,las=1)
#dev.off()

####################################
# DIFFERENTIAL INVESTMENT IN SHAPE #
####################################

## First, let's use the procD.lm() function to get eh full ANOVA table with
## F-values for all our variables. 

model1<-procD.lm(coords~Csize*species*sex,data=agdt,
                 RRPP=T,iter=999,logsz=T)
summary(model1)

## Now, for pairwise comparisons it is a bit more complicated. This analysis tests one model (f1)
## against another more complex model (f2), and does the pairwise comparisons based on that analysis.
## Thus, you have to be very careful on what you add in both models to make sure you are interpreting
## correctly. For more info on this: 
# https://github.com/geomorphR/geomorph/wiki/advanced.procD.lm-for-pairwise-tests-and-model-comparisons #

## In my case, I want to test if sexes and species are investing differently in shape. This means that I am
## interested in the interactions of the factors with the size co-variable (i.e. size * species and 
## size * sex). My first model is thus a model without slopes (i.e. size is not interacting with the factors)
## and my second model is the full model. Although it is not a significant triple interaction (as shown in
## the procD.lm analysis above), it is what I want to test. Morever, I chose to analyse the angles in degrees
## simply because they seem more straightforward to me.

model2<-advanced.procD.lm(f1=coords~log(Csize)+species*sex,f2=~log(Csize)*species*sex,
                          data=agdt,slope=~log(Csize),groups=~species*sex,
                          RRPP=T,iter=999,angle.type = "deg")
summary(model2)

#############################
# TESTING SHAPE VARIABILITY #
#############################

## Here, we will use morphological disparity to test which sex and species has the most
## variable shape. Again, this requires a bit of prior reading to understand what you are 
## really testing. The function morphol.disparity() calculates Procrustes variance based on
## a linear fit, which then procedes onto doing the pairwise comparisons. Therefore, if you
## are interested in analyzing the difference between species, the species factor should not
## go in the formula of the analysis, but rather in the groups= argument. If you add species 
## in the formula you are actually accounting for the difference in species and then testing
## if they are different. 
## In my case, my linear fit is accounting for a common allometric slope, and then doing the
## pairwise comparisons between species and sex. That way, I am ensuring that all the variance 
## left is related to species and sex prior to testing.

mdispar<-morphol.disparity(coords~log(Csize),groups=~species*sex,data=agdt)
mdispar


## The trajectory analysis is simpler. What you are testing needs to go into the formula.
## However, it does not accept continuous co-variables. Therefore, I am running one analysis
## without correcting for size and another one with the size-corrected variables (similar to)
## what I have done in the PCA analyses.
## In this analysis, we are testing how much shape is varying and in which direction. 
## The magnitude of shape change is the length of the line in the plot, and the direction 
## is the angle of the line. For more information, read: 
## Collyer and Adams (2015) - Phenotypic trajectory analysis: comparison of shape change 
## patterns in evolution and ecology. Hystrix, 24: 75-83.

# First, the size uncorrected analysis

traject<-trajectory.analysis(coords~species*sex,data=agdt)
summary(traject)
plot(traject)


# Now, the size corrected analysis

traject2<-trajectory.analysis(shape.resid~species*sex,data=agdt)
summary(traject2)
plot(traject2)

## For publication figures, Michael Collyer suggested messing around their plot function 
## manually, which is what I did. I am not going to post that code, but if someone wants
## more info, just send me and e-mail.

