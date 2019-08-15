# data analysis of Phytophthora infestans population in Peru
library(poppr)
citation("poppr")
vignette("poppr_manual", "poppr")

#import data using genealex format(excel table): make sure the missing data is coded as 0
##make two tables: one with the genotypic data where the first column has the names of genotypes and second column has the population info.
## the second table has all other data.
Pinf_Peru <- read.genalex("Pinf_926E.csv", ploidy=3, geo = FALSE)
Pinf_Peru

#add other information, in this case the rest of the data:
demographics <- read.csv("926E_other3.csv", header=T)
other(Pinf_Peru)$xy <- demographics
other(Pinf_Peru)

#add repeats of the SSR loci
repeats <- c(2,2,2,2,3,2,3,2,2,2,2,2) #number of repeats of loci inthe same order: D13,SSR8,SSR4,Pi04,Pi70,SSR6,Pi63,G11,SSR3_Pi02,SSR11,SSR2,4B
repeats
other(Pinf_Peru)$repeat_lengths <- repeats
other(Pinf_Peru)

#add the dataframe from the other slot into the starta slot
strata(Pinf_Peru) <- other(Pinf_Peru)$xy[-1]
Pinf_Peru

#remove zeros from the dataset (Bruvo dist must be run without 0s)
Pinf_Peru_rc <- recode_polyploids (Pinf_Peru, newploidy = TRUE)

#######################################################################################

#Diversity statistics table by clonal lineage
poppr(Pinf_Peru_rc, sample=999, method=4)

#############################################################################################
#DAPC
library(adegenet)
adegenetTutorial("dapc")
vignette("adegenet-dapc")
#http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
#https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

#Determine the number of clusters in the dataset
K.clust <- find.clusters(Pinf_Peru_rc, clust = NULL, n.pca = NULL, n.clust = NULL,
                         method = "kmeans", stat = "WSS",
                         choose.n.clust = FALSE, criterion = "diffNgroup", 
                         max.n.clust = 10,
                         n.iter = 1e5, n.start = 2, scale = FALSE, truenames = TRUE)

#check the number of clusters(=size)
summary(K.clust)

plot (K.clust$Kstat, ylab="WSS", xlab="K", col="blue")

#look for the appropriate number of n.pca
set.seed(999)
#first identify the range of optimal principal component axes
Pinfx <- xvalDapc(tab(Pinf_Peru_rc, NA.method = "mean"), pop(Pinf_Peru_rc))
#narrow down the number of axes to a range that produced best results and run the analysis
system.time(Pinfx <- xvalDapc(tab(Pinf_Peru_rc, NA.method = "mean"), pop(Pinf_Peru_rc),
                              n.pca = 20:60, n.rep = 1000,
                              parallel = "multicore", ncpus = 4L))

names(Pinfx)
#to view the number of PCs achieving lowest MSE
Pinfx[-1] 

#PCA
dapc.Pinf <- dapc(Pinf_Peru_rc, var.contrib = TRUE, scale = FALSE, n.pca = 42, 
                  n.da = nPop(Pinf_Peru_rc) - 1)
par(mfrow=c(1,1))
scatter(dapc.Pinf, pch = 15:19, cstar = TRUE, mstree = FALSE, legend =TRUE,
        scree.pca = TRUE, scree.da = FALSE, posi.da = "topright",
        posi.pca = "topleft", clab=FALSE, cex.lab=0.5, cleg = 0.5,
        cellipse =2,
        posi.leg = "topright", xax = 1, yax = 2, inset.solid = 1, cex=1)

summary(dapc.Pinf)

scatter(dapc.Pinf, cell=2, pch=15:19, cstar=FALSE, scree.pca=TRUE, posi.pca="topleft",
        scree.da = FALSE,
        axesel=FALSE, legend=TRUE, cleg=1, clab=FALSE, posi.leg = "bottomright")

#######################################################################################
#Rarefaction curve to demonstrate genotypic richness on alternative host species
library("vegan")
#define host species as populations
setPop(Pinf_Peru_rc) <-~hostA
popNames(Pinf_Peru_rc)
#subset the data to include species with more than one sample (and excluding S.tuberosum)
alt_host <- popsub(Pinf_Peru_rc, 
                   sublist=c("S. zahlbruckneri","S. ochrantum",
                             "S. lycopersicum","S. peruvianum", 
                             "S. habrochaites", "S. caripense", "S. billhoockerii",
                             "S. chrysotrichum","S. medians",
                             "S. wittmackii","S. hypacrarthrum","S. pennellii","S. sogarandinum",
                             "S. bulbocastanum X S.tuberosum",
                             "S. muricatum","Iochroma grandiflorum","S. grandidentatum","S. huancabambense"))  
#number of MLG by host specie
BB <- mlg.table(alt_host)
min_sample <- min(rowSums(BB))
#draw the rarefaction plot using vegan package
rarecurve2(BB, step =1, sample = min_sample, col="blue", label = T,
           xlab = "Sample Size", ylab = "Number of MLGs",
           selected.labels = c(rep(T,8), rep(F,4),T,F,T,F,F,T))

#Rarefaction curve to demonstrate genotypic richness at the field level
#define potato fields as populations
setPop(Pinf_Peru_rc) <-~field_number
popNames(Pinf_Peru_rc)
#number of MLG in each field
BB <- mlg.table(Pinf_Peru_rc)
min_sample <- min(rowSums(Pinf_Peru_rc))
rarecurve3(BB, step =1, sample = min_sample, col="blue", cex=0.6, label = T,
           xlab = "Sample Size", ylab = "Number of MLGs")

#######################################################################################
##Minimum Spanning Network (MSN)
library(RColorBrewer)
#create RColorBrewer palette
twelve_colors <-brewer.pal(n=12, name="Paired")

#MSN for EC1 lineage
#group the data by old and current samples
setPop(EC1) <-~time
pcmsn <- bruvo.msn(EC1, replen=repeats, palette=twelve_colors)
plot_poppr_msn(EC1, pcmsn, size.leg = FALSE, inds=0, beforecut = TRUE, cutoff = 0.076)


######################################################################################
#Analysis of Molecular Variance, AMOVA 

#Hierarchical AMOVA by clonal lineage (EC1, PE3, PE7, US1), region (North, Center, South) and 
#political department (Amazonas, Cajamarca, La Libertad, Piura, Ancash, 
#Huancavlica, Huanuco, Junin, Lima, Pasco, Apurimac, Arequipa, Cusco, Puno, Tacna)
amova.result <-poppr.amova(Pinf_Peru_rc, ~lineage/Region/Department, clonecorrect =TRUE,
                           method="ade4", within=FALSE)


#Hierarchical AMOVA by geographical location (potato field) and host specie 
#among the isolates of the current collection (collected 2016-17)
amova.result <-poppr.amova(current, ~field_number/hostA, clonecorrect =TRUE,
                           method="ade4", within=FALSE)  
#test for significanse
amova.result
amova.test <- randtest(amova.result, nrepet=999)
plot(amova.test)
amova.test

#########################################################################################

