# The purpose of this script is to generate Figure 3 from the manuscript Orkney et al., 2026
# 'Declining rates of evolution and limited convergence in bat
# sternum shape can be explained by a ratchet of specialisation'
# The analysis will determine whether roosting in enclosed spaces represents a constraint on bat thoracic skeleton
# evolution. We hypothesise that additional functional constraints upon sternum evolution introduced by enclosed roosting behaviours
# have the effect of reducing the available degrees of freedom for neutral evolution, overall restricting the evolutionary variance
# of thoracic features.  
# The available data are sternal shape landmark constellations across a diversity of bat species.
# An estimate of the shared ancestry between species is available from Shi & Rabosky 2015 (https://doi.org/10.1111/evo.12681)
# A species-wise binary bag of words is available from Orkney et al., 2025 (https://doi.org/10.1038/s41559-024-02572-9) 
# describing roosting ecology, which can be collapsed into a dichotomous
# 'enclosed/exposed' categorisation. (We use the language 'compressed/free' interchangeably in comments here.) 
# This script conducts the core analyses and produces Figure 3 of the main manuscript Orkney et al., 2026. 

#
# The analyses consist of 
# 
# 1) A Brown-Forsythe Test for different evolutionary variances in Sternum shape between bat lineages that roost in exposed/enclosed spaces
#
# 2) A stochastic character map set produced using an all-rates-different Markov Model of exposed/enclosed roosting ecology
#
# 3) A set of model estimates to determine whether Brownian Motion, Early Burst or Roosting-contingent evolutionary variance represent 
# likelier parameterizations of our data. 
#
# 4) A Chi-squared test to determine whether the ancestral estimate of roosting ecology as exposed/enclosed influences model fit preferences
#
# 5) Visualizations. 
#
# We predicted that the Brown-Forsythe test would reveal that enclosed roosting is associated with reduced variance in sternal shape evolution
# because we conceptualized transitions to roosting in enclosed environments as a second-order adaptive demand placed on bat thoracic osteology
# with the potential to limit the range of sternal shapes that satisfy adaptive demands. 
#
# We predicted that transitions to enclosed roosting should be more likely than transitions to exposed roosting. 
#
# We predicted that sternal shape variety in bats will tend to be explained as an early burst of diversification. We predicted
# That roosting ecology should explain this structure, with exposed roosting associated with early branches of the tree with high rates of change.
# We anticipate that the accumulation of enclosed roosting ecologies through time imposes a constraint on bat sternal evolutionary variance. 
#
# We predicted that reconstructions of exposed/enclosed roosting ecology in which bats are estimated to have ancestrally roosted in exposed
# environments would better explain decaying rates of bat sternal evolutionary variance through time, being more commensurate with existing knowledge
# of the bat fossil record. 

# This script was written by Dr. A Orkney and the final version was compiled on July 19th 2026. 


setwd()
# Set the directory to the location of the landmark data
# This will change if you download the data to a personal computer

load('sternum_array_sep_27_2024.RData')
# load the data

taxa <- dimnames(sternum.array)[[3]]
# This is a vector of available bat species. 

metadata <- read.csv('Bat_CT_process_list_Andrew_only.csv')
families <- metadata$Family[ match(taxa,metadata$Shi_match) ]
names(families) <- taxa
original.names <- paste(metadata$Genus[match(taxa,metadata$Shi_match)], metadata$Species[match(taxa,metadata$Shi_match)],sep='_')
names(original.names)<-taxa
# Substitute names in our collected taxa with congeners on the phylogeny of Shi & Rabosky 2015, if 
# a direct match is not available. 

library(ape) # v 5.8-1
# Package for managing family trees

library(phytools) # v 2.4-4
# Package for managing family trees

bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# This phylogeny far exceeds the number of taxa for which we have landmark constellations
# we must therefore prune the phylogeny 

pruned.tree <- keep.tip(bat.tree,dimnames(sternum.array)[[3]])
# prune the bat tree to the taxa of interest

pruned.tree <- drop.tip(pruned.tree, c('Micropteropus_pusillus','Molossops_temminckii'))
# I don't have ecological metadata for Micropteropus pussilus, and the sternum measurements for Molossops temminckii were erroneous
# Therefore prune these taxa

pruned.tree <- phytools::force.ultrametric(pruned.tree)
# Pruned tree constructed 

# 70 bat species. Thank you Elizabeth Augustin, thank you Beyonca Akers for your hard work 

library(geomorph) # v 4.0.10 
# Package for shape data management 

# We need to align the landmark constellations into a common reference frame:
sliders<- rbind( geomorph::define.sliders(c(7, 17:20,9), write.file=F) ,
                 geomorph::define.sliders(c(8, 21:24,10), write.file=F) ,
                 geomorph::define.sliders(c(3, 25:31,4), write.file=F) ,
                 geomorph::define.sliders(c(1, 32:39,2), write.file=F) ,
                 geomorph::define.sliders(c(15, 40:48,13), write.file=F) ,
                 geomorph::define.sliders(c(16, 49:57,14), write.file=F) ,
                 geomorph::define.sliders(c(15, 58:67,12), write.file=F))
# These landmark curves will be allowed to slide along tangent vectors

GPA.sternum<-geomorph::gpagen(sternum.array, curves=sliders, approxBE=T )
# Align the constellations of landmarks into a common reference frame. 

coords <- GPA.sternum$coords[,,pruned.tree$tip]
# These are the aligned coordinates. 

x<- geomorph::two.d.array(coords)
# It will be convenient to treat the coordinates as a 2-D array.

# https://blog.phytools.org/2012/03/testing-whether-pagels-lambda-is.html
lambda <- phytools::phylosig(x,tree=pruned.tree,method='lambda', test=T, start=1)
# Estimate the sorting of trait values by shared ancestry using Pagel's Lambda 

C <- ape::vcv.phylo(pruned.tree)
C[which(upper.tri(C)==T | lower.tri(C)==T)] <- C[which(upper.tri(C)==T | lower.tri(C)==T)]*lambda$lambda
# This is the matrix of shared ancestry among bat species.
# The off-diagonals are scaled by the estimate of Pagel's lambda- an index of the autoregressive quality of traits
# upon the family tree that relates species. 
# Lambda is below one, so this reflects a contribution of a white-noise process that disturbs tip values. 

# We are now going to investigate whether the vectors between the ancestral condition, and observed tip-state, 
# embedded in a Euclidean geometry sympathetic with sample non-independence, demonstrates that
# the variance of path lengths in free-roosting bats is greater than in tight-roosting bats. 
#
# This statement is equivalent to an ivnestigation of the variance of evolutionary change associated with each 
# species of bat. Are some lineages more variabel than others and why might that be? 

# The linear algebraic reasoning is as follows: 
# Path lengths = sqrt(sum((UiSi)^2)) 
# where i is defined such that Si/sum(diag(S)) > 0.05
# and USV' = Z
# We decompose the phylogenetically-sphered trait measurements by Singular Value Decomposition 
# and Z = (Qsqrt(L)Q')^-1(x-a)
# Z is defined as the Sphered differences in observations from a phylogenetic weighted estimate of the mean
# and QLQ' = C
# Eigendecomposition of C, the matrix of pairwise shared ancestry between species 
# and a = 1'Cx/sum(C^-1)

ones <- matrix(1,dim(x)[1],1)
a <- (t(ones)%*%solve(C)%*%x)/sum(solve(C))
# This is the vectorized estimated ancestral state; 'phylogenetic mean'. 

eigen_C <- eigen(C) 
Q <- eigen_C$vectors
L <- eigen_C$values
# Eigen decomposition of C; the basis vectors of shared ancestry

x<-x[rownames(C),]
Z <- solve(Q%*%diag(sqrt(L))%*%t(Q))%*%(x-ones%*%a)
# Phylogenetic transform performed on x

USV <- svd(Z)
# Singular value decomposition of the rectangular matrix xc
i <- max(which((USV$d/sum(USV$d))>0.05))
# Which eigenvectors explain more than 5% of variance?
U <- USV$u[,1:i]
S <- diag(USV$d[1:i])
# Distribution of species over their eigenvectors, and the explained variance

magxc <- sqrt( rowSums( (U%*%S)^2 ) )
names(magxc)<-rownames(C)
# The magnitude of the vectors from the ancestral condition to the observed species states, 
# in a basis that is sympathetic with the evolutionary history of bats, omitting noisy processes.

# We need to load roosting style categorizations 
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney')
metadata <- read.csv('Bat_eco_metadata.csv')
# Load ecological metadata that documents roosting ecology across a sample of bat species. 

eco <- metadata

binary.roost.styles <- eco[,c(6:16)]
# We actually have 112 taxa we can use for imputation of exposed roosting and enclosed roosting evolutionary history
# even though we only have 70 valid sterna. We're going to use as many as we can to reconsuct the history.

binary.roost.styles<-as.matrix(binary.roost.styles)
# Prepare the data as a matrix to make it easy to index. 

binary.roost.styles[which(binary.roost.styles=='?' | binary.roost.styles=='' )] <- 0 
binary.roost.styles <- apply(binary.roost.styles,2,FUN=as.numeric)
# Ensure the data class is numeric

binary.roost.styles[which(eco$Shi=='Myotis_grisescens'),9]<-1
# Myotis grisescens is observed in photographs to hide in crevices. A manual correction is made


compression <- rowSums(binary.roost.styles[,c(1,2,3,4,7,9)])
compression[which(compression>1)]<-1
# Roosting inside termite mounds, animal burrows, inside trees, 
# and crevices etc is considered 'enclosed' 
# I do not consider caves or artificial roosts to represent enclosed environments unless they are
# accompanied by the word 'crevice'.
# This is because many Pteropodid species enjoy cave entrances but do not actually venture
# inside their dark recesses, 
# Most bats roost in enclosed environments.

names(compression)<- eco$Shi
# Name the vector

# A Brown-Forsythe version of Levene's test will be performed to determine whether roosting
# in enclosed environments is associated with reduced variability in the magnitude of sternal evolution across
# bat history. 

# car package dependency; car v 3.1-5
car::leveneTest(lm(magxc~factor(compression)[rownames(x)]),center='median')
# The magnitude of vectors of sternum shape evolution from the ancestral state are significantly more variable
# in the bat species that roost in exposed environments, compared to bats that roost in enclosed spaces. 

mad(magxc[names(which(compression[rownames(x)]=='1'))])
mad(magxc[names(which(compression[rownames(x)]=='0'))])
# Calculations of the median average differences of both groups
# Observe a value in enclosed-roosting bats that is less than half the value observed in exposed-roosting species. 


# An attempt will now be made to determine whether estimating the history of roosting ecology evolution across bats
# substantially improves explained variance in sternal shape, contrasted against Brownian motion and Early Burst dynamics. 
# In particular, out-performing an Early-burst model would imply that roosting ecology is a substantial determinant
# of variation in sternum evolutionary rates through bat evolutionary history. 

roost<-setNames(compression,names(compression))
roost[which(roost==1)]<-'compressed'
roost[which(roost==0)]<-'free'
roost<-factor(roost)
# Transform roosting into a factor variable.

library(phytools) # v 2.4-4
# Package for managing phylogenetic trees.

eco.tree <- phytools::force.ultrametric(ape::keep.tip(bat.tree,names(roost)))
# It is necessary that the tree is ultrametric for the next steps of analysis. 

roost <- roost[eco.tree$tip]
# Pruned the roosting factor variable to the taxa of interest.

ard_model<-fitMk(eco.tree,roost,model="ARD",pi='equal')
# Fit a Markov Model of roosting ecology in which transitions between enclosed and exposed roosting can occur asymmetrically 
AIC(ard_model)
# The prior on roosting ecology, given our data, cannot be constrained.
# the root could conceivably include roosting in compressed environments. 
# Fossil limb proportions of bats are quite long and resemble arboreal clambering mammals 
# according to qualitative comments by Nancy Simmons in 2008, other work etc. See our published manuscript. 
# Hence, more confidence should be invested in the stochastic character maps which resolve
# free-roosting root states. 

set.seed(1)
# For reproducibility
mcmc100 <- simmap(ard_model, nsim=100, method='mcmc')
# Markov-Chain Monte-Carlo estimate of possible histories of roosting ecology across a sample of bats, 
# which encompasses those for which sternum data is available as well as their outgroups. 
# This is, in essence, a variety of guesses about which branches occurred under which roosting ecology. 
 
obj <- summary(mcmc100)
# Create a summary of the 100 stochastic character maps 
crop.mcmc100 <- list()
for(i in 1:length(mcmc100)){
	crop.mcmc100[[i]] <- keep.tip.simmap(phy=mcmc100[[i]],rownames(x)) 
}
class(crop.mcmc100)<-class(mcmc100)
crop.obj <- summary(crop.mcmc100)
# Cropping to taxa for which sternum data are available 

crop.obj$ace[1,]
# The root node is 50% compressed and 50% free for this seed
# This reflects the prior choice

node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Pteropodidae'))),rownames(crop.obj$ace))
crop.obj$ace[node,]
# 94% probability of free-roosting ancestor for Pteropodidae

node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Hipposideridae' | families=='Rhinolophidae'))),rownames(crop.obj$ace))
crop.obj$ace[node,]
# 53% probability of free-roosting ancestor for the combined clade of Rhilophidae and Hipposideridae


node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Hipposideridae' | families=='Rhinolophidae' | families=='Pteropodidae'))),rownames(crop.obj$ace))
crop.obj$ace[node,]
# 54% probability of free-roosting ancestor for the leap between Horshoe and Flying fox

roster <-names(which(families=='Molossidae' | families=='Vespertilionidae' ))
roster<-roster[which(roster %in% pruned.tree$tip==T)]

node <- match(getMRCA(crop.mcmc100[[1]], roster ),rownames(crop.obj$ace))
crop.obj$ace[node,]
# 79% probability of specialised ancestor in leap between Vespertilionidae and Molossidae

# We might ask how these probabilities change on the condition that the root is non specialized




pca<-prcomp(x)
i <- max(which( ((pca$sdev^2)/sum(pca$sdev^2)) >0.05))
data <- pca$x[,1:i]
# The first four axes of sternum shape variety are selected as a truncation. 
# Brownian Motion, Early Burst, and stochastic-map contingent multiple-rates Brownian Motion models will be fit.
# AICc values will be calculated for these models. 
# This will reveal whether estimates of the 

favour <- matrix(NA,100,3)
root <- matrix(NA,100,1)
# Matrices to store model outcomes, and the inferred root state of roosting ecology 

library(mvMORPH) # v 1.2.1
# Package for evolutionary model fitting 

j<-1
	# Brownian motion
	BM <-  tryCatch({mvMORPH::mvBM(tree=crop.mcmc100[[j]], data=data, model="BM1", optimization='Nelder-Mead')},
	error=function(e){})
	# Other optimization procedures reach solutions with the same AIC, but do not necessarily support local concavity
	
 	# Early Burst
	EB <- tryCatch({mvMORPH::mvEB(tree=crop.mcmc100[[j]], param=list(up=0 ), data=data, optimization='Nelder-Mead')},
	error=function(e){}) # -1014.653 
	# Other optimization procedures either fail to converge or converge on beta=0 (Brownian walk)
	
	favour[,1]<-BM$AICc ; favour[,3]<-EB$AICc


# Compute multiple rate Brownian motion fits for all 100 stochastic character estimates
# Depending on your hardware this may take a while 
for(j in 1:100){
	BMM <- tryCatch({mvMORPH::mvBM(tree=crop.mcmc100[[j]], data=data, model="BMM", optimization= "Nelder-Mead")},
	error=function(e){}) 
	if(is.null(BMM)==F){
		favour[j,2] <- BMM$AICc
		rm(BMM)
	}
	print(j)
	
}

for(j in 1:100){
	root[j] <- getStates(crop.mcmc100[[j]])[1]
}
# Record the roosting ecology that is estimated at the root of the family tree 


obj.2 <- summary(mcmc100[which(root=='free')])
crop.mcmc.free <- list()
for(i in 1:length(mcmc100[which(root=='free')])){
	crop.mcmc.free[[i]] <- keep.tip.simmap(phy=mcmc100[which(root=='free')][[i]],rownames(x)) 
}
class(crop.mcmc.free)<-class(mcmc100[which(root=='free')])
crop.obj.free <- summary(crop.mcmc.free)
# Cropping to taxa for which sternum data are available 

crop.obj.free$ace[1,]

node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Pteropodidae'))),rownames(crop.obj$ace))
crop.obj.free$ace[node,]
# a non specialized root all but guarantees that the ancestors of pteropodidae lacked specializations 

node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Hipposideridae' | families=='Rhinolophidae'))),rownames(crop.obj$ace))
crop.obj.free$ace[node,]
# a non specialized root makes it highly likely the ancestors of Rhinolophoidea lacked specializations 


node <- match(getMRCA(crop.mcmc100[[1]], c('Pteronotus_davyi',names(which(families=='Phyllostomidae')) )),rownames(crop.obj$ace))
crop.obj.free$ace[node,]
# a non specialized root makes it highly likely the ancestors of Rhinolophoidea lacked specializations 


node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Hipposideridae' | families=='Rhinolophidae' | families=='Pteropodidae'))),rownames(crop.obj$ace))
crop.obj.free$ace[node,]
# Division between Pteropodidae and Rhinolophoidea definitely without specializations


roster <-names(which(families=='Molossidae' | families=='Vespertilionidae' ))
roster<-roster[which(roster %in% pruned.tree$tip==T)]
node <- match(getMRCA(crop.mcmc100[[1]], roster ),rownames(crop.obj$ace))
crop.obj.free$ace[node,]
# Equivocal evidence of roosting specializations at the root of the division between Vespertilionidae and Molossidae 


roster <-c('Pteronotus_davyi',names(which(families=='Molossidae' )))
roster<-roster[which(roster %in% pruned.tree$tip==T)]
node <- match(getMRCA(crop.mcmc100[[1]], roster ),rownames(crop.obj$ace))
crop.obj.free$ace[node,]
# a non specialized root makes it highly likely the ancestors of Rhinolophoidea lacked specializations 



# A chi-squared test on a contingency table of the co-occurrences of model support and
# root state will reveal whether the assumed root state influences model preference

support <- (favour[,2]-favour[,3]) < -2
# A delta AICc of 2 indicates a substantial difference in model preference

outcomes <- table(paste(root,support))
Xi <- as.matrix( rbind(c(outcomes['compressed FALSE'], outcomes['compressed TRUE']  ), c(outcomes['free FALSE'], outcomes['free TRUE']  ) ))
Xi[is.na(Xi)] <-0

colnames(Xi) <- c('No','BMM')
rownames(Xi) <- c('compressed','free')
chisq.test(Xi,simulate.p.value = TRUE )
# When the root value is 'free/exposed' there is a stronger support for a faster rate of evolution in
# the free-roosting group. This is probably because there are strong inter-familial differences in 
# shape that are established early in the radiation of bats. 
# Decelleration of phenotypic evolution might actually causally result from the acquisiton of constraints 
# associated with roosting in tight recessess. 



# The results will now be presented visually; the following will produce Figure 2 in the main manuscript. 
library( ggplot2 ) # v 3.4.1
library(cowplot) # v 1.1.3
# Plotting function package
library( ggdendro ) # v 0.1.23
library( dendextend ) # v 1.17.1
library(zoo) # v 1.8-12
library(dplyr) # v 1.1.1
library(viridis) # v 0.6.5 
# Packages for plotting and data curation to prepare it for plotting 

	# Prepare phylogenetic tree for plotting 
	dendr <- ggdendro::dendro_data(stats::as.dendrogram(crop.obj$tree[[1]]))
	lab.dat <- dendr$labels

	fit <- as.matrix(crop.obj$ace[,1]/rowSums(crop.obj$ace))
	dendr.mod<-dendr$segments/2 # This is the tree structure 
	dendr.mod$z <- rep(NA,dim(dendr.mod)[1])
	tips <- which(dendr$segments$yend==0)
	nodes <- which(dendr$segments$y==dendr$segments$yend)
	internal.branches <- which(dendr$segments$x==dendr$segments$xend & dendr$segments$yend!=0)
	dendr.mod$z[tips] <- fit[match(lab.dat$label,rownames(fit))][match(dendr.mod$x[tips]*2,lab.dat$x)]
	dendr.mod[nodes[seq(1,length(nodes),2)] ,]$z <- fit[1:((length(fit)-1)/2)]
	dendr.mod$z[seq(3,dim(dendr.mod)[1],by=4)]<-dendr.mod$z[seq(1,dim(dendr.mod)[1],by=4)]
	dendr.mod$z[which(is.na(dendr.mod$z)==T)]<-dendr.mod$z[which(is.na(dendr.mod$z)==T)-1]

	# Prepare a dataset of evolutionary path lengths for species across the family tree
	bars <- data.frame( cbind(lab.dat[,1:2]/2, (magxc/max(magxc))[lab.dat$label] ))
	colnames(bars)[3]<-'bar'


 Subplot.title.font.size <- 10 
 Subplot.axis.font.size <- 9 
 Legend.axis.font.size <- 7.5 
 Legend.title.font.size <- 9 
 Legend.tick.size <- 1/4
 Subplot.module.linewidth <- 1/4
 Subplot.tick.size <-1/4
 Subplot.border.linewidth <- 1/4
 legend.key.width <- 1/2
 legend.key.height <- 0.5
# Plot parameters




breaky <- quantile(bars$bar[which(compression[rownames(bars)]==1)]^2,probs=(0.95))
extreme <- rownames(bars)[which((bars$bar^2)>breaky)]
# Which bats are in the top 5% of most unusual sternum shape given their position in the bat family tree?

breaky <- c(mean(bars$bar)^2,(mean(bars$bar)+2*sd(bars$bar))^2)
# What are the sternum divergence magnitudes which are average, and what is the level which is 2 standard deviations beyond this? 

library(extrafont) # 0.19 
# Used to plot with Times New Roman font as per ProcB guidelines 


# We're going to produce a plot of a family tree with our summary stochastic character map of roosting ecology super-imposed on it.
# The tips of the family tree with mesh with a bar plot which documents how divergent each species' sternum shape is, compared to its position in the family tree. 
# We need to extract sub-components of the legend from this plot, and a direct way to achieve this is to simply re-compile the plot 3 times and harvest
# the components independently:

bargram <- 
	ggplot()+
	geom_segment(data=dendr.mod, aes(x = -y/2, y = x, xend = -yend/2, yend =xend), lwd=4*Subplot.module.linewidth,lineend='round',col='black')+ # Subplot.module.linewidth*4
	geom_segment(data=dendr.mod, aes(x = -y/2, y = x, xend = -yend/2, yend =xend, col=z ), lwd=2*Subplot.module.linewidth,lineend='round')+# Subplot.module.linewidth*2
	
	scale_colour_viridis( labels=c('Exposed','Enclosed'), breaks=c(0,1),direction=-1,option='D')+
	labs(color = "Roost:")+
	 ggnewscale::new_scale_color() +
	geom_segment(data=bars, aes(x = .2, y = x, xend = .2+20*bar^2, yend =x , col=bar^2), lwd=4*Subplot.module.linewidth, show.legend = F)+ # Subplot.module.linewidth*4
 	scale_colour_gradient(high = "blue", low = "cyan",breaks=breaky, label=c(expression(bar(sqrt(D))),expression(paste("2", sigma))))+
	labs(color = "D:")+
	geom_text(data=bars[extreme,],aes(x= 1.2+20*bar^2, y=x+0.5), label= 12-rank(bars[extreme,]$bar) )+
	ylab('Mya')+
	scale_y_continuous(breaks=c(-60,-50,-40,-30,-20,-10,0))+
    	theme(legend.title=element_text(size=Legend.title.font.size,family = "Times New Roman"),
	legend.text=element_text(size=Legend.title.font.size,family = "Times New Roman") )

# Extract legend distinguishing the colour palette of exposed/enclosed roosting 
Roost.legend <- cowplot::get_legend(bargram)

bargram <- 
	ggplot()+
	geom_segment(data=dendr.mod, aes(x = -y/2, y = x, xend = -yend/2, yend =xend), lwd=4*Subplot.module.linewidth,lineend='round',col='black')+ # Subplot.module.linewidth*4
	#geom_segment(data=dendr.mod, aes(x = -y/2, y = x, xend = -yend/2, yend =xend, col=z), lwd=2*Subplot.module.linewidth,lineend='round', show.legend = F)+# Subplot.module.linewidth*2
	scale_colour_viridis( labels=c('Exposed','Enclosed'), breaks=c(0,1),direction=-1,option='D')+
	labs(color = "Roost:")+
	 ggnewscale::new_scale_color() +
	geom_segment(data=bars, aes(x = .2, y = x, xend = .2+20*bar^2, yend =x , col=bar^2), lwd=4*Subplot.module.linewidth, show.legend = T)+ # Subplot.module.linewidth*4
 	scale_colour_gradient(high = "blue", low = "cyan",breaks=breaky, label=c(expression(bar(sqrt(D))),expression(paste("2", sigma))))+
	labs(color = "D:")+
	geom_text(data=bars[extreme,],aes(x= 1.2+20*bar^2, y=x+0.5), label= 12-rank(bars[extreme,]$bar) )+
	ylab('Mya')+
	scale_y_continuous(breaks=c(-60,-50,-40,-30,-20,-10,0))+
    	theme(legend.title=element_text(size=Legend.title.font.size,family = "Times New Roman"),
	legend.text=element_text(size=Legend.title.font.size,family = "Times New Roman") )

# Extract the legend of evolutionary path lengths
D.legend <- cowplot::get_legend(bargram)


bargram <- 
	ggplot()+
	geom_segment(data=dendr.mod, aes(x = -y/2, y = x, xend = -yend/2, yend =xend ), lwd=4*Subplot.module.linewidth,lineend='round',col='black')+
	geom_segment(data=dendr.mod, aes(x = -y/2, y = x, xend = -yend/2, yend =xend, col=z ), lwd=2*Subplot.module.linewidth,lineend='round', show.legend=F)+
	scale_colour_viridis( labels=c('Exposed','Enclosed'), breaks=c(0,1),direction=-1,option='D')+
	labs(color = "Roost:")+
	ggnewscale::new_scale_color() +
	geom_segment(data=bars, aes(x = .2, y = x, xend = .2+20*bar^2, yend =x , col=bar^2), lwd=4*Subplot.module.linewidth, show.legend = F)+
 	scale_colour_gradient(high = "blue", low = "cyan",breaks=breaky, label=c(expression(bar(sqrt(D))),expression(paste("2", sigma))))+
	labs(color = "D:")+
	geom_text(data=bars[extreme,],aes(x= 2+20*bar^2, y=x+0.5), label= 12-rank(bars[extreme,]$bar),family = "Times New Roman" )+
	ylab('Mya')+
	scale_y_continuous(breaks=c(-60,-50,-40,-30,-20,-10,0))+
   	theme(legend.position='right',
	legend.title=element_text(size=Legend.title.font.size),
	legend.text=element_text(size=Legend.title.font.size),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='white', color='white'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.text=element_blank(),
	axis.ticks=element_blank(),
	axis.title=element_blank(),
	legend.spacing.x = unit(1.5,'cm'),
	legend.spacing.y = unit(1.5,'cm'),
	legend.key.width=unit(.25,'cm'),
	legend.key.height=unit(.25,'cm'), 
	legend.key.spacing.x=unit(-.0,'cm'),
	legend.key.spacing.y=unit(-.0,'cm'),
	legend.background = element_rect(fill='transparent'),
	text = element_text(family = "Times New Roman"))
# A version of the plot of evolutionary path lengths which omits the legends. We will arrange all components together in an aesthetic way.


# We're now going to visualise a possible evolutionary history of sternal morphology trait values, reconstructed assuming a Brownian-process variance-covariance structure,
# of sternum shape over the leading two principal components. 

# Wrapper for ancestral state estimation at nodes: 
anc.recon <-function(tree,data,mod,axis){
	dendr <- dendro_data(as.dendrogram(tree))
	lab.dat <- dendr$labels
	fit<-mvMORPH::estim(tree, data= data, object=mod, asr=TRUE)
	dendr.mod<-dendr$segments/2 
	# This is the tree structure 
	tips <- which(dendr$segments$yend==0)
	nodes <- which(dendr$segments$y==dendr$segments$yend)
	internal.branches <- which(dendr$segments$x==dendr$segments$xend & dendr$segments$yend!=0)
	dendr.mod$x[tips] <- fit$estim[ match( tree$edge[,1][match(dendr$segments$xend[tips],tree$edge[,2])],rownames(fit$estim) ),axis ] 
	dendr.mod$xend[tips] <- data[dendr$segments$xend[tips],axis]
	# This has aligned the tips so that their heads match observation and their tails match ancestral estimates
	dendr.mod[nodes[seq(1,length(nodes),2)] ,]$x <- fit$estim[,axis]
	dendr.mod[nodes[seq(1,length(nodes),2)]+1 ,]$x <- fit$estim[,axis]
	dendr.mod$xend[internal.branches] <- dendr.mod$x[internal.branches+3]
	# This gets all of the tips, their nodes and half of the internal branches in place
	dendr.mod[nodes[seq(1,length(nodes),2)]+3 ,]$x <- fit$estim[,axis]
	dendr.mod$x[nodes] <- dendr.mod$xend[nodes] <-NA
	for(i in internal.branches){
		# Conditions
		x.match <- which(dendr$segments$x==dendr$segments$xend[i] )
		# which tails of other lines have the same horizontal position as the head of the branch?
		y.match <- which(dendr$segments$y==dendr$segments$yend[i] )
		# which tails of other lines have the same vertical position as the head of the branch?
		head.nodes <- intersect(x.match,y.match)
		# These are the nodes that lead to the next vertical branches
		head.nodes.x.match <- match(dendr$segments$xend[head.nodes], dendr$segments$x)
		head.nodes.x.match <- head.nodes.x.match[which(head.nodes.x.match%in%y.match==T)]
		# which tails of vertical branches meet the head nodes?
		dendr.mod$xend[i] <- dendr.mod$x[head.nodes.x.match][1]	
	}
	return(dendr.mod)
}

# Initialise a Brownian model process estimate: 
j <-1 
BM <- tryCatch({mvMORPH::mvBM(tree=crop.mcmc100[[j]], data=data, model="BM1", optimization="Nelder-Mead")},
error=function(e){}) 

i<-4
PC1 <- anc.recon(tree=crop.mcmc100[[j]],data=pca$x[,1:i],mod=BM,axis=1)
PC2 <- anc.recon(tree=crop.mcmc100[[j]],data=pca$x[,1:i],mod=BM,axis=2)
# Estimate the projection of our node-wise ancestral state imputations upon PC1 and PC2 of sternum shape variety. 

# Wrapper function to colour the family tree according to our summary of roosting ecology evolution. 
paint.colors <-function(tree,fit){
	dendr <- ggdendro::dendro_data(as.dendrogram(tree))
	lab.dat <- dendr$labels
	dendr.mod<-dendr$segments/2 # This is the tree structure 
	dendr.mod$z <- rep(NA,dim(dendr.mod)[1])
	tips <- which(dendr$segments$yend==0)
	nodes <- which(dendr$segments$y==dendr$segments$yend)
	internal.branches <- which(dendr$segments$x==dendr$segments$xend & dendr$segments$yend!=0)
	dendr.mod$z[tips] <- fit[match(lab.dat$label,rownames(fit))][match(dendr.mod$x[tips]*2,lab.dat$x)]
	dendr.mod[nodes[seq(1,length(nodes),2)] ,]$z <- fit[1:((length(fit)-1)/2)] # This line good 
	dendr.mod$z[seq(3,dim(dendr.mod)[1],by=4)]<-dendr.mod$z[seq(1,dim(dendr.mod)[1],by=4)]
	dendr.mod$z[which(is.na(dendr.mod$z)==T)]<-dendr.mod$z[which(is.na(dendr.mod$z)==T)-1]
	branch.color <- dendr.mod$z
	return(branch.color)
}

# Colour the tree
branch.color <- paint.colors(tree=force.ultrametric(crop.mcmc100[[j]]), fit=as.matrix(crop.obj$ace[,1]/rowSums(crop.obj$ace)) )

#Euclidean distance function
euclid<-function(camera,stack){
	ones <- matrix(1,dim(stack)[1],1)
	camera <- ones%*%camera
 	sqrt(rowSums((stack[,1:3]-camera)^2)) 
}

# You can ignore the mechanics of this function; we're going to view the ancestral state imputation from a cool angle. 
make.3D <- function(x,y, scale, alpha, gamma, camera){
	starts <- cbind(x[complete.cases(x),1],y[complete.cases(y),1:2])
	stops <- cbind(x[complete.cases(x),3],y[complete.cases(y),3:4])

	Zc <-( diff(range(c(starts[,3],stops[,3])))/diff(range(c(starts[,1],stops[,1]))) )
	Yc <-( diff(range(c(starts[,2],stops[,2])))/diff(range(c(starts[,1],stops[,1]))) )
	Zc<-Zc*scale

	starts[,2]<-starts[,2]/Yc
	stops[,2]<-stops[,2]/Yc
	starts[,3]<-starts[,3]/Zc 
	stops[,3]<-stops[,3]/Zc 

	starts.shadow<-starts
	stops.shadow<-stops

	x1 <- min(c(starts[,1],stops[,1]))
	x2 <- max(c(starts[,1],stops[,1]))
	y1 <- min(c(starts[,2],stops[,2]))
	y2 <- max(c(starts[,2],stops[,2]))
	z1 <- min(c(starts[,3],stops[,3]))
	z2 <- max(c(starts[,3],stops[,3]))

	axes.starts <- cbind(c(x1,x2,x2,x1, x2, x1,x2,x2,x1),c(y1,y1,y2,y2, y2, y1,y1,y2,y2),c(z2,z2,z2,z2, z1,  z1,z1,z1,z1))
	axes.stops <- cbind(c(x2,x2,x1,x1, x2, x2,x2,x1,x1),c(y1,y2,y2,y1, y2, y1,y2,y2,y1),c(z2,z2,z2,z2, z2, z1,z1,z1,z1))

	alpha <- (alpha*pi)/180
	gamma <- (gamma*pi)/180

	rot.x <- matrix(NA,3,3)
	rot.x[1:9] <- c(1,0,0,0,cos(gamma),-sin(gamma),0,sin(gamma),cos(gamma))
	rot.z<- matrix(NA,3,3)
	rot.z[1:9] <- c(cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1)

	starts <- as.matrix(starts)[,1:3]%*%rot.z%*%rot.x
	stops <- as.matrix(stops)[,1:3]%*%rot.z%*%rot.x

	rotated <- data.frame(cbind(starts[,c(1,2,3)],stops[,c(1,2,3)]))
	colnames(rotated) <- c('x','y','z','xend','yend','zend')

	axes.starts <- as.matrix(axes.starts)[,1:3]%*%rot.z%*%rot.x
	axes.stops <- as.matrix(axes.stops)[,1:3]%*%rot.z%*%rot.x

	axes.rotated <- data.frame(cbind(axes.starts[,c(1,2,3)],axes.stops[,c(1,2,3)]))
	colnames(axes.rotated) <- c('x','y','z','xend','yend','zend')

	stoppies <- c(0.1734556, 0.3469112, 0.5203668, 0.6938224, 0.8672780)*1

	time.starts <- cbind(rep(x1,7),rep(y1,7),c(z1,z1,z1+(stoppies[1]*z2),z1+(stoppies[2]*z2),z1+(stoppies[3]*z2),z1+(stoppies[4]*z2),z1+(stoppies[5]*z2) ))
	time.stops <- cbind(c(x1,rep(x1+x1/20,6)),c(y1,rep(y1+y1/20,6)),c(z2,z1,z1+(stoppies[1]*z2),z1+(stoppies[2]*z2),z1+(stoppies[3]*z2),z1+(stoppies[4]*z2),z1+(stoppies[5]*z2) ))
	time.starts <- as.matrix(time.starts)[,1:3]%*%rot.z%*%rot.x
	time.stops <- as.matrix(time.stops)[,1:3]%*%rot.z%*%rot.x

	time.rotated <- data.frame(cbind(time.starts[,c(1,2,3)],time.stops[,c(1,2,3)]))
	colnames(time.rotated) <- c('x','y','z','xend','yend','zend')

	arrow.starts <- cbind(c(x2*.95,x2),c(y1,y1+(y2*.05)),c( z1-(z2*0.1),z1-(z2*0.1) ) )
	arrow.stops <- cbind(c(x1+(x2*0.1),x2),c(y1,(y2*.90)),c( z1-(z2*0.1),z1-(z2*0.1) ) )
	arrow.starts <- as.matrix(arrow.starts)[,1:3]%*%rot.z%*%rot.x
	arrow.stops <- as.matrix(arrow.stops)[,1:3]%*%rot.z%*%rot.x

	arrow.rotated <- data.frame(cbind(arrow.starts[,c(1,2,3)],arrow.stops[,c(1,2,3)]))
	colnames(arrow.rotated) <- c('x','y','z','xend','yend','zend')

	camera <- t(matrix(cbind(x1,y2,0)))

	camera.rotated <- data.frame(camera%*%rot.z%*%rot.x)
	rotated.euclid <- euclid(as.matrix(camera.rotated),rotated[,4:6]) 

	starts.shadow[,3] <- max(c(starts.shadow[,3],stops.shadow[,3]))
	stops.shadow[,3] <- max(c(starts.shadow[,3],stops.shadow[,3]))

	starts.shadow <- as.matrix(starts.shadow)[,1:3]%*%rot.z%*%rot.x
	stops.shadow <- as.matrix(stops.shadow)[,1:3]%*%rot.z%*%rot.x
	rotated.shadow <- data.frame(cbind(starts.shadow[,c(1,2,3)],stops.shadow[,c(1,2,3)]))
	colnames(rotated.shadow) <- c('x','y','z','xend','yend','zend')

	base <- t(matrix(cbind(0,0,z2)))
	base.rotated <- data.frame(base%*%rot.z%*%rot.x)
	base.euclid <- euclid(as.matrix(base.rotated),rotated.shadow[,1:3])


	return <- list()
	return[[1]] <- rotated.shadow
	return[[2]] <- rotated.euclid
	return[[3]] <- axes.rotated
	return[[4]] <- rotated
	return[[5]] <- camera.rotated
	return[[6]] <- base.rotated
	return[[7]] <- time.rotated
	return[[8]] <- arrow.rotated
	names(return) <- c('rotated.shadow','rotated.euclid','axes.rotated','rotated',
	'camera.rotated','base.rotated','time.rotated','arrow.rotated')
	return(return)
}


# Angles to rotate our data by 
alpha<-150; gamma<-65
alpha <- (alpha*pi)/180
gamma <- (gamma*pi)/180

# Various plotting angle statements
rot.x <- matrix(NA,3,3)
rot.x[1:9] <- c(1,0,0,0,cos(gamma),-sin(gamma),0,sin(gamma),cos(gamma))
rot.z<- matrix(NA,3,3)
rot.z[1:9] <- c(cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1)
rot.points <- cbind(data[,1:2],1)%*%rot.z%*%rot.x

# Perform data rotation
plotty <- make.3D(scale=5,alpha=150,gamma=65,x=PC1,y=PC2 ) 

# Various aesthetic plotting housekeeping jobs
shade <- (plotty$rotated.euclid^.2)
shade <- shade/max(shade)
base.shadow <- rgb(shade,shade,shade)
super.shadow <- rgb(shade*.9,shade*.9,shade*.9)
col.v<-branch.color[seq(2,dim(PC1)[1],by=2)][rev(order(plotty$rotated.euclid))] 
plot.ord <- plotty$rotated[rev(order(plotty$rotated.euclid)),]
col.ord<-rep(1,length(col.v)*2)
col.ord[seq(2,length(col.ord),by=2)]<- branch.color[seq(2,dim(PC1)[1],by=2)][rev(order(plotty$rotated.euclid))] 
line.ord <- rep(2*plotty$rotated.euclid[(order(plotty$rotated.euclid))]^.6,each=2)
line.ord[seq(1,length(col.ord),by=2)] <- line.ord[seq(1,length(col.ord),by=2)]+.5
testy<-make.3D(scale=5,alpha=150,gamma=65,x=PC1[match(data[,1],PC1$xend),],y=PC2[match(data[,1],PC1$xend),] ) 

# Make the plot of the reconstructed evolutionary history
# We expect to see that divergent taxa tend to roost in exposed environments, and that as enclosed roosting ecologies accumulate through time
# the manifest variance of sternal evolution tends to decline and taxa distributions will begin to look bunchy. 
phenogram<-
ggplot()+
	geom_segment(data=plotty$rotated.shadow[rev(order(plotty$rotated.euclid)),],
  	aes(x = x, y = y, xend = xend, yend =yend), linejoin = "round", lineend = "round",col=base.shadow[rev(order(plotty$rotated.euclid))],
	lwd= .5*4*plotty$rotated.euclid[(order(plotty$rotated.euclid))]^.25 )+ # 2
	geom_segment(data=plotty$rotated.shadow[rev(order(plotty$rotated.euclid)),],
  	aes(x = x, y = y, xend = xend, yend =yend), linejoin = "round", lineend = "round",
	col=super.shadow[rev(order(plotty$rotated.euclid))],lwd= .5*3*plotty$rotated.euclid[(order(plotty$rotated.euclid))]^.5 )+ # 1*
 	geom_segment(data=plotty$axes.rotated,
	aes(x = x, y = y, xend = xend, yend =yend), col=c('grey','grey',NA,NA, NA,  NA,NA,NA,NA),lwd=Subplot.module.linewidth,lineend='round')+ #1/4
 	geom_segment(data=plotty$time.rotated,
	aes(x = x, y = y, xend = xend, yend =yend), col='black',lwd=Subplot.module.linewidth,lineend='round')+ # 1/4
 	geom_text(data=plotty$time.rotated,
	aes(x = xend+0.01, y = yend),label=c('','0','10','20','30','40','50'), size=Subplot.axis.font.size/.pt,family = "Times New Roman")+
	geom_text(aes(x = 0.2308668, y = 0.022),label=c('Mya'), size=Subplot.axis.font.size/.pt,family = "Times New Roman")+
 	geom_segment(data=plotty$arrow.rotated,
	aes(x = x, y = y, xend = xend, yend =yend), col='black',lwd=Subplot.module.linewidth,lineend='round',arrow = arrow(length = unit(0.3, "cm")))+ #1/4
 	geom_text(data=plotty$arrow.rotated,
	aes(x = -0.005+(x+xend)/2, y = 0.008+(y+yend)/2), label=c('PC1','PC2'), angle=c(-25,50), size=Subplot.axis.font.size/.pt ,family = "Times New Roman")+
	geom_segment(
  	aes(x = rep(plot.ord$x,each=2), y = rep(plot.ord$y,each=2), xend = rep(plot.ord$xend,each=2), yend = rep(plot.ord$yend,each=2), col=col.ord )  ,
	 lwd= line.ord*1 ,  linejoin = "round", lineend = "round",show.legend=F)+ # line.ord/2
	scale_colour_viridis( labels=c('free','tight'), breaks=c(0,1),direction=-1,option='D')+
	geom_segment(data=plotty$axes.rotated,
  	aes(x = x, y = y, xend = xend, yend =yend), col=c(NA,NA,'black'  ,'black'  , NA  ,  'black','black','black','black'),lwd=Subplot.module.linewidth,lineend='round')+ #1/4
	geom_point(data=testy$rotated[match(extreme,rownames(data)),], aes(x=xend,y=yend), col='black',size=6, shape=21, stroke=1/2, fill='white')+
	# We need to put text in the bubbles:
	geom_text(data=testy$rotated[match(extreme,rownames(data)),], aes(x=xend,y=yend), label=12-rank(bars[extreme,]$bar),size=Subplot.axis.font.size/.pt,fontface='bold',family = "Times New Roman")+
	theme( panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color='white'),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
	   axis.ticks=element_blank(),
		axis.text=element_blank(),
		axis.title=element_blank(),
	text = element_text(family = "Times New Roman")
	)


# Various sundry plotting housekeeping follows

# Make a blank background:
blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_rect(fill='white',color='white'),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())




# Arrange all plot sub-components
  ggdraw() +
  draw_plot(blank) +
  draw_plot(phenogram, x = 0.0, y = 0.0, width = .6, height = 1)+
  draw_plot(bargram, x = 0.66, y = 0.0, width = .35, height = 1)+
  draw_grob(D.legend, x=0.6, y=0.75, width=.05, height=.05)+
  draw_grob(Roost.legend, x=0.58, y=0.2, width=.05, height=.05)+
  draw_text(text=c('a','b'), x=c(0.05,0.65), y=c(0.95,0.95), size=12,family = "Times New Roman")+ # fontface='italic'
  theme(plot.background = element_rect(fill='white', color = NA))


# Save the output; adjust the filepath as necessary 
ggsave(filename='Whatever_name_you_want.png',width=18,height=11,unit='cm',dpi=600,device=png, type='cairo')



