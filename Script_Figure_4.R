# The purpose of this script is to generate Figure 4 from the manuscript Orkney et al., 2026
# 'Declining rates of evolution and limited convergence in bat
# sternum shape can be explained by a ratchet of specialisation'
# The analyses will determine whether roosting in enclosed spaces represents a constraint on bat thoracic skeleton
# evolution and produce Figure 3 in the main manuscript of Orkney et al., 2026. 
# The available data are sternal shape landmark constellations across a diversity of bat species.
# An estimate of the shared ancestry between species.
# A species-wise binary bag of words describing roosting ecology, which can be collapsed into a dichotomous
# 'exposed/enclosed' categorization.  (We use the language 'compressed/free' interchangeably in comments here.) 
# An estimate of the shared ancestry between species is available (Shi & Rabosky 2015 https://doi.org/10.1111/evo.12681) 
# A species-wise binary bag of words is available from Orkney et al., 2025 (https://doi.org/10.1038/s41559-024-02572-9)
# This script was written by Dr. A Orkney and the final version was compiled on July 19th 2026. 
#
# We will, explicitly, fit a multiple Ornstein-Uhlenbeck model to determine whether gravitating selection 
# towards distinct combinations of sternum shape traits and flight-style ecological properties 
# is a likely explanation of bat sternum shape evolution through time. 
#
# We hypothesize that the gradual settlement of bats within these adaptive zones, consistent with a Simpsonian
# Quantum evolutionary dynamic of adaptive radiation, is influences by second-order adaptive demands imposed
# by the relationship between sternum shape variety and exposed/enclosed roosting ecology. 
#
# Specifically, we expect that enclosed roosting ecologies decrease the capacity of sternum shape evolutionary change
# and that this changes the gradient descent of bat lineages towards adaptive optima- allowing them to be more easily
# captured by the field of stabilizing selection surrounding different flight-style adaptive zones. 
# 
# We will further an estimate of ancestral sternum shape states as a series of discrete quantum leaps between
# different states attested in living bat species- with the hypothesis that the a reconstruction of this series of 
# evolutionary events is likelier if we assume stationary variances of the evolutinary process associated
# with exposed roosting ecologies. 



setwd()
# Set the directory to the location of the landmark data
# This will change if you download the data to a personal computer

load('sternum_array_sep_27_2024.RData')
# load the data

taxa <- dimnames(sternum.array)[[3]]
# This is a vector of available bat species. 

# Set the work directory to the location with the bat metadata. 
metadata <- read.csv('Bat_CT_process_list_Andrew_only.csv')
families <- metadata$Family[ match(taxa,metadata$Shi_match) ]
names(families) <- taxa
original.names <- paste(metadata$Genus[match(taxa,metadata$Shi_match)], metadata$Species[match(taxa,metadata$Shi_match)],sep='_')
names(original.names)<-taxa
# Substitute names in our collected taxa with congeners on the phylogeny of Shi & Rabosky 2015, if 
# a direct match is not available. 

library(ape) # v 5.8-1
# Package for managing family trees

library(phytools) #  v 2.4-4
# Package for managing family trees

bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# This phylogeny far exceeds the number of taxa for which we have landmark constellations
# we must therefore prune the phylogeny 

pruned.tree <- keep.tip(bat.tree,dimnames(sternum.array)[[3]])
# prune the bat tree to the taxa of interest

pruned.tree <- drop.tip(pruned.tree, c('Micropteropus_pusillus','Molossops_temminckii'))
# I don't have ecological metadata for Micropteropus pussilus
# It is a nectar and fruit loving bat, that seems to prefer tropical rainforest
# margins in subsaharan Africa, but I did not find a clear description of
# its prefered roosts. 
pruned.tree <- phytools::force.ultrametric(pruned.tree)

# 70 bat species. Thank you Elizabeth Augustin, thank you Beyonca Akers. 

# We need to align the landmark constellations into a common reference frame:


library(geomorph) # v 4.0.10 
# Package for shape data management 

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
# These are the aligned coordinates

x<- geomorph::two.d.array(coords)
# It will be convenient to treat the coordinates as a 2-D array

pca<-prcomp(x)
i <- max(which( ((pca$sdev^2)/sum(pca$sdev^2)) >0.05))
data <- pca$x[,1:i]
# A singular value decomposition is performed and truncated 

# We need to load roosting style categorizations 
metadata <- read.csv('Bat_eco_metadata.csv')

eco <- metadata
# The ecological data has been subset to the species of interest

binary.roost.styles <- eco[,c(6:16)]
# We actually have 112 taxa we can use for imputation of free-roosting and tight-roosting evolutionary history
# even though we oly have 71 valid sterna. 

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
# and crevices etc is considered 'tight' 
# I do not consider caves or artificial roosts to represent compressed environments unless they are
# accompanied by the word 'crevice'. 
# Most bats roost in compressed environments.

names(compression)<- eco$Shi
# Name the vector


roost<-setNames(compression,names(compression))
roost[which(roost==1)]<-'compressed'
roost[which(roost==0)]<-'free'
roost<-factor(roost)
# Transform roosting into a factor variable

library(phytools)
eco.tree <- phytools::force.ultrametric(ape::keep.tip(bat.tree,names(roost)))
# It is necessary that the tree is ultrametric for the next steps of analysis 
roost <- roost[eco.tree$tip]


library(mvMORPH) # v 1.2.1
# Package for evolutionary model fitting 


# Assign different bat groups to different hypothesized adaptive zones sympathetic with flight-style and sternum shape properties
# across the group

painted.tree <- pruned.tree
node <- findMRCA(pruned.tree, tips=names(families[which(families=='Rhinolophidae' | families=='Hipposideridae' | families=='Megadermatidae' | families=='Craseonycteridae')]) )
painted.tree <- paintSubTree(painted.tree,node,state='Horseshoe',stem=T)

node <- findMRCA(pruned.tree, tips=names(families[which(families=='Emballonuridae' | families=='Noctilionidae'| families=='Mormoopoidae' | families=='Thyropteridae' | families=='Phyllostomidae')]) )
painted.tree <- paintSubTree(painted.tree,node,state='Grade',stem=T)

node <- findMRCA(pruned.tree, tips=names(families[which(families=='Phyllostomidae')]) )
painted.tree <- paintSubTree(painted.tree,node,state='Phyllostomidae',stem=T)

node <- findMRCA(pruned.tree, tips=names(families[which(families=='Molossidae')])[1:3] )
painted.tree <- paintSubTree(painted.tree,node,state='Molossidae',stem=T)

node <- findMRCA(pruned.tree, tips=names(families[which(families=='Vespertilionidae')]) )
painted.tree <- paintSubTree(painted.tree,node,state='Vespertilionidae',stem=T)

plot(painted.tree)
# >> Approach
# We will estimate the Optima for the different families as their ancestral state

unique.states <- unique(getStates(painted.tree, type='tips')) 
anc <- matrix(NA, length(unique.states),dim(data)[2]) 
rownames(anc)<-unique.states
for(i in 1:length(unique.states)){
	species <- names(which(getStates(painted.tree, type='tips')==unique.states[i]))
	ones<- matrix(1,length(species),1)
	Csub <- ape::vcv.phylo(keep.tip(pruned.tree,species))
	anc[i,] <- (t(ones)%*%solve(Csub)%*%data[species,])/sum(solve(Csub))
}
# We shall allow ourselves to assume that the ancestral states of each adaptive zone are good approximations of the 
# trait optima associated with them 

# Most recent common ancestor of each adaptive zone we hypothesize.
MRCAs <- matrix(NA,length(unique.states),1)
for(i in 1:length(unique.states)){
	MRCAs[i,] <- getMRCA(pruned.tree, names(which(getStates(painted.tree, type='tips')==unique.states[i])) )
}

excluded.nodes <- 
 apply(MRCAs, 1, getDescendants, tree=pruned.tree)
excluded.nodes <-  c(MRCAs, unique(unlist(excluded.nodes)))
all.nodes <- unique(as.vector(pruned.tree$edge))
all.nodes[!(all.nodes %in% unique(unlist(excluded.nodes)))]
# We actually only need to estimate 4 nodes. 

# For each total descendant group of each MRCA you should calculate the ancestral states within those trees
# as normal 

deep.nodes <- all.nodes[!(all.nodes %in% unique(unlist(excluded.nodes)))]

intercept <- anc[getStates(painted.tree, type='tips'),]
rownames(intercept)<-rownames(data); colnames(intercept)<-colnames(data)
# We have a matrix of the optima specific to each animal

data.centered <- data-intercept
# This matrix contains the stationary variance around the optima of the adaptive zones

painted.tree.compressed <- drop.tip.simmap(painted.tree,names(roost)[which(roost=='free')])
painted.tree.free <- drop.tip.simmap(painted.tree,names(roost)[which(roost=='compressed')])
# pruned trees of tight and free roosting bats

tight.cov <- (t(data.centered[painted.tree.compressed$tip,])%*%data.centered[painted.tree.compressed$tip,])/(length(painted.tree.compressed$tip)-1)
free.cov <- (t(data.centered[painted.tree.free$tip,])%*%data.centered[painted.tree.free$tip,])/(length(painted.tree.free$tip)-1)
# Variance-covariance matrices for tight and free roosting bats, estimated as var = (x'x)/(n-1)

OUM <- mvMORPH::mvOU(tree=painted.tree, model='OUM', data=data[,1:4], optimization='subplex',  param = list(vcv="randomRoot", alpha = "diagonal", sigma='diagonal'))
# A likelihood maximized estimate of a multiple OU model, which includes stationary variance as a latent feature


favour <- matrix(NA,100,3)

# Estimate Markov Model and produce Stochastic character maps of Exposed/Enclosed roosting ecology
ard_model<-fitMk(eco.tree,roost,model="ARD",pi='equal')
AIC(ard_model)
set.seed(1)
# Set for reproducibility
mcmc100 <- simmap(ard_model, nsim=100, method='mcmc')
obj <- summary(mcmc100)
crop.mcmc100 <- list()
for(i in 1:length(mcmc100)){
	crop.mcmc100[[i]] <- keep.tip.simmap(phy=mcmc100[[i]],rownames(x)) 
}
class(crop.mcmc100)<-class(mcmc100)
crop.obj <- summary(crop.mcmc100)
# Cropping to taxa for which sternum data are available 

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
# This may, depending on your hardware, take a long time
for(j in 1:100){
	BMM <- tryCatch({mvMORPH::mvBM(tree=crop.mcmc100[[j]], data=data, model="BMM", optimization= "Nelder-Mead")},
	error=function(e){}) 
	if(is.null(BMM)==F){
		favour[j,2] <- BMM$AICc
		rm(BMM)
	}
	print(j)
	
}

root <- matrix(NA,100,1)
# Matrices to store model outcomes, and the inferred root state of roosting ecology 
for(j in 1:100){
	root[j] <- getStates(crop.mcmc100[[j]])[1]
}


library(ggplot2)# v 3.5.2 
library(cowplot)# v 1.1.3
# Plotting packages

# Model AICc outputs, which will be the comparison metric, are going to be stored for expression
# as a histogram: 
histy <- data.frame(cbind(favour,root))
for(i in 1:3){
	histy[,i]<-as.numeric(histy[,i])
}


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

library(extrafont) # 0.19 
# Package for plotting Times New Roman fonts 
# May take time when called

# Make histogram of model performance metrics, identifying early burst (Beta), multiple adaptive optimum (Theta)
# models, as well as subsetting the histogram depending on whether the ancestor of all bats is estimated to have roosted
# in exposed or enclosed settings: 

AIC_histogram<-
ggplot()+
geom_histogram(data=histy,binwidth = 3, aes(X2,fill=X4),position="stack", colour='black')+
viridis::scale_fill_viridis(labels=c('Enclosed','Exposed'), direction=1,option='D',discrete=T)+
 coord_flip()+
labs(x='AICc',y='Counts',fill='Root state')+
geom_segment(aes(x=OUM$AICc,xend=OUM$AICc,y=10,yend=0.3),linewidth=4*Subplot.module.linewidth,arrow = arrow(length=unit(3, "mm")))+
geom_segment(aes(x=OUM$AICc,xend=OUM$AICc,y=10,yend=0.3),linewidth=3*Subplot.module.linewidth/2,arrow = arrow(length=unit(3, "mm")),color='red')+
geom_segment(aes(x=EB$AICc,xend=EB$AICc,y=10,yend=0.3),linewidth=4*Subplot.module.linewidth,arrow = arrow(length=unit(3, "mm")))+
geom_segment(aes(x=EB$AICc,xend=EB$AICc,y=10,yend=0.3),linewidth=3*Subplot.module.linewidth/2,arrow = arrow(length=unit(3, "mm")),color='red')+
geom_text(aes(x=OUM$AICc-1,y=11.5),label=expression(Theta),size=2*Subplot.axis.font.size/.pt, family='Times New Roman')+
geom_text(aes(x=EB$AICc-1,y=11.5),label=expression(beta),size=2*Subplot.axis.font.size/.pt, family='Times New Roman')+
    theme(legend.position="bottom",
	legend.title=element_text(size=Subplot.axis.font.size, family='Times New Roman'),
	legend.text=element_text(size=Legend.axis.font.size, family='Times New Roman'),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		axis.text=element_text(size=Subplot.axis.font.size, family='Times New Roman'),
		axis.ticks=element_blank(),
		axis.title=element_text(size=Subplot.axis.font.size, family='Times New Roman'),
		legend.spacing.x = unit(-.0,'cm'),
		legend.spacing.y = unit(-.0,'cm'),
		legend.key.width=unit(1,'cm'), 
		legend.key.height=unit(.5,'cm'), 
		legend.key.spacing.x=unit(.1,'cm'),
		legend.key.spacing.y=unit(-.0,'cm'),
		axis.line = element_line(color = "black", linewidth = Subplot.module.linewidth))

AIC.legend <-get_plot_component(AIC_histogram, "guide-box-bottom")
AIC_histogram <- AIC_histogram + theme(legend.position = "none")
# Obtain legend and remove it from the main histogram, we will assemble plot components later. 



# We are now going to estimate evolutionary 'Quantum leaps' in sternum shape across the history of bats. 

modded.data <- intercept
modded.data[names(which(getStates(painted.tree, type='tips')=='Horseshoe')),] <- matrix(1,length(names(which(getStates(painted.tree, type='tips')=='Horseshoe'))),1)%*%anc[6,]
modded.data[names(which(getStates(painted.tree, type='tips')=='Phyllostomidae')),] <- matrix(1,length(names(which(getStates(painted.tree, type='tips')=='Phyllostomidae'))),1)%*%anc[4,]
modded.data[names(which(getStates(painted.tree, type='tips')=='Molossidae')),] <- matrix(1,length(names(which(getStates(painted.tree, type='tips')=='Molossidae'))),1)%*%anc[1,]
# Assign estimates of the adaptive optima of shape values of each major branch of bats to their ancestral nodes in the family tree

pruned.tree$node.label <- (length(pruned.tree$tip) + 1):((length(pruned.tree$tip))+ pruned.tree$Nnode)
anc.e<-list()
anc.e <- matrix(NA, length(all.nodes) ,dim(data)[2]) 
rownames(anc.e) <- 1:dim(anc.e)[1]

for(i in MRCAs){
	a <- keep.tip(pruned.tree,pruned.tree$tip[getDescendants(pruned.tree, i)][complete.cases(pruned.tree$tip[getDescendants(pruned.tree, i)])])
	for(k in 1:dim(anc.e)[2]){
		estimates <- ace(data[a$tip,k], a)$ace
		anc.e[as.character(names(estimates)),k] <- estimates
	}
}

color <-names(unlist(lapply(painted.tree$maps, head, n = 1)))

# We are going to propagate a Mahalanobis-distance minimizing optimizer down the family tree, assigning unresolved nodes
# to the trait values of the optimum that a phylogenetic mean estimate is closest to, and then we will iterate until we reach the root, producing
# one possible history of bat sternum shape evolution through 'Quantum Leaps' while minimizing the overall amount of inferred
# evolutionary change that has elapsed:

modded.data <- intercept
for(i in deep.nodes){
	a <- keep.tip(pruned.tree,pruned.tree$tip[getDescendants(pruned.tree, i)][complete.cases(pruned.tree$tip[getDescendants(pruned.tree, i)])])
	# Sub tree
	species <- a$tip
	ones<- matrix(1,length(species),1)
	Csub <- ape::vcv.phylo(keep.tip(pruned.tree,species))
	anc.e[i,] <- (t(ones)%*%solve(Csub)%*%modded.data[species,])/sum(solve(Csub))
	# Sub tree phylogenetic mean 

	tippies <- getDescendants(pruned.tree, i)
	tippies <- tippies[which(tippies<=length(pruned.tree$tip))]
	children <-pruned.tree$edge[,2][which(pruned.tree$edge[,1] == i)]
	# Adjacent lineages

	dist.daughter.1 <- mahalanobis(x=anc, center=anc.e[children[1],], cov= tight.cov)^.5
	dist.daughter.2 <- mahalanobis(x=anc, center=anc.e[children[2],], cov= tight.cov)^.5
	# Distances to adjacent lineages' optima in trait space 

	represented <- c(which(dist.daughter.1==min(dist.daughter.1)),which(dist.daughter.2==min(dist.daughter.2)))
	# We want to restrict to only the optima available in the descendants 

	distances <- mahalanobis(x=anc[represented,], center=anc.e[i,], cov= tight.cov)^.5 # OUM$alpha
	# Which optimum is closest? 

	closest <- which(distances==min(distances))[1]
	anc.e[i,] <- anc[represented[closest],]
	modded.data[tippies,] <- ones%*%anc[represented[closest],]

	color[which(pruned.tree$edge[,2]==children)] <- names(represented[-closest])
	print(which(pruned.tree$edge[,2]==children))
	print(which(pruned.tree$edge[,1] == i))
}

anc.e[deep.nodes,]

anc.e[89,]<- anc['Grade',]
anc.e[91,]<- anc['Grade',]
anc.e[92,]<- anc['Phyllostomidae',]
# Resolving a paraphyletic relationship; the Quantum Leap topology is stated a priori 
# (Phyllostomid type sternum assumed a derivative of the para-phyllostomid grade sternum shape)

nodals <- anc.e
nodals[1:length(pruned.tree$tip),] <- data[pruned.tree$tip,]
rownames(nodals)[1:length(pruned.tree$tip)]<-pruned.tree$tip

N <- length(pruned.tree$tip.label)
edges <- as.matrix(pruned.tree$edge)
# Harvesting various plotting parameters

for (i in 1:nrow(edges)) {
	if(i==1){
		pts <- c(nodals[edges[i, ][1], ],nodals[edges[i, ][2], ])
	} else {
		pts <- rbind(pts,c(nodals[edges[i, ][1], ],nodals[edges[i, ][2], ]))
	}
}
# This conditional loop builds the data frame for the family tree plotting 

dat.pts<-data.frame(as.matrix(pts))
colnames(dat.pts) <- c('PC1','PC2','PC3','PC4','PC1.2','PC2.2','PC3.2','PC4.2')
colnames(nodals)<-c('PC1','PC2','PC3','PC4')
# Various plot metadata 

# Plot the likely topology of evolutionary Quantum leaps in Sternum morphology across bats: 
Leaps<-
ggplot()+
	geom_segment(data=dat.pts,aes(x=PC1,xend=PC1.2,y=PC2,yend=PC2.2), lwd=Subplot.module.linewidth, color=as.numeric(as.factor(color)),alpha=.75)+
	geom_point(data=data.frame(nodals[1:70,1:2]), aes(x=PC1, y=PC2, col= (factor(getStates(painted.tree, type='tips'))), fill= (factor(getStates(painted.tree, type='tips'))),shape= (factor(getStates(painted.tree, type='tips'))) ), size=1 )+
	scale_fill_manual(labels = c('Pteropodidae','Para-phyllostomid','Rhinolophoidea','Molossidae','Phyllostomidae','Vespertilionidae'),
	values=levels(factor(1:6)) )+
	scale_color_manual(labels = c('Pteropodidae','Para-phyllostomid','Rhinolophoidea','Molossidae','Phyllostomidae','Vespertilionidae'),
	values=levels(factor(1:6)) )+
	scale_shape_manual(labels = c('Pteropodidae','Para-phyllostomid','Rhinolophoidea','Molossidae','Phyllostomidae','Vespertilionidae'),
	values=c(8,21:25) )+
	geom_segment(data=dat.pts[c(97,36,134,135,138),],aes(x=PC1,xend=PC1.2,y=PC2,yend=PC2.2), lwd=6*Subplot.module.linewidth, 
	arrow = arrow(length = unit(4, "mm")),color='black')+
	geom_segment(data=dat.pts[c(97),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.5,0.5)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.5,0.5))), lwd=4*Subplot.module.linewidth,color=2,alpha=.65)+
	geom_segment(data=dat.pts[c(97),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.65,0.35)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.65,0.35))), lwd=4*Subplot.module.linewidth,color=2,alpha=.75)+
	geom_segment(data=dat.pts[c(97),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.75,0.25)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.75,0.25))), lwd=4*Subplot.module.linewidth,color=2,alpha=1)+
	geom_segment(data=dat.pts[c(97),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.35,0.65)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.35,0.65)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=5,alpha=.8)+
	geom_segment(data=dat.pts[c(97),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.25,0.75)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.25,0.75)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=5,alpha=1,
	arrow = arrow(length = unit(3.75, "mm")))+

	geom_segment(data=dat.pts[c(36),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.5,0.5)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.5,0.5))), lwd=4*Subplot.module.linewidth,color=3,alpha=.5)+
	geom_segment(data=dat.pts[c(36),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.65,0.35)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.65,0.35))), lwd=4*Subplot.module.linewidth,color=3,alpha=.65)+
	geom_segment(data=dat.pts[c(36),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.75,0.25)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.75,0.25))), lwd=4*Subplot.module.linewidth,color=3,alpha=.75)+
	geom_segment(data=dat.pts[c(36),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.35,0.65)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.35,0.65)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color='white',alpha=.8)+
	geom_segment(data=dat.pts[c(36),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.25,0.75)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.25,0.75)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color='white',alpha=1,
	arrow = arrow(length = unit(3.75, "mm")))+

	geom_segment(data=dat.pts[c(134),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.5,0.5)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.5,0.5))), lwd=4*Subplot.module.linewidth,color=6,alpha=.5)+
	geom_segment(data=dat.pts[c(134),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.65,0.35)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.65,0.35))), lwd=4*Subplot.module.linewidth,color=6,alpha=.65)+
	geom_segment(data=dat.pts[c(134),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.75,0.25)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.75,0.25))), lwd=4*Subplot.module.linewidth,color=6,alpha=.75)+
	geom_segment(data=dat.pts[c(134),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.35,0.65)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.35,0.65)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=4,alpha=.8)+
	geom_segment(data=dat.pts[c(134),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.25,0.75)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.25,0.75)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=4,alpha=1,
	arrow = arrow(length = unit(3.75, "mm")))+

	geom_segment(data=dat.pts[c(135),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.5,0.5)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.5,0.5))), lwd=4*Subplot.module.linewidth,color=2,alpha=.5)+
	geom_segment(data=dat.pts[c(135),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.65,0.35)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.65,0.35))), lwd=4*Subplot.module.linewidth,color=2,alpha=.65)+
	geom_segment(data=dat.pts[c(135),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.75,0.25)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.75,0.25))), lwd=4*Subplot.module.linewidth,color=2,alpha=.75)+
	geom_segment(data=dat.pts[c(135),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.35,0.65)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.35,0.65)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=6,alpha=.8)+
	geom_segment(data=dat.pts[c(135),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.25,0.75)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.25,0.75)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=6,alpha=1,
	arrow = arrow(length = unit(3.75, "mm")))+

	geom_segment(data=dat.pts[c(138),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.5,0.5)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.5,0.5))), lwd=4*Subplot.module.linewidth,color=2,alpha=.5)+
	geom_segment(data=dat.pts[c(138),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.65,0.35)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.65,0.35))), lwd=4*Subplot.module.linewidth,color=2,alpha=.65)+
	geom_segment(data=dat.pts[c(138),],aes(x=PC1,xend=weighted.mean(c(PC1,PC1.2), w=c(0.75,0.25)),y=PC2,yend=weighted.mean(c(PC2,PC2.2), w=c(0.75,0.25))), lwd=4*Subplot.module.linewidth,color=2,alpha=.75)+
	geom_segment(data=dat.pts[c(138),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.35,0.65)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.35,0.65)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=3,alpha=.8)+
	geom_segment(data=dat.pts[c(138),],aes(x=weighted.mean(c(PC1,PC1.2), w=c(0.25,0.75)),xend=PC1.2,y=weighted.mean(c(PC2,PC2.2), w=c(0.25,0.75)),yend=PC2.2), lwd=4*Subplot.module.linewidth,color=3,alpha=1,
	arrow = arrow(length =  unit(3.75, "mm")))+
   	 theme(
	legend.position='bottom',
	legend.justification = c("left", "bottom"),
	legend.title=element_blank(),
	legend.text=element_text(size=Legend.title.font.size, family='Times New Roman'),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		axis.text=element_blank(),
		axis.ticks=element_blank(),
		axis.title=element_text(size=Legend.title.font.size, family='Times New Roman'),
		axis.title.x=element_text(size=Subplot.axis.font.size, family='Times New Roman'),
		axis.ticks.x = element_blank(),
		axis.text.x=element_blank(),
		legend.spacing.x = unit(-.0,'cm'),
		legend.spacing.y = unit(-.0,'cm'),
		legend.key.width=unit(1,'cm'), # 0.2
		legend.key.height=unit(.5,'cm'), #0.2
		legend.key.spacing.x=unit(-.0,'cm'),
		legend.key.spacing.y=unit(-.0,'cm'),
		legend.background=element_rect(fill='transparent',color=NA),
		axis.line = element_line(color = "black", linewidth = Subplot.module.linewidth ))


Leaps.legend <-get_plot_component(Leaps, "guide-box-bottom")
Leaps <- Leaps + theme(legend.position = "none")
# Divide the main plot and legend so we can rearrange them later 


# We also want to have a figure illustrating the leap topology with explicit distances in 4-D trait space labelled and flatenned: 
# We will need to know the positions of the points. 
# We will determine the Y-axis position scaled sympathetically with the Mahalanobis distance under a tight covariance matrix
# The X-axis position should be determined by the node depth 
# Start with position 'Grade'/para-phyllostomids at zero, this is root state. 
# Then positive Y towards Rhinolophoidea
# Then positive Y towards Pteropodidae


node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Hipposideridae' | families=='Rhinolophidae' |families=='Pteropodidae' ))),rownames(crop.obj$ace))
Rhinolophoidea.X<-nodeheight(painted.tree,node+length(painted.tree$tip))
Rhinolophoidea.Y<-(t(as.matrix(anc[5,]-anc[4,]))%*%solve(tight.cov)%*%(as.matrix(anc[5,]-anc[4,])))^.5

node <- match(getMRCA(crop.mcmc100[[1]], names(which(families=='Pteropodidae'))),rownames(crop.obj$ace))
nodeheight(painted.tree,node+length(painted.tree$tip))
Pteropodidae.X<-nodeheight(painted.tree,node+length(painted.tree$tip))
Pteropodidae.Y<-(t(as.matrix(anc[6,]-anc[5,]))%*%solve(tight.cov)%*%(as.matrix(anc[6,]-anc[5,])))^.5 + Rhinolophoidea.Y

node <- match(getMRCA(crop.mcmc100[[1]], c('Pteronotus_davyi',names(which(families=='Phyllostomidae'| families=='Vespertilionidae')) )),rownames(crop.obj$ace))
nodeheight(painted.tree,node+length(painted.tree$tip))
Grade.X<-nodeheight(painted.tree,node+length(painted.tree$tip))
Grade.Y<-0

node <- match(getMRCA(crop.mcmc100[[1]], c('Pteronotus_davyi',names(which(families=='Phyllostomidae') ))),rownames(crop.obj$ace))
nodeheight(painted.tree,node+length(painted.tree$tip))
Phyllostomidae.0<-nodeheight(painted.tree,node+length(painted.tree$tip))
node <- match(getMRCA(crop.mcmc100[[1]], c(names(which(families=='Phyllostomidae') ))),rownames(crop.obj$ace))
Phyllostomidae.X<-nodeheight(painted.tree,node+length(painted.tree$tip))
Phyllostomidae.Y<-(t(as.matrix(anc[3,]-anc[4,]))%*%solve(tight.cov)%*%(as.matrix(anc[3,]-anc[4,])))^.5

roster <-names(which(families=='Molossidae' | families=='Vespertilionidae' ))
roster<-roster[which(roster %in% pruned.tree$tip==T)]
roster<-c(roster,names(which(families=='Vespertilionidae') ))
node <- match(getMRCA(crop.mcmc100[[1]], roster ),rownames(crop.obj$ace))
nodeheight(painted.tree,node+length(painted.tree$tip))
Vespertilionidae.X<-nodeheight(painted.tree,node+length(painted.tree$tip))
Vespertilionidae.Y<-(t(as.matrix(anc[1,]-anc[4,]))%*%solve(tight.cov)%*%(as.matrix(anc[1,]-anc[4,])))^.5

roster <-names(which(families=='Molossidae' ))
roster<-roster[which(roster %in% pruned.tree$tip==T)]
node <- match(getMRCA(crop.mcmc100[[1]], roster ),rownames(crop.obj$ace))
Molossidae.X<-nodeheight(painted.tree,node+length(painted.tree$tip))
Molossidae.Y<-(t(as.matrix(anc[2,]-anc[1,]))%*%solve(tight.cov)%*%(as.matrix(anc[2,]-anc[1,])))^.5 +Vespertilionidae.Y

graph<-cbind(x=c(0,Rhinolophoidea.X,0,Phyllostomidae.0,0,Vespertilionidae.X),
xend=c(Rhinolophoidea.X,Pteropodidae.X,Grade.X,Phyllostomidae.X,Vespertilionidae.X,Molossidae.X),
y=c(0,Rhinolophoidea.Y,0,0,0,-Vespertilionidae.Y),
yend=c(Rhinolophoidea.Y,Pteropodidae.Y,0,Phyllostomidae.Y,-Vespertilionidae.Y,-Molossidae.Y))
graph<-data.frame(graph)


Graph<-
ggplot()+
geom_curve(data=graph[1,],aes(x=x,xend=xend,y=y,yend=yend), 
arrow = arrow(length = unit(0.0, "cm")),
curvature = -0.1, linewidth=3*Subplot.module.linewidth)+
geom_curve(data=graph[2,],aes(x=x,xend=xend,y=y,yend=yend), 
arrow = arrow(length = unit(0.0, "cm")),
curvature = -0.4, linewidth=3*Subplot.module.linewidth)+
geom_curve(data=graph[3,],aes(x=x,xend=xend,y=y,yend=yend), 
curvature = -0.0, linewidth=3*Subplot.module.linewidth)+
geom_curve(data=graph[4,],aes(x=x,xend=xend,y=y,yend=yend), 
arrow = arrow(length = unit(0.0, "cm")),
curvature = -0.3, linewidth=3*Subplot.module.linewidth)+
geom_curve(data=graph[5,],aes(x=x,xend=xend,y=y,yend=yend), 
arrow = arrow(length = unit(0.0, "cm")),
curvature = 0.4, linewidth=3*Subplot.module.linewidth)+
geom_curve(data=graph[6,],aes(x=x,xend=xend,y=y,yend=yend), 
arrow = arrow(length = unit(0.0, "cm")),
curvature = 0.1, linewidth=3*Subplot.module.linewidth)+
geom_curve(aes(x=Grade.X,xend=Phyllostomidae.0,y=0,yend=0), 
arrow = arrow(length = unit(0.0, "cm")),
curvature = 0.0, linewidth=3*Subplot.module.linewidth)+
geom_point(data=graph,aes(x=xend,y=yend),size=6, stroke=1,fill=c(3,1,2,5,6,4),shape=c(22,8,21,24,25,23))+
geom_text(aes(x=2+(graph[4,]$x+graph[4,]$xend)/2,
y=(graph[4,]$y+graph[4,]$yend)/2,),
label= round(log(pchisq(q=(t(as.matrix(anc[3,]-anc[4,]))%*%solve(free.cov)%*%(as.matrix(anc[3,]-anc[4,]))), df=4, ncp = 0, lower.tail = F))-
log(pchisq(q=(t(as.matrix(anc[3,]-anc[4,]))%*%solve(tight.cov)%*%(as.matrix(anc[3,]-anc[4,]))), df=4, ncp = 0, lower.tail = F)),digits=2),size=Subplot.axis.font.size/.pt, family='Times New Roman' )+
geom_text(aes(x=3+(graph[1,]$x+graph[1,]$xend)/2,
y=(graph[1,]$y+graph[1,]$yend)/2,),
label= round(log(pchisq(q=(t(as.matrix(anc[5,]-anc[4,]))%*%solve(free.cov)%*%(as.matrix(anc[5,]-anc[4,]))), df=4, ncp = 0, lower.tail = F))-
log(pchisq(q=(t(as.matrix(anc[5,]-anc[4,]))%*%solve(tight.cov)%*%(as.matrix(anc[5,]-anc[4,]))), df=4, ncp = 0, lower.tail = F)),digits=2),size=Subplot.axis.font.size/.pt, family='Times New Roman' )+
geom_text(aes(x=-1+(graph[2,]$x+graph[2,]$xend)/2,
y=(graph[2,]$y+graph[2,]$yend)/2,),
label= round(log(pchisq(q=(t(as.matrix(anc[5,]-anc[6,]))%*%solve(free.cov)%*%(as.matrix(anc[5,]-anc[6,]))), df=4, ncp = 0, lower.tail = F))-
log(pchisq(q=(t(as.matrix(anc[5,]-anc[6,]))%*%solve(tight.cov)%*%(as.matrix(anc[5,]-anc[6,]))), df=4, ncp = 0, lower.tail = F)),digits=2),size=Subplot.axis.font.size/.pt, family='Times New Roman' )+
geom_text(aes(x=3+(graph[5,]$x+graph[5,]$xend)/2,
y=(graph[5,]$y+graph[5,]$yend)/2,),
label= round(log(pchisq(q=(t(as.matrix(anc[1,]-anc[4,]))%*%solve(free.cov)%*%(as.matrix(anc[1,]-anc[4,]))), df=4, ncp = 0, lower.tail = F))-
log(pchisq(q=(t(as.matrix(anc[1,]-anc[4,]))%*%solve(tight.cov)%*%(as.matrix(anc[1,]-anc[4,]))), df=4, ncp = 0, lower.tail = F)),digits=2),size=Subplot.axis.font.size/.pt, family='Times New Roman' )+
geom_text(aes(x=3+(graph[6,]$x+graph[6,]$xend)/2,
y=(graph[6,]$y+graph[6,]$yend)/2,),
label= round(log(pchisq(q=(t(as.matrix(anc[1,]-anc[2,]))%*%solve(free.cov)%*%(as.matrix(anc[1,]-anc[2,]))), df=4, ncp = 0, lower.tail = F))-
log(pchisq(q=(t(as.matrix(anc[1,]-anc[2,]))%*%solve(tight.cov)%*%(as.matrix(anc[1,]-anc[2,]))), df=4, ncp = 0, lower.tail = F)),digits=2),size=Subplot.axis.font.size/.pt, family='Times New Roman' )+
lims(x=c(-1,30),y=c(-6,12))+
labs(x='Time (My)',y='Leap Distance')+
    theme(legend.position="bottom",
	legend.title=element_blank(),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		axis.text=element_text(size=Subplot.axis.font.size,color='black', family='Times New Roman'),
		axis.ticks=element_line(),
		axis.title=element_text(size=Subplot.axis.font.size, family='Times New Roman'),
		legend.spacing.x = unit(-.0,'cm'),
		legend.spacing.y = unit(-.0,'cm'),
		legend.key.width=unit(1,'cm'), # 0.2
		legend.key.height=unit(.5,'cm'), #0.2
		legend.key.spacing.x=unit(-.0,'cm'),
		legend.key.spacing.y=unit(-.0,'cm'),
		axis.line = element_line(color = "black", linewidth = Subplot.module.linewidth))


library(Cairo) # v 1.7-0
cowplot::set_null_device("cairo")
# Load package and set plot device 

# Prepare a blank plot background 
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

v<-0.12
v2<-0.145
ggdraw()+
draw_plot(blank)+
draw_plot(Leaps.legend,x=1.2/3,y=0.01,width=1/3,height=1-v)+
draw_plot(AIC.legend,x=0.02,y=0.02,width=1/3,height=0.1)+
draw_plot(AIC_histogram,x=0.01,y=v-0.01,width=1/3,height=1-v)+
draw_plot(Leaps,x=1/3,y=v2-0.01,width=1/3,height=1-v2)+
draw_plot(Graph,x=(2/3)-0.01,y=v-0.01,width=1/3,height=1-v)+
draw_text(text=c('a','b','c'),x=c(0.025,(1/3)-0.01,(2/3)+0.01),y=c(rep(.96,3)),
size=12,family = "Times New Roman") # , fontface='italic'
# Assemble the plots 

# Save the output; adjust the filepath as necessary 
ggsave(filename='Whatever_name_you_want.png',width=18,height=11,unit='cm',dpi=600,device=png, type='cairo')



