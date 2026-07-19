# The purpose of this script is to generate Figure 1 from the manuscript Orkney et al., 2026
# 'Declining rates of evolution and limited convergence in bat
# sternum shape can be explained by a ratchet of specialisation'
# The Figure will portray an anatomically labelled 'average bat sternum' that is 
# calculated as the mean of all available sterna. 
# A spectral decomposition of the feature matrix will be visualized, with different bat groups coloured, 
# and the latent axes of variation will be visualised as warps of the average sternum, with changeable
# features colourised. 
# An estimate of the shared ancestry between species is available from Shi & Rabosky 2015 (https://doi.org/10.1111/evo.12681)
# A species-wise binary bag of words is available from Orkney et al., 2025 (https://doi.org/10.1038/s41559-024-02572-9) 
# describing bat ecological characteristics. 

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

library(extrafont) # v 0.20

bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# This phylogeny can be obtained from Shi & Rabosky 2015: https://doi.org/10.1111/evo.12681
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
# Embed as Principal Components

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

# Assigning subclades with different hypothesized means
#plot(painted.tree) # Optional plotting

library( ggplot2 ) # v 3.4.1
library(cowplot) # v 1.1.3
# Plotting function package
library( ggdendro ) # v 0.1.23
library( dendextend ) # v 1.17.1
library(zoo) # v 1.8-12
library(dplyr) # v 1.1.1
library(viridis) # v 0.6.5 
# Packages for plotting and data curation to prepare it for plotting 

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
	dendr.mod$x[tips] <- fit$estim[ match( tree$edge[,1][match(dendr$segments$xend[tips],tree$edge[,2])],rownames(fit	$estim) ),axis ] 
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

BM <- tryCatch({mvMORPH::mvBM(tree=painted.tree, data=data[,1:4], model="BM1", optimization="Nelder-Mead")},
error=function(e){}) 

#data<-pca$x[,1:4]
#OUM <- mvMORPH::mvOU(tree=painted.tree, model='OUM', data=data[,1:4], optimization='subplex',  param = list(vcv="randomRoot", alpha = #"diagonal", sigma='diagonal'),
#control = list(maxit = 1000000))
# Optional Ornstein Uhlenbeck multiple optima alternative 

i<-4
# PC1 <- anc.recon(tree=painted.tree,data=pca$x[,1:i],mod=OUM,axis=1)
# PC2 <- anc.recon(tree=painted.tree,data=pca$x[,1:i],mod=OUM,axis=2)
# Estimate the projection of our node-wise ancestral state imputations upon PC1 and PC2 of sternum shape variety. 


PC1 <- anc.recon(tree=painted.tree,data=pca$x[,1:i],mod=BM,axis=1)
PC2 <- anc.recon(tree=painted.tree,data=pca$x[,1:i],mod=BM,axis=2)
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

x<-PC1
y<-PC2
starts <- cbind(x[complete.cases(x),1],y[complete.cases(y),1:2])
stops <- cbind(x[complete.cases(x),3],y[complete.cases(y),3:4])
df<-data.frame(cbind(starts[,c(1,2,3)],stops[,c(1,2,3)]))
	colnames(df) <- c('x','y','z','xend','yend','zend')
# Building a data frame to draw the painted phylogeny in the embedded space


# colour the tree
color <- names(unlist(lapply(painted.tree$maps, head, n = 1)))

	dendr <- ggdendro::dendro_data(stats::as.dendrogram(painted.tree))
	lab.dat <- dendr$labels
	dendr.mod<-dendr$segments/2 # This is the tree structure 
	dendr.mod$z <- rep(NA,dim(dendr.mod)[1])
	tips <- which(dendr$segments$yend==0)
	nodes <- which(dendr$segments$y==dendr$segments$yend)
	internal.branches <- which(dendr$segments$x==dendr$segments$xend & dendr$segments$yend!=0)

	
	match.nodes<- match(dendr.mod[nodes ,]$y,max(nodeHeights(painted.tree))-nodeHeights(painted.tree)[,1] )
	dendr.mod[nodes[seq(1,length(nodes),2)] ,]$z <- color[match.nodes][1:((length(color))/2)]


	dendr.mod[nodes ,]$z <- color[match.nodes]
	# This colors the nodes 
	dendr.mod[internal.branches,]$z <- dendr.mod[nodes,]$z[match(dendr.mod[internal.branches,]$yend,dendr.mod[nodes,]$yend)]
	# This colors the stems of clades
	dendr.mod$z[which(is.na(dendr.mod$z)==T)]<-dendr.mod$z[which(is.na(dendr.mod$z)==T)-1]
	# This colors downstream nodes

	dendr.mod$z[c(1,2,3,4,5,6,7,8,9,11,205,207)]<-NA
	# Practice to determine how I change the color of the rooty parts of the tree.  




branch.splatter <- dendr.mod$z[complete.cases(x)]
xsplain <-round(((pca$sdev^2)/sum(pca$sdev^2))[1],digits=2)*100
ysplain <-round(((pca$sdev^2)/sum(pca$sdev^2))[2],digits=2)*100
# Explained variances in the embedded space

# Assemble plot of embedding of phylogeny along first two principal components;
# the aim here is just to recognize that different bat subclades are likely to conform to distinct 
# optima
branchies<-
ggplot()+
geom_segment(data=df, aes(x=x,xend=xend,y=y,yend=yend), lwd=2, col='black',lineend='round')+
geom_segment(data=df, aes(x=x,xend=xend,y=y,yend=yend), lwd=1, col='grey',lineend='round')+
geom_segment(data=df, aes(x=x,xend=xend,y=y,yend=yend, col= factor(branch.splatter)), lwd=1,lineend='round')+
	scale_color_manual(labels = c('Pteropodidae','Para-phyllostomid','Rhinolophoidea','Molossidae','Phyllostomidae','Vespertilionidae'),
	values=levels(factor(1:6)),na.translate=F )+
labs(x=paste('PC1 ',xsplain,' % explained variance',sep=''),y=paste('PC1 ',ysplain,' % explained variance',sep=''))+
	theme(legend.position='bottom',
	legend.text=element_text(size=7.5,family='Times New Roman'),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.key.spacing.x = unit(0.0, 'cm'),
		legend.key.spacing.y = unit(0, 'cm'),
		legend.title=element_blank(),
		axis.title.x=element_text(size=9.5,family='Times New Roman'),
		axis.title.y=element_text(size=9.5,family='Times New Roman'),
		axis.line.x=element_line(linewidth = 1,color='black'),
		axis.line.y=element_line(linewidth = 1,color='black'),
		legend.background = element_rect(fill = 'transparent', color = NA))




# Now we will visualize the changes in the shape of the sternum along the 
library(Morpho) # v 2.13
# Package for handling meshes
mesh <- ply2mesh('sternum.surf.ply')
# A random sample mesh topology
library(rgl) # v 1.3.36
# 3D viewing package

gpca<- gm.prcomp(coords)
preds <- shape.predictor(coords, x = gpca$x[,1:4], Intercept = FALSE, 
					preds0 =c(0,0,0,0),
					pred1 = c(min(gpca$x[,1]),0,0,0), 
					pred2 = c(max(gpca$x[,1]),0,0,0),
					pred3 = c(0,min(gpca$x[,2]),0,0 ),
					pred4 = c(0,max(gpca$x[,2]),0,0 ),
					pred5 = c(0,0,min(gpca$x[,3]),0 ),
					pred6 = c(0,0,max(gpca$x[,3]),0 ),
					pred7 = c(0,0,0,min(gpca$x[,3]) ),
					pred8 = c(0,0,0,max(gpca$x[,3]) ) )
# Projections of sternum landmarks along the minima and maxima of the PC axes, and the null point

# The following lines warp the sample mesh topology to the estimated
# landmark constellations: 
start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[1]]
mean_object <- tps3d(mesh, start_points, end_points)


start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[2]]
warped_PC1_low <- tps3d(mesh, start_points, end_points)
start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[3]]
warped_PC1_high <- tps3d(mesh, start_points, end_points)


start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[4]]
warped_PC2_low <- tps3d(mesh, start_points, end_points)
start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[5]]
warped_PC2_high <- tps3d(mesh, start_points, end_points)

start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[6]]
warped_PC3_low <- tps3d(mesh, start_points, end_points)
start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[7]]
warped_PC3_high <- tps3d(mesh, start_points, end_points)

start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[8]]
warped_PC4_low <- tps3d(mesh, start_points, end_points)
start_points <- sternum.array[,,grep('Saccop',dimnames(sternum.array)[[3]])]
end_points <- preds[[9]]
warped_PC4_high <- tps3d(mesh, start_points, end_points)

# Generate visualizations
PC1_HIGHWARP <- meshDist(warped_PC1_high,warped_PC1_low, sign = FALSE, rampcolors=c('white','red','black'))
PC1_LOWWARP <- meshDist(warped_PC1_low,warped_PC1_high, sign = FALSE, rampcolors=c('white','red','black'))
PC2_HIGHWARP  <- meshDist(warped_PC2_high,warped_PC2_low, sign = FALSE, rampcolors=c('white','red','black'))
PC2_LOWWARP <- meshDist(warped_PC2_low,warped_PC2_high, sign = FALSE, rampcolors=c('white','red','black'))

mesh_distances <- meshDist(warped_PC3_high,warped_PC3_low, sign = FALSE, rampcolors=c('white','red','black'))
mesh_distances <- meshDist(warped_PC3_low,warped_PC3_high, sign = FALSE, rampcolors=c('white','red','black'))
mesh_distances <- meshDist(warped_PC4_high,warped_PC4_low, sign = FALSE, rampcolors=c('white','red','black'))
mesh_distances <- meshDist(warped_PC4_low,warped_PC4_high, sign = FALSE, rampcolors=c('white','red','black'))

# Take screenshots:
# (User may need to open and size a window appropriately before hand)

rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC1_HIGHWARP$colMesh,specular='black')
snapshot3d('PC1_HIGHWARP.png', fmt = "png") #,width = 900, height = 400)
rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC1_LOWWARP$colMesh,specular='black')
snapshot3d('PC1_LOWWARP.png', fmt = "png") #,width = 900, height = 400)

rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC1_HIGHWARP$colMesh,specular='black')
snapshot3d('PC1_HIGHWARP_lateral.png', fmt = "png") #,width = 900, height = 400)
rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC1_LOWWARP$colMesh,specular='black')
snapshot3d('PC1_LOWWARP_lateral.png', fmt = "png") #,width = 900, height = 400)


rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC2_HIGHWARP$colMesh,specular='black')
snapshot3d('PC2_HIGHWARP.png', fmt = "png") #,width = 900, height = 400)
rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC2_LOWWARP$colMesh,specular='black')
snapshot3d('PC2_LOWWARP.png', fmt = "png") #,width = 900, height = 400)

rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC2_HIGHWARP$colMesh,specular='black')
snapshot3d('PC2_HIGHWARP_lateral.png', fmt = "png") #,width = 900, height = 400)
rgl::open3d(windowRect = c(0, 0, 900, 500));rgl::shade3d(PC2_LOWWARP$colMesh,specular='black')
snapshot3d('PC2_LOWWARP_lateral.png', fmt = "png") #,width = 900, height = 400)


rgl::open3d(windowRect = c(0, 0, 3*900, 3*500));rgl::shade3d(mean_object,specular='black',color='grey')
snapshot3d('Mean_sternum.png', fmt = "png") #,width = 900, height = 400)

snapshot3d('Mean_sternum2.png', fmt = "png") #,width = 900, height = 400)


library(png) # v 0.1-9
library(grid) # base package R v 4.60

# We will define a function to prepare our images of 3D meshes. 
# Users replicating our analysis can skip this step 
prep.image <- function(string){
	string <- readPNG(string)
	white <- which(string[,,1]==1 & string[,,2]==1 & string[,,3]==1)
	dim<-dim(string); dim[3] <-4
	temporary <- array(NA, dim = dim)
	temporary[,,1]<-string[,,1]; temporary[,,2]<-string[,,2]; temporary[,,3]<-string[,,3]
	string<-temporary
	string[,,4]<-matrix(1,dim(string[,,1])[1],dim(string[,,1])[2])
	string[,,4][white]<-0
	string.raster <- rasterGrob(string, interpolate = TRUE)
	return(string.raster)
}

# Plot a large image of the 'mean' estimated sternum shape
big_mean <-
ggplot()+
annotation_custom(prep.image('Mean_sternum2.png'), xmin = -	Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
	theme(legend.position='bottom',
	legend.title=element_blank(),
	legend.text=element_text(size=9.5,family='Times New Roman'),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.key.spacing.x = unit(0.5, 'cm'),
		legend.key.spacing.y = unit(0, 'cm'),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.line.y.left=element_blank(),
		axis.ticks.y.left=element_blank(),
		legend.background = element_rect(fill = 'transparent', color = NA))


# Plot images of the different minima and maxima of the embedded space, 
# colourized to indicate which features are most changeable 
warps<-
ggplot()+
annotation_custom(prep.image('PC1_HIGHWARP.png'), xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1)+
theme(panel.background=element_rect(fill=NA,color=NA))+
annotation_custom(prep.image('PC1_LOWWARP.png'), xmin = 0, xmax = 0.25, ymin = 0, ymax = 0.5)+
theme(panel.background=element_rect(fill=NA,color=NA))+
annotation_custom(prep.image('PC1_HIGHWARP_lateral.png'), xmin = 0.25, xmax = 0.5, ymin = 0.5, ymax = 1)+
theme(panel.background=element_rect(fill=NA,color=NA))+
annotation_custom(prep.image('PC1_LOWWARP_lateral.png'), xmin = 0.25, xmax = 0.5, ymin = 0, ymax = 0.5)+
theme(panel.background=element_rect(fill=NA,color=NA))+
annotation_custom(prep.image('PC2_HIGHWARP.png'), xmin = 0.5, xmax = 0.75, ymin = 0.5, ymax = 1)+
theme(panel.background=element_rect(fill=NA,color=NA))+
annotation_custom(prep.image('PC2_LOWWARP.png'), xmin = 0.5, xmax = 0.75, ymin = 0, ymax = 0.5)+
theme(panel.background=element_rect(fill=NA,color=NA))+
annotation_custom(prep.image('PC2_HIGHWARP_lateral.png'), xmin = 0.75, xmax = 1, ymin = 0.5, ymax = 1)+
theme(panel.background=element_rect(fill=NA,color=NA))+
annotation_custom(prep.image('PC2_LOWWARP_lateral.png'), xmin = 0.75, xmax = 1, ymin = 0, ymax = 0.5)+
theme(panel.background=element_rect(fill=NA,color=NA))+
geom_rect(aes(xmin=0,xmax=0.5,ymin=0.5,ymax=1),fill=NA,lwd=1/2,color='darkgrey',lineend='round')+
geom_rect(aes(xmin=0.5,xmax=1,ymin=0.0,ymax=0.5),fill=NA,lwd=1/2,color='darkgrey',lineend='round')+
geom_rect(aes(xmin=0.5,xmax=1,ymin=0.5,ymax=1),fill=NA,lwd=1/2,color='darkgrey',lineend='round')+
geom_rect(aes(xmin=0,xmax=0.5,ymin=0.0,ymax=0.5),fill=NA,lwd=1/2,color='darkgrey',lineend='round')+
geom_text(aes(x=c(0.25,0.25,0.75,0.75),y=c(0.6,0.1,0.6,0.1)),angle=0,label=c('PC1+','PC1-','PC2+','PC2-'),family='Times New Roman')+
	theme(legend.position='bottom',
	legend.title=element_blank(),
	legend.text=element_text(size=9.5,family='Times New Roman'),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.key.spacing.x = unit(0.5, 'cm'),
		legend.key.spacing.y = unit(0, 'cm'),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.line.y.left=element_blank(),
		axis.ticks.y.left=element_blank(),
		legend.background = element_rect(fill = 'transparent', color = NA))

# Load the fonts that the journal requested (takes time)
font_import() 
loadfonts(device = "win") 

library(cowplot) # v 1.2.0
# For combining plots

# Assemble subplots:
row2<- plot_grid(branchies,warps,ncol=2)

dev.new()
plot_grid(big_mean, row2, ncol = 1)+ 
  theme(plot.background = element_rect(fill = "white", colour = NA))+
draw_text(text=c('a','b','c'),x=c(0.05,0.05,0.55),y=c(0.975,0.55,0.55),
size=16.5, family='Times New Roman')+
ggpubr::geom_bracket(
    xmin = c(0.7), xmax = c(0.95),
    y.position = c(0.935), label = c("Manubrium"),
    tip.length = 0.01, family='Times New Roman',label.size=12.5/.pt
  )+
ggpubr::geom_bracket(
    xmin = c(0.26), xmax = c(0.68),
    y.position = c(0.935), label = c("Mesosternum"),
    tip.length = 0.01, family='Times New Roman',label.size=12.5/.pt
  )+
ggpubr::geom_bracket(
    xmin = c(0.02), xmax = c(0.24),
    y.position = c(0.935), label = c("Xiphoid"),
    tip.length = 0.01, family='Times New Roman',label.size=12.5/.pt
  )+
geom_curve(curvature = 0.25,aes(x=0.78,xend=0.7,y=0.82,yend=0.9),color='black')+
geom_curve(curvature = -0.25,aes(x=0.85,xend=0.69,y=0.88,yend=0.8),color='red',size=1)+
geom_text(aes(x=0.7,y=0.91),label=expression(italic('M.pectoralis')), size=9.5/.pt,family='Times New Roman')+
geom_curve(curvature = -0.25,aes(x=0.5,xend=0.7,y=0.82,yend=0.9),color='black')+
geom_curve(curvature = -0.15,aes(x=0.3,xend=0.7,y=0.8,yend=0.9),color='black')+
geom_curve(curvature = 0.15,aes(x=0.67,xend=0.28,y=0.78,yend=0.79),color='red',size=1)+
geom_curve(curvature = -0.15,aes(x=0.05,xend=0.28,y=0.71,yend=0.73),color='red',size=1)+
geom_curve(curvature = 0.25,aes(x=0.1,xend=0.2,y=0.72,yend=0.65),color='black')+
geom_text(aes(x=0.2,y=0.64),label=expression(italic('Fasc. abdominalis')), size=9.5/.pt,family='Times New Roman')+
geom_curve(curvature = 0.15,aes(x=0.935,xend=0.93,y=0.88,yend=0.75),color='red',size=1,linetype='dashed')+
geom_curve(curvature = 0.25,aes(x=0.92,xend=0.8,y=0.82,yend=0.65),color='black')+
geom_text(aes(x=0.75,y=0.64),label=expression(italic('M.clavodeltoideus')), size=9.5/.pt,family='Times New Roman')+
geom_curve(curvature = -0.10,aes(x=0.3,xend=0.67,y=0.72,yend=0.73),color='red',size=1,linetype='dashed')+
geom_curve(curvature = 0.25,aes(x=0.5,xend=0.45,y=0.74,yend=0.65),color='black')+
geom_text(aes(x=0.45,y=0.64),label=expression(italic('Art. costosterna')), size=9.5/.pt,family='Times New Roman')+
geom_curve(curvature = -0.15,aes(x=0.92,xend=0.75,y=0.6,yend=0.58),color='black')+
geom_text(aes(x=0.67,y=0.58),label=expression(italic('Art. costochondral')), size=9.5/.pt,family='Times New Roman')+
geom_curve(curvature = -0.30,aes(x=0.97,xend=0.76,y=0.62,yend=0.54),color='black')+
geom_text(aes(x=0.67,y=0.54),label=expression(italic('Art. sternoclavicular')), size=9.5/.pt,family='Times New Roman')


# Save plot 
ggsave(filename='Figure_1.png',width=18,height=18,unit='cm',dpi=600,device=png, type='cairo')

