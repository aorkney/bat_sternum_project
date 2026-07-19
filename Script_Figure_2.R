# The purpose of this script is to produce Figure 1 in the Main manuscript of Orkney et al., 2026. 
# 'Declining rates of evolution and limited convergence in bat
# sternum shape can be explained by a ratchet of specialisation'
# The figure (Figure 2) will display a family tree of the bat species studied. 
# Distinct groups of bats that practice different combinations of 
# flight-style and possess characteristically different sternal morphology will be identified on the tree. 
# Representative 3D images of sterna and our geometric measurement scheme used to quantify sternal shape will be presented in
# conjunction with this.
# An estimate of the shared ancestry between species is available from Shi & Rabosky 2015 https://doi.org/10.1111/evo.12681 

# This script was written by Dr. A Orkney and the final version was compiled on July 19th 2026. 
# Users aiming to replicate this figure should observe that the 3D mesh files required to render 3D images are not supplied 
# for public use because of file size constraints. Users should disable these components of the code. 


setwd()
# Set the directory to the location of the landmark data
# This will change if you download the data to a personal computer

load('sternum_array_sep_27_2024.RData')
# load the sternum shape data data

taxa <- dimnames(sternum.array)[[3]]
# This is a vector of all the bat species that Elizabeth and Beyonca have landmarked

metadata <- read.csv('Bat_CT_process_list_Andrew_only.csv')
families <- metadata$Family[ match(taxa,metadata$Shi_match) ]
names(families) <- taxa
original.names <- paste(metadata$Genus[match(taxa,metadata$Shi_match)], metadata$Species[match(taxa,metadata$Shi_match)],sep='_')
names(original.names)<-taxa
# Substitute names in our collected taxa with congeners on the phylogeny of Shi & Rabosky 2015, if 
# a direct match is not available. 


library(ape) # v 5.8-1
# Package for managing family trees

bat.tree <- read.tree('chiroptera.no_outgroups.absolute.tre')
# Phylogeny of Shi & Rabosky 2015: https://doi.org/10.1111/evo.12681
# This phylogeny far exceeds the number of taxa for which we have landmark constellations
# we must therefore prune the phylogeny to the taxa of interest.

pruned.tree <- keep.tip(bat.tree,dimnames(sternum.array)[[3]])
# prune the bat tree to the taxa of interest

pruned.tree <- drop.tip(pruned.tree, c('Micropteropus_pusillus','Molossops_temminckii'))
# I don't have ecological metadata for Micropteropus pussilus
# It is a nectar and fruit loving bat, that seems to prefer tropical rainforest
# margins in subsaharan Africa, but I did not find a clear description of
# its prefered roosts. 

library(phytools) # v 2.4-4
# Package for managing phylogenetic trees

pruned.tree <- phytools::force.ultrametric(pruned.tree)
# Render tree ultrametric 

library(geomorph) # v 4.0.10
# Package to manage shape data 

sliders<- rbind( geomorph::define.sliders(c(7, 17:20,9), write.file=F) ,
                 geomorph::define.sliders(c(8, 21:24,10), write.file=F) ,
                 geomorph::define.sliders(c(3, 25:31,4), write.file=F) ,
                 geomorph::define.sliders(c(1, 32:39,2), write.file=F) ,
                 geomorph::define.sliders(c(15, 40:48,13), write.file=F) ,
                 geomorph::define.sliders(c(16, 49:57,14), write.file=F) ,
                 geomorph::define.sliders(c(15, 58:67,12), write.file=F))
# These landmarks, which represent curves, will be allowed to slide along their tangent vectors
# to achieve a superior alignment of shape data prior to analysis.

GPA.sternum<-geomorph::gpagen(sternum.array, curves=sliders, approxBE=T )
# Align the constellations of landmarks into a common reference frame. 

coords <- GPA.sternum$coords[,,pruned.tree$tip]
# These are the aligned coordinates

consensus.shape <- GPA.sternum$consensus
col<-rep('red',length(consensus.shape[,1]))
col[unlist(sliders[,2])] <-'yellow'
# Landmark constellation used to record sternum shape 

# The following is a function which retrieves 3D mesh files to prepare them to make images 
# Morpho 2.12 is a dependency 
# Users replicating our analysis do not need to process meshes and may skip past this section

prepare.mesh <- function(species){
	setwd('H:/Andrew_backup/Avizo_projects/')
	path <- dir()[grep(species, dir())]
	setwd(paste(getwd(),path,sep='/'))
	mesh <- Morpho::ply2mesh("Sternum.surf.ply")
	target <- GPA.sternum$coords[,,gsub(' ','_',species)][-unlist(sliders[,2]),]
	scaled.points <- sternum.array[,,gsub(' ','_',species)]/GPA.sternum$Csize[gsub(' ','_',species)]
	scaled.points.save <- scaled.points
	# We have to center the data
	scaled.points[,1] <- scaled.points[,1]-mean(scaled.points[,1])
	scaled.points[,2] <- scaled.points[,2]-mean(scaled.points[,2])
	scaled.points[,3] <- scaled.points[,3]-mean(scaled.points[,3])

	H <-  t(scaled.points[-unlist(sliders[,2]),])%*%target
	svd <- svd(H)
	R <- svd$v%*%diag(1,3)%*%t(svd$u)	
	# Find the rotation matrix that aligns the specimen and target

	point.cloud <- mesh$vb[1:3,]/GPA.sternum$Csize[gsub(' ','_',species)]
	# we need to center the point.cloud
	point.cloud[1,] <- point.cloud[1,] -mean(scaled.points.save[,1])
	point.cloud[2,] <- point.cloud[2,] -mean(scaled.points.save[,2])
	point.cloud[3,] <- point.cloud[3,] -mean(scaled.points.save[,3])

	scaled.points<- t(R%*%t(scaled.points))%*%diag(1.01,3)
	mesh$vb[1:3,] <- (R%*%(point.cloud))
	# Perform rotation

	point.cloud <- mesh$normals[1:3,]
	# we probably don't need to center the point.cloud
	#point.cloud[1,] <- point.cloud[1,] -mean(scaled.points.save[,1])
	#point.cloud[2,] <- point.cloud[2,] -mean(scaled.points.save[,2])
	#point.cloud[3,] <- point.cloud[3,] -mean(scaled.points.save[,3])
	mesh$normals[1:3,] <- (R%*%(point.cloud))
	return(list(mesh,scaled.points))
}

# We're going to take screenshots of our annotated mesh files
# Readers emulating our analysis will not have these mesh files, 
# so are advised to skip these operations. 

library(rgl) # v 1.3.17
open3d(windowRect = c(0, 0, 900, 500))
mesh <- prepare.mesh('Aethalops alecto')
library(Morpho); library(rgl)
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
#clear3d(type = "lights")
#light3d(theta = 50, phi = 250, diffuse = "gray90", specular = "white") # 45
#light3d(theta = 0, phi = -10, diffuse = "gray90", specular = "white")
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
theta=0;phi=-25
view3d(theta, phi,zoom = .48)
snapshot3d('Aethalops_alecto_08_21.png', fmt = "png") #,width = 900, height = 400)

open3d(windowRect = c(0, 0, 900, 500))
mesh <- prepare.mesh('Hipposideros pratti')
library(Morpho); library(rgl)
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
#clear3d(type = "lights")
#light3d(theta = 50, phi = 40, diffuse = "gray90", specular = "white") # 45
#light3d(theta = 0, phi = -10, diffuse = "gray90", specular = "white")
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
theta=0;phi=-25
view3d(theta, phi,zoom = .52)
snapshot3d('Hipposideros_pratti_08_21.png', fmt = "png") 

open3d(windowRect = c(0, 0, 900, 500))
mesh <- prepare.mesh('Pteronotus davyi')
library(Morpho); library(rgl)
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
#clear3d(type = "lights")
#light3d(theta = 50, phi = 00, diffuse = "gray90", specular = "white") # 45
#light3d(theta = 0, phi = -10, diffuse = "gray90", specular = "white")
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
theta=0;phi=-25
view3d(theta, phi,zoom = .48)
snapshot3d('Pteronotus_davyi_08_21.png', fmt = "png") 

open3d(windowRect = c(0, 0, 900, 500))
mesh <- prepare.mesh('Carollia sowelli')
library(Morpho); library(rgl)
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
#clear3d(type = "lights")
#light3d(theta = 50, phi = 90, diffuse = "gray90", specular = "white") # 45
#light3d(theta = 0, phi = -10, diffuse = "gray90", specular = "white")
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
theta=0;phi=-25
view3d(theta, phi,zoom = .48)
snapshot3d('Carollia sowelli_08_21.png', fmt = "png") 

open3d(windowRect = c(0, 0, 900, 500))
mesh <- prepare.mesh('Eumops auripendulus')
library(Morpho); library(rgl)
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
#clear3d(type = "lights")
#light3d(theta = 50, phi = 00, diffuse = "gray90", specular = "white") # 45
#light3d(theta = 0, phi = -10, diffuse = "gray90", specular = "white")
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
theta=0;phi=-25
view3d(theta, phi,zoom = .50)
snapshot3d('Eumops_abrasus_08_21.png', fmt = "png") 

open3d(windowRect = c(0, 0, 900, 500))
mesh <- prepare.mesh('Myotis horsfieldii')
library(Morpho); library(rgl)
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
#clear3d(type = "lights")
#light3d(theta = 50, phi = 00, diffuse = "gray90", specular = "white") # 45
#light3d(theta = 0, phi = -10, diffuse = "gray90", specular = "white")
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
theta=0;phi=-25
view3d(theta, phi,zoom = .50)
snapshot3d('Myotis horsfieldii_08_21.png', fmt = "png") 

mesh <- prepare.mesh('Myotis horsfieldii')
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
open3d(windowRect = c(0, 0, 900, 500))
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
theta=0;phi=-25
view3d(theta, phi,zoom = .52)
snapshot3d('Myotis horsfieldii_08_21.png', fmt = "png") 
close3d()
open3d(windowRect = c(0, 0, 900, 500))
#clear3d(type = "lights")
#light3d(x = 180, y = 100, z = 100, ambient = "black")
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
theta=0;phi=-90
view3d(theta, phi,zoom = .52)
snapshot3d('Myotis horsfieldii lateral_08_21.png', fmt = "png") 

open3d(windowRect = c(0, 0, 900, 500))
mesh <- prepare.mesh('Tylonycteris pachypus')
library(Morpho); library(rgl)
shade3d(mesh[[1]], col = "gray90")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
setwd('C:/Users/HedrickLab/Documents/Andy_CUTS/Documents/AOrkney/Sternum_paper_July_2025')
theta=0;phi=-25
view3d(theta, phi,zoom = .52)
snapshot3d('Tylonycteris pachypus_08_21.png', fmt = "png") 
close3d()
open3d(windowRect = c(0, 0, 900, 500))
shade3d(mesh[[1]], col = "white")
points3d(mesh[[2]], col = col, size = 10)
points3d(mesh[[2]], col = 'black', size = 15)
points3d(mesh[[2]], col = col, size = 10)
theta=0;phi=-90
view3d(theta, phi,zoom = .52)
snapshot3d('Tylonycteris pachypus lateral_08_21.png', fmt = "png") 

# We have made all the lovely screenshots that we need


library(png)# v 0.1-8
library(ggplot2) # v 3.5.2
library(cowplot) # v 1.1.3
# Various packages with plotting functionalities 

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

library(grid) # base R 4.60

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


library( ggdendro ) # v 0.1.23
library( dendextend ) # v 1.17.1
library(zoo) # v 1.8-12
library(dplyr) # v 1.1.1
# Various packages for plotting phylogenetic trees in the ggplot2 environment 

# Let's call our lovely labelled images.
# Readers emulating our analysis can skip these operations. 

Aethal <- 
ggplot()+
annotation_custom(prep.image('Aethalops_alecto_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))

Hipposideros <- 
ggplot()+
annotation_custom(prep.image('Hipposideros_pratti_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))

Pteronotus <- 
ggplot()+
annotation_custom(prep.image('Pteronotus_davyi_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))

Carollia <- 
ggplot()+
annotation_custom(prep.image('Carollia sowelli_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))

Eumops <- 
ggplot()+
annotation_custom(prep.image('Eumops_abrasus_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))

Myotis <- 
ggplot()+
annotation_custom(prep.image('Myotis horsfieldii_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))

Myotis_lateral <- 
ggplot()+
annotation_custom(prep.image('Myotis horsfieldii lateral_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))

Tylonycteris_lateral <- 
ggplot()+
annotation_custom(prep.image('Tylonycteris pachypus lateral_08_21.png'), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
theme(panel.background=element_rect(fill=NA,color=NA))


# We want to paint a phylogeny to identify distinct groups of bats which share 
# similar sternum properties and flight-style ecologies.
# We hypothesise that these represent distinct lineages belonging to discrete adaptive zones. 

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
# These lines have annotated the tree



	# The following lines will render the phylogenetic tree in a format that we can plot in a ggplot2 environment 
	dendr <- ggdendro::dendro_data(stats::as.dendrogram(painted.tree))
	lab.dat <- dendr$labels
	dendr.mod<-dendr$segments/2 # This is the tree structure 
	dendr.mod$z <- rep(NA,dim(dendr.mod)[1])
	tips <- which(dendr$segments$yend==0)
	nodes <- which(dendr$segments$y==dendr$segments$yend)
	internal.branches <- which(dendr$segments$x==dendr$segments$xend & dendr$segments$yend!=0)
	
	color <-names(unlist(lapply(painted.tree$maps, head, n = 1)))
	
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


# We are going to annotate the family tree with a tile plot illustrating shape properties
# as the projection of shapes upon the leading eigenvectors of a spectral decomposition of the shape matrix

x<- geomorph::two.d.array(coords)
# It will be convenient to treat the coordinates as a 2-D array
pca<-prcomp(x)
i<-2
data <- pca$x[,1:i]
# Prepare dataframe of PC scores 

tile.data <- reshape2::melt(data) # reshape2 1.4.4
levels(tile.data$Var1)<- lab.dat$label
# Organise sternum shape morphology leading variables 

metadata <- read.csv('Bat_eco_metadata.csv')
# Load ecological metadata to produce a stacked barchart of bat flight style ecology across the family tree 

eco <- metadata
# We have re-ordered the ecological data and subsetted it to only those bats that have landmarked sternum bones

binary.flight.scores <- eco[,c(21:35)]
binary.flight.scores<-as.matrix(binary.flight.scores)
# Prepare the data as a matrix to make it easy to index. 
binary.flight.scores[which(binary.flight.scores=='?' | binary.flight.scores=='' )] <- 0 
binary.flight.scores[grep(binary.flight.scores,pattern='\\?')] <- 0 
binary.flight.scores <- apply(binary.flight.scores,2,FUN=as.numeric)
# Set the row names of the flight style matrix to our taxa 
binary.flight.scores <- binary.flight.scores[,-c(13,15)]
# The flight style scores for bats have been compiled. 1 indicates a bat is recorded as practising the flight style. 

library(extrafont) # v 0.19
# Package for Times New Roman font

df2 <- data.frame(binary.flight.scores)
# coerce matric to dataframe in preparation for plotting 
df3<-cbind(eco$Shi,df2)
colnames(df3)[1]<-'taxon'
df3<- df3[match(lab.dat$label,df3$taxon),]
df3<-reshape2::melt(df3)
df3<-df3[-which(df3$value=='0'),]
df3$variable <- factor(df3$variable, levels= c(names(table(df3$variable))[order(table(df3$variable))]))
df3$taxon <- factor(df3$taxon, levels= lab.dat$label)
# Dataframe of flight style properties prepared 

# Build a stacked bar chart of flight style ecology and store it:
check <- ggplot()+
geom_bar(position = "stack",aes( x=taxon, fill=variable ),data=df3)+
 	scale_fill_manual(values= c("gray30", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#0000FF","grey","#FFFF00","#000000","#FF0000"),
	labels = (names(table(df3$variable)))  )+
	labs(fill='')+
 	guides(fill = guide_legend(nrow = 3))+
	theme(legend.title=element_text(size=9.5,family='Times New Roman'),
	legend.text=element_text(size=9.5,family='Times New Roman'))


manip <- ggplot_build(check)$data[[1]]
manip$x <- manip$x +12.5
manip$ymin <- manip$ymin +2.6
manip$ymax <- manip$ymax +2.6
# Extract the plot data so that we can manipulate it 

library(Cairo) #1.7-0
# Package for the aesthetic plotting export PNG device 
cowplot::set_null_device("cairo")


fill.legend <-  get_legend(check)
# Extract the legend of the stacked bar chart for later use 


ticks <- data.frame(x=rep(-1,7),xend=rep(0,7),y=0.35-seq(0,60,by=10)/7,yend=0.35-seq(0,60,by=10)/7)
# A tiny data frame of 10 million year interval marks

# Now it is time to assemble the whole plot 
main.block<-
	ggplot()+
	geom_tile(data=tile.data, aes(y=Var2, x=Var1, fill = value))+
	viridis::scale_fill_viridis(option="A", limits=c(range(c(tile.data$value))))+
	geom_segment(data=dendr.mod, aes(y = .35-y/7, x = (2+dim(dendr.mod)[1]/4)-x*2, yend = .35-yend/7, xend =(2+dim(dendr.mod)[1]/4)-xend*2), lwd=8*Subplot.module.linewidth,col='black',lineend = "round" )+
	geom_segment(data=dendr.mod, aes(y = .35-y/7, x = (2+dim(dendr.mod)[1]/4)-x*2, yend = .35-yend/7, xend =(2+dim(dendr.mod)[1]/4)-xend*2), lwd=6*Subplot.module.linewidth/2,col='grey',lineend = "round" )+
	geom_segment(data=dendr.mod, aes(y = .35-y/7, x = (2+dim(dendr.mod)[1]/4)-x*2, yend = .35-yend/7, xend =(2+dim(dendr.mod)[1]/4)-xend*2, col= z), lwd=6*Subplot.module.linewidth/2,lineend = "round" )+
	scale_color_manual(labels = c('Pteropodidae','Para-phyllostomid','Rhinolophoidea','Molossidae','Phyllostomidae','Vespertilionidae'),
	values=levels(factor(1:6)),na.translate=F )+
	labs(y='Mya', fill='PCA', color='')+
	scale_x_discrete(limits = rev, expand = expansion(add = c(4.5, 0)) )+
	scale_y_discrete(expand=expansion(add=c(9.5,0)))+
	geom_segment(data=ticks,aes(x=x,xend=xend,y=y,yend=yend),lwd=Legend.tick.size,col='black',lineend='round')+
	geom_text(data=ticks,aes(x=x-1.1,y=y,label=seq(0,60,by=10)), size=Subplot.axis.font.size/.pt, family='Times New Roman')+
	geom_text(aes(x=c(-1,-1),y=c(1,2),label=c('PC1','PC2')), size=Subplot.axis.font.size/.pt, family='Times New Roman')+
	labs(y='Mya')+
	ggnewscale::new_scale_fill()+
	geom_rect(data=manip, aes(xmin=1+length(unique(manip$xmin))-xmin,ymin=ymin,xmax=1+length(unique(manip$xmin))-xmax,ymax=ymax,fill=factor(fill)), fill=manip$fill )+
	theme(legend.position='bottom',
	legend.title=element_text(size=9.5,family='Times New Roman'),
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
		axis.title.y=element_text(size=9.5,family='Times New Roman'),
		axis.ticks.x=element_blank(),
		axis.line.y.left=element_blank(),
		axis.ticks.y.left=element_blank(),
		legend.background = element_rect(fill = 'transparent', color = NA))


 main.legend <-get_plot_component(main.block, "guide-box-bottom")
# I would like to suppress the 'Mya' and 'PC1/PC2' labels and draw them on top with text
# I would like to replace the fonts in all cases with Times New Roman 

main.block <- main.block + theme(legend.position = "none")
# Temporarily remove legend so we can place it in a different location 

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


# Assemble all plot elements 
vert<-0.9
upset<-0.03
ggdraw()+
draw_plot(blank)+
draw_plot(fill.legend, x = 0.41, y = -0.068, width = 0.2, height = 0.4*vert)+
draw_plot(main.legend, x = 0.35, y = 0.02, width = 0.2, height = 0.4*vert)+
draw_plot(main.block,x=0.005,y=0.215,width=0.85,height=0.6*vert)+

# Users replicating our analysis should not run the below lines: 
draw_plot(Aethal+ theme_half_open(12),x=0.0,y=0.6+upset,width=.25,height=0.4*vert)+
draw_plot(Hipposideros+ theme_half_open(12),x=0.127,y=0.6+upset,width=.25,height=0.4*vert)+
draw_plot(Pteronotus+ theme_half_open(12),x=0.23,y=0.62+upset,width=.25,height=0.4*vert)+
draw_plot(Carollia+ theme_half_open(12),x=0.4,y=0.62+upset,width=.25,height=0.4*vert)+
draw_plot(Eumops+ theme_half_open(12),x=0.56,y=0.62+upset,width=.25,height=0.37*vert)+
draw_plot(Myotis+ theme_half_open(12),x=0.66,y=0.62+upset,width=.25,height=0.4*vert)+
draw_plot(Myotis_lateral+ theme_half_open(12),x=0.8,y=0.5+upset,width=.25,height=0.5*vert)+
draw_plot(Tylonycteris_lateral+ theme_half_open(12),x=0.82,y=0.08+upset,width=.25,height=0.5*vert)+

# Users replicating our analysis resume here: 
draw_text(text=c('a','b','c','d','e','f','f','g'),x=c(0.095,0.21,0.33,0.5,0.66,0.76,0.90,0.90),y=c(0.75,0.75,0.75,0.75,0.75,0.75,0.6,0.2)+upset,
size=9.5, family='Times New Roman')

# Save the output; adjust the filepath as necessary 
ggsave(filename='Whatever_name_you_want.png',width=18,height=14,unit='cm',dpi=600,device=png, type='cairo')




