# bat_sternum_project
[![DOI](https://zenodo.org/badge/1305842103.svg)](https://doi.org/10.5281/zenodo.21443692)

A supporting repository for a manuscript by Orkney and colleagues investigating the evolution of bat sternum shape features.
The manuscript is entitled 'Declining rates of evolution and limited convergence in bat sternum shape can be explained by a ratchet of specialisation'

# Abstract
Strong clustering in animal trait variety between subclades can be explained by evolutionary
bursts early in clade history followed by subsequent decay in rates of change. However, the mechanisms
underlying evolutionary burstiness and the unevenness of diversity across the Tree of Life are still debated. 
Suggested explanations appeal to founder effects and constraint release followed by stabilising
selection. Bat evolution, which evokes founder effects given bats’ innovation of features that resemble
developmental pathologies, is a model system for helping to disentangle whether and how evolutionary
burstiness relates to stabilising selection and adaptive constraint. We examine the evolution of the bat
sternum, a component of the skeleton integral to bat flight. We investigate the relationship between
sternum shape, flight-style and roosting ecology across a representative dataset of bat species and find
that sternum shape differs strongly between taxonomic groups pursuing distinct flight-style ecologies.
Further, variance within adaptive zones is reduced in species that roost in enclosed spaces. We explore
further, estimating the history of roosting ecology across bats, and find that initial bursts in sternum
evolution may decelerate as ecological constraints from roosting accumulate. We hypothesise that roost
ing adaptation constitutes a ratchet of specialisation, explaining bursty rate variation and perhaps bats’
ultimate settlement within adaptive zones structured by flight-style ecology.

# Data
# Data_source_document.csv
  Sources of all original uCT scan data
# sternum.surf.ply 
  A sample mesh topology used to visualise sternum shape and variation
# Bat_CT_process_list_Andrew_only.csv
  A metadata spreadsheet
# Bat_eco_metadata.csv
  Supporting ecological metadata, represented as binary presence/absence data across a wide
sample of bats. Originally aggregated by Dr. Orkney for another project which containsa full file of citations: https://doi.org/10.5281/zenodo.13742209
# '.png'
  Various image files used to decorate plots
# sternum_array_sep_27_2024.RData
  Landmark data describing bat sternum shape, collected by Jeanne Elizabeth Augustin and Beyonca Akers

# Scripts
# Script_Figure_1.R
The purpose of this script is to generate Figure 1 from the manuscript Orkney et al., 2026
'Declining rates of evolution and limited convergence in bat
sternum shape can be explained by a ratchet of specialisation'
The Figure will portray an anatomically labelled 'average bat sternum' that is 
calculated as the mean of all available sterna. 
A spectral decomposition of the feature matrix will be visualized, with different bat groups coloured, 
and the latent axes of variation will be visualised as warps of the average sternum, with changeable
features colourised. 
An estimate of the shared ancestry between species is available from Shi & Rabosky 2015 (https://doi.org/10.1111/evo.12681)
A species-wise binary bag of words is available from Orkney et al., 2025 (https://doi.org/10.1038/s41559-024-02572-9) 
describing bat ecological characteristics. 

# Script_Figure_2.R
The purpose of this script is to produce Figure 1 in the Main manuscript of Orkney et al., 2026. 
'Declining rates of evolution and limited convergence in bat
sternum shape can be explained by a ratchet of specialisation'
The figure (Figure 2) will display a family tree of the bat species studied. 
Distinct groups of bats that practice different combinations of 
flight-style and possess characteristically different sternal morphology will be identified on the tree. 
Representative 3D images of sterna and our geometric measurement scheme used to quantify sternal shape will be presented in
conjunction with this.
An estimate of the shared ancestry between species is available from Shi & Rabosky 2015 https://doi.org/10.1111/evo.12681 
This script was written by Dr. A Orkney and the final version was compiled on July 19th 2026. 
Users aiming to replicate this figure should observe that the 3D mesh files required to render 3D images are not supplied 
for public use because of file size constraints. Users should disable these components of the code. 

# Script_Figure_3.R
The purpose of this script is to generate Figure 3 from the manuscript Orkney et al., 2026
'Declining rates of evolution and limited convergence in bat
sternum shape can be explained by a ratchet of specialisation'
The analysis will determine whether roosting in enclosed spaces represents a constraint on bat thoracic skeleton
evolution. We hypothesise that additional functional constraints upon sternum evolution introduced by enclosed roosting behaviours
have the effect of reducing the available degrees of freedom for neutral evolution, overall restricting the evolutionary variance
of thoracic features.  
The available data are sternal shape landmark constellations across a diversity of bat species.
An estimate of the shared ancestry between species is available from Shi & Rabosky 2015 (https://doi.org/10.1111/evo.12681)
A species-wise binary bag of words is available from Orkney et al., 2025 (https://doi.org/10.1038/s41559-024-02572-9) 
describing roosting ecology, which can be collapsed into a dichotomous
'enclosed/exposed' categorisation. (We use the language 'compressed/free' interchangeably in comments here.) 
This script conducts the core analyses and produces Figure 3 of the main manuscript Orkney et al., 2026. 

# Script_Figure_4.R
The purpose of this script is to generate Figure 4 from the manuscript Orkney et al., 2026
'Declining rates of evolution and limited convergence in bat
sternum shape can be explained by a ratchet of specialisation'
The analyses will determine whether roosting in enclosed spaces represents a constraint on bat thoracic skeleton
evolution and produce Figure 3 in the main manuscript of Orkney et al., 2026. 
The available data are sternal shape landmark constellations across a diversity of bat species.
An estimate of the shared ancestry between species.
A species-wise binary bag of words describing roosting ecology, which can be collapsed into a dichotomous
'exposed/enclosed' categorization.  (We use the language 'compressed/free' interchangeably in comments here.) 
An estimate of the shared ancestry between species is available (Shi & Rabosky 2015 https://doi.org/10.1111/evo.12681) 
A species-wise binary bag of words is available from Orkney et al., 2025 (https://doi.org/10.1038/s41559-024-02572-9)
This script was written by Dr. A Orkney and the final version was compiled on July 19th 2026. 
We will, explicitly, fit a multiple Ornstein-Uhlenbeck model to determine whether gravitating selection 
towards distinct combinations of sternum shape traits and flight-style ecological properties 
is a likely explanation of bat sternum shape evolution through time. 
We hypothesize that the gradual settlement of bats within these adaptive zones, consistent with a Simpsonian
Quantum evolutionary dynamic of adaptive radiation, is influences by second-order adaptive demands imposed
by the relationship between sternum shape variety and exposed/enclosed roosting ecology. 
Specifically, we expect that enclosed roosting ecologies decrease the capacity of sternum shape evolutionary change
and that this changes the gradient descent of bat lineages towards adaptive optima- allowing them to be more easily
captured by the field of stabilizing selection surrounding different flight-style adaptive zones. 
We will further an estimate of ancestral sternum shape states as a series of discrete quantum leaps between
different states attested in living bat species- with the hypothesis that the a reconstruction of this series of 
evolutionary events is likelier if we assume stationary variances of the evolutinary process associated
with exposed roosting ecologies. 

  
