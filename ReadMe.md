#####Ontogenetic history of woody cambial cells in the _Pedilanthus_ clade 

This repository contains the cell file/lineage history of _Pedilanthus_ clade species of the genus Euphorbia. Files with the cell file/lineage are in the Data folder. 
Within the scripts folder executables for the analysis of those files

In the script files, the [wordanalysis.py] script contains several functions used in other scripts such as [count_transitions.py]

## Description

This repository contains the data and code used to produce the results figures of the paper,"The evolution of ontogenetic _decision-making_ in the wood of a clade of tropical plants", by Emilio Petrone-Mendoza, Mariana Benítez, María Elena Lárraga and Mark E. Olson (UNAM, Mexico):

#The order of analysis can proceed in the following order:

	1) count_longitudes.py
	2) count_transitions.py
	3) count_transitionswithR.py
	4) homogenity_index.py
	5) wordcomparisonmethods.py####
	6) wordcomplexity.py
	7) wordcomplexity_withR.py
	

# Read me


### R script

* Calculating root architecture traits ('RootArchitecture_Calc.R')

* Field competition experiment datananalysis ('Field_SelAnalysis_LinearMixedModels')

* Rhizotron greenhouse study analysis ('Rhizotron_2017_LinearMixedModel_Figures.R')

* Useful R script sourced in the main analysis above to compute summary statistics ('SummarySE.R')

### Raw Data 

* Coordinate data of individuals measured in the greenhouse rhizotron study ('RhizotronExpCoordinates.csv')

### Clean Data

* Field experiment root trait data output ('BRTdata.csv')

* Field experiment fitness data ('FitnessData.csv')

* Rhizotron root architecture data ('Rhiz_root_traits.csv')

# Key:

# Rhizotron experiment

Id* = plant identification number
Experiment* = temporal block 
Population* = plant population of origin
Code* = code of plant maternal line and population combined
Species* = Species code, Ip = I. purpurea and Ihed = Ihederacea
ML* = Maternal line
PrimaryRootLength* = calculated primary root length
Angle_* = Angle root measurement for the left '_1' and the right '_2' root angles
AvAng* = Average root angle
Total.Area* = Total root system area


# Field experiment

MaxWidth* = Root system width
AreaConvexHull* = Root size
AvAng* = Average root angle
Of_exp* = Cohort; O*= cohort 1, R*= cohort 2
ML* = Maternal line
Combos* = species by maternal line by species maternal line competition pairing
Important note: maternal lines are explicitly nested within population--they only occur within specific populations and species.


