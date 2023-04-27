# Read me

#### The evolution of ontogenetic _decision-making_ in the wood of a clade of tropical plants. 

## Description

This repository contains the data and code used for reproducing the results and figures of the paper/manuscript, "The evolution of ontogenetic _decision-making_ in the wood of a clade of tropical plants," by Emilio Petrone-Mendoza, Mariana Benítez, María Elena Lárraga, and Mark E. Olson (UNAM, Mexico).

The repository is divided in: 

1. `Data`: Data that we generated from the wood cell lineages from *Pedilanthus* clade species of the genus *Euphorbia*. 
   - `Raw_data/`
   - `Cell_files_data/` 
   - `euclidean_distance_all/`
   - `euclidean_distance_morethanone/`
   - `words_count_all/`
   - `words_count_morethanone/`
   
2. `Figures`: figures of the manuscript that we generated with scripts and data 
3. `meta`: additional tables and information from the samples.
4. `scripts`: executables for the analysis of cell lineage data. Scripts are written in python and other in R.

##### Description of the `Data` folder.

* `Data`:
  - `Raw_data/`:
    - `892_edited.txt`
    - `896_edited.txt`
    - `939_edited.txt`

The `sample_edited.txt` files are tab delimited files containing three columns: the first one specifying the species of the sample, the second has the cell lineage number, and the third one has the sequence of coded cells.  

| Species       | Cell file number | Cell lineage sequence |
| -------------- | :----:    | :----:   |
| *E. calcarata* |    C1     |  FFFFVPFFFPF... |
| *E. calcarata* |    C2     |  FFFFFFPFPFP... |

In addition to the cell-to-letter code from wood cells, some cell lineage sequences contain three additional metadata with the following syntax and biological meaning: 

* Wood cells delimited by parenthesis and a ^ character (*i.e.* ^P). This code represent cells that, potentially, are not part of the coded lineage but some they were coded because some part of the cell is interrupting the continuous series of cells derived from the coded lineage.  
  
* A Converge- word delimited by parenthesis (Converge-). Two lineages that converge in one cell by anticlinal divisions (2 cell lineages converge in one cell lineage) were marked with the (Converge-). Only one of the two derived cell lineages was coded in the same cell file. 

* Hyphens (-) at the beginning of some cell lineages. Some wood cell lineages were visible not from the last differentiated wood cells near the vascular cambium. When cell lineages were cleary identified nearby other coded cell lineages, we coded cells starting at the relativeposition of the other cells.  

#### Data preprocessing
Using awk commands we removed the first two columns of each file and we redirected the file too a new file within the Cell_files_data folder. We used Sed commands to remove the meta coding because we do not analize these data for the present work.  
Details of the commands use for preprocessing 
``` bash
awk '{print $4'} Data/Raw_data/845_edited.txt > Data/Raw_data/Cell_files_data/845_edited_cells.txt 
sed -i 's/(Converge-)//g' Data/Raw_data/P_macrocarpus/EPM13_edited_cells.txt 
#to add cells as a new cell file after codifying a convergence event add:
sed -Ei 's/\(Converge-\)/&\n/g' 845_edited_cells_NotConverge.txt
```
After editing, we moved the files to the Cell_files_data. Files with convergence events coded as separate cell lineages were saved in `Data/Cell_files_data/ConvergeAssOtherLineage`

#### Workflow

Global functions were created and located in the `wordanalysis.py` python script. The script has functions to count number of cells (letters), length of lineages, word counting, estimating euclidean distance, and measuring Lempel-Ziv compression and Shannon-Entropy estimation methods. The `l-system.py` script generates the virtual wood cell lineages. Previous to starting cell lineage analysis the `l-system.py` is run to include the virtual cell lineages with the rest of the _Pedilanthus_ lineages

##### Cell lineage sequence analysis

To determine the length of cell lineages and cell type frequencies: 

1. count_longitudes.py generates a file named cell_lengths.csv having the lengths of the cell lineages placed in Data/
2. count_transitions.py counts number of cell types and number of words (words from length 2 to 31) and number of words appearing more than ones. Output of the word counts are named as wordcounts*n*.csv and are placed on the Data/word_counts_all folder. The *n* denotes the length of the word. Words counts appearing more than one are placed on the word_counts_morethanone folder.The count_transitionswithR.py makes the same as count_transitions.py but it includes the ray system. 
3. The homogeneity_index.py script generates evaluates the homogeneity indexes for the cell lineages and creates a data frame with the value measure for each cell lineage across samples. 

### Clean Data

* Field experiment root trait data output ('BRTdata.csv')

* Field experiment fitness data ('FitnessData.csv')

* Rhizotron root architecture data ('Rhiz_root_traits.csv')


## Worwflow. 

In the script files, the [wordanalysis.py] script contains several functions used in other scripts such as [count_transitions.py]

#The order of analysis can proceed in the following order:

	1) count_longitudes.py
	2) count_transitions.py
	3) count_transitionswithR.py
	4) homogenity_index.py
	5) wordcomparisonmethods.py####
	6) wordcomplexity.py
	7) wordcomplexity_withR.py
	

### R script

* Calculating root architecture traits ('RootArchitecture_Calc.R')

* Field competition experiment datananalysis ('Field_SelAnalysis_LinearMixedModels')

* Rhizotron greenhouse study analysis ('Rhizotron_2017_LinearMixedModel_Figures.R')

* Useful R script sourced in the main analysis above to compute summary statistics ('SummarySE.R')


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


