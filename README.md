---
output:
  pdf_document: default
  html_document: default
---
# Read me

#### The evolution of ontogenetic _decision-making_ in the wood of a clade of tropical plants. 

## Description of the data and file structure

This repository contains the data and code used for reproducing the results and figures of the manuscript titled, "The evolution of ontogenetic _decision-making_ in the wood of a clade of tropical plants," by Emilio Petrone-Mendoza, Mariana Benítez, María Elena Lárraga, and Mark E. Olson (UNAM, Mexico).

The repository is divided in: 

1. `Data`: Data that we generated from the wood cell lineages from *Pedilanthus* clade species of the genus *Euphorbia*. 
   - `Raw_data/`
   - `Cell_files_data/` 
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

Each file within the `Raw_data` folder comes from a different sample. The alphanumeric code before the `_edited.txt` is from the sample. The `sample_edited.txt` files are tab delimited files containing three columns: the first one specifying the species of the sample, the second has the cell lineage number, and the third one has the sequence of coded cells.  

| Species       | Cell file number | Cell lineage sequence |
| -------------- | :----:    | :----:   |
| *E. calcarata* |    C1     |  FFFFVPFFFPF... |
| *E. calcarata* |    C2     |  FFFFFFPFPFP... |

In addition to the cell-to-letter code from wood cells, some cell lineage sequences contain three additional metadata with the following syntax and biological meaning: 

* Wood cells delimited by parenthesis and a ^ character (*i.e.* ^P). This code represent cells that, can be intrusive gorwth fibers or other cells interrupting the continuous series of cells derived from the coded lineage.  
  
* A Converge- word delimited by parenthesis (Converge-). Two lineages that converge in one cell by anticlinal divisions (2 cell lineages converge in one cell lineage) were marked with the (Converge-). Only one of the two derived cell lineages was coded in the same cell file, and then processed to be another cell lineage.

* Hyphens (-) at the beginning of some cell lineages. Some wood cell lineages were visible not from the last differentiated wood cells at the vascular cambium. Series of cell lineages identified nearby other coded cell lineages, we coded cells starting at the relative position of the other cells.  

### Data preprocessing
Using awk commands we removed the first two columns of each file and we redirected the file too a new file within the Cell_files_data folder. We used Sed commands to remove the meta coding because we did not analyze these data for the present work.  
Here are some details of the commands use for processing  the text files:
``` bash
awk '{print $4'} Data/Raw_data/845_edited.txt > Data/Raw_data/Cell_files_data/845_edited_cells.txt 
sed -i 's/(Converge-)//g' Data/Raw_data/P_macrocarpus/EPM13_edited_cells.txt 
#to add cells as a new cell file after codifying a convergence event add:
sed -Ei 's/\(Converge-\)/&\n/g' 845_edited_cells_NotConverge.txt
```
We moved the edited files to the `Data/Cell_files_data`. Files with convergence events coded as separate cell lineages were saved in `Data/Cell_files_data/ConvergeAssOtherLineage`

### Workflow
Global functions were created and located in the `wordanalysis.py` python script. The script has functions to count number and type of cells, length of lineages, word counting, estimating euclidean distance, and measuring Lempel-Ziv compression and Shannon-Entropy estimation methods. The `l-system.py` script generates the virtual wood cell lineages. Previous to starting the cell lineage analysis we run the `l-system.py` to include the virtual cell lineages with the rest of the _Pedilanthus_ lineages.
The workflow can be divided in two parts, the first part includes python scripts returning data frames with values from the cell lineage analysis, such as word counts, homogeneity indexes values or cell lineage lengths. The second part includes r scripts returning figures and performing statistical tests. 

#### Cell lineage analysis

The order in which we analyzed data is the following: 

##### First part. Python scripts:

1. To determine the length of cell lineages and cell type frequencies we run the `count_longitudes.py` which generated a file named `cell_lengths_notConverge.csv` containing the lengths of each cell lineage from each individual. The `cell_lengths_notConverge.csv` is placed in Data/. Also, the script returns a file named `cell_lengths_withoutR.csv` with the cell length of the lineages without including ray cells.

2. To determine the total number of words at different _k-mer_ lengths we run the `count_transitions.py` script, which counts number of cell types and number of words and number of words appearing more than once (words from _k-mer_ length 2 to 34). Output of the word counts at each _k-mer_ length are named as wordcounts*n*.csv, where _n_ is the _k-mer_. Files are located on the Data/word_counts_all and in the Data/word_counts_morethanone. Additionaly the script generates two files resuming the total number of words observed at each _k-mer_ length, from 2 to 34, for each individual. One file for all thw eords and the other for the words appearing more than once. The files are named as wordcounts_all.csv.  

4. To determine the Lempel-Ziv compression algorithm values and the Shannon entropy values of the cv lineages from our individuals we run the `wordcomplexity.py` script which generates one file named `shannonentropy.csv` and another named `lemplzivbyfile.csv`.

5. To determine the homogeneity indexes values for each cell lineage we run the `homogeneity_index.py` script which generates returns a data frame with the values for all cell lineages and individuals. The name of the file is `homogentiy_index.csv` and is located in the Data folder.

##### Second part. Rscripts

The second part performs statistical tests and generates figures. In the following list we describe the content and what does each different script: 

- `information_theoryMetrics.R`: we make plots from the lempel-ziv compression algorithm and the Shannon entropy values. 
- `homogenity_index.r`: we make plots for the homogenity indexes. 
- `distance_plots.R`: we make plots related to word counts, cell lengths, and dissimilitud analysis. We run the Bray curtis dissimilitud metric using vegan package, and the NDMS.
- `geographic_info.R`: in this script we extract climate information based on the geographic coordinates of the sampled individuals.



