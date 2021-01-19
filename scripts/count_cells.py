#Create dictionary of dictionaries
all_cells = {}
# Open the file and evaluate safe every line in a list of P. bracteatus
f=open('../Data/Pedilanthus/P_bracteatus/845_edited_cells.txt')
S845 = count_cells(f)
#Make a data frame with pandas
#S845 = pd.DataFrame.from_dict(counts, orient='index', columns=['Counts'])
#S845.insert(0, "Sample", np.repeat("845", np.shape(S845)[0]), True)
#S845.insert(0, "Species", np.repeat("Euphorbia bracteata", np.shape(S845)[0]), True)
#df.to_csv('../Data/Pedilanthus/P_bracteatus/bracteatus_length_filecells.csv', index = False)

all_cells["S845"] = S845
#############Open files of P. calcaratus
f=open('../Data/Pedilanthus/P_calcaratus/892_edited_cells.txt')
S892 = count_cells(f)

all_cells["S892"] = S892

f=open('../Data/Pedilanthus/P_calcaratus/896_edited_cells.txt')
S896 = count_cells(f)

all_cells["S896"] = S896
#############Open files of P. coalcomanensis
f=open('../Data/Pedilanthus/P_coalcomanensis/883_edited_cells.txt')
S883 = count_cells(f)

all_cells["S883"] = S883

##############Open files of P. colligata connatus
f=open('../Data/Pedilanthus/P_colligata/867_edited_cells.txt')
S867 = count_cells(f)

all_cells["S867"] = S867

##############Open files of P. cymbiferus
f=open('../Data/Pedilanthus/P_cymbiferus/979_edited_cells.txt')
S979 = count_cells(f)

all_cells["S979"] = S979
##############Open files of P. diazluna
f=open('../Data/Pedilanthus/P_diazluna/EPM10_edited_cells.txt')
EPM10 = count_cells(f)

all_cells["EPM10"] = EPM10

f=open('../Data/Pedilanthus/P_diazluna/EPM11_edited_cells.txt')
EPM11 = count_cells(f)

all_cells["EPM11"] = EPM11
f=open('../Data/Pedilanthus/P_diazluna/EPM12_edited_cells.txt')
EPM12 = count_cells(f)

all_cells["EPM12"] = EPM12
#########################Files of E. finkii
f=open('../Data/Pedilanthus/P_finkii/917_edited_cells.txt')
S917 = count_cells(f)
all_cells["S917"] = S917

#########################Files of E. lomelii
f=open('../Data/Pedilanthus/P_macrocarpus/853_edited_cells.txt')
S853 = count_cells(f)

all_cells["S853"] = S853
#########################Files of E. personata
f=open('../Data/Pedilanthus/P_personata/EPM7_edited_cells.txt')
EPM7 = count_cells(f)

all_cells["EPM7"] = EPM7

f=open('../Data/Pedilanthus/P_personata/EPM9_edited_cells.txt')
EPM9 = count_cells(f)

all_cells["EPM9"] = EPM9
#########################Files of E. peritropoides
f=open('../Data/Pedilanthus/P_peritropoides/974_edited_cells.txt')
S974 = count_cells(f)

all_cells["S974"] = S974
#########################Files of E. conzattii
f=open('../Data/Pedilanthus/P_pulchellus/971a_edited_cells.txt')
S917a = count_cells(f)

all_cells["S917a"] = S917a
#########################Files of E. tehuacana
f=open('../Data/Pedilanthus/P_tehuacanus/981_edited_cells.txt')
S981 = count_cells(f)

all_cells["S981"] = S981
#########################Files of E. tithymaloides
f=open('../Data/Pedilanthus/P_tithymaloides/EPM6_S2-1_edited_cells.txt')
EPM6 = count_cells(f)
all_cells["EPM6"] = EPM6
f=open('../Data/Pedilanthus/P_tithymaloides/EPM5_edited_cells.txt')
EPM5 = count_cells(f)
all_cells["EPM5"] = EPM5
#########################Files of E. cyri
f=open('../Data/Pedilanthus/P_tomentellus/973_edited_cells.txt')
S973 = count_cells(f)
all_cells["S973"] = S973

all_counts = pd.DataFrame.from_dict(all_cells, orient='index')
all_counts.to_csv('../Data/Pedilanthus/cell_counts.csv', index = True)


