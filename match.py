import pandas as pd
import numpy as np
from sys import argv

# Matching the spectrum mass list againts the database:
df = pd.read_csv(argv[3])

picklist = open(argv[2], 'r').read()
mass = picklist.split('\n')

# One could define a specific tolerance range (t), default: 1
def matching_spectrum(dataframe, picklist, t):
	candidates = []
	count = 0
	index = []
	for n in dataframe['Peptide Mass']:
		for m in picklist:
			left = n-t
			right = n+t
			if float(m) >= left and float(m)<= right:
				index.append(count)
				break
		count += 1 
	return df.loc[index,]
	#return index

matches = matching_spectrum(df, mass, 1)
# List of matches 
matches.to_csv('matches.csv') 

# Showing the top 10 proteins
top10 = pd.Series.value_counts(matches['Protein'])[:10]
if argv[1] == '--no-scoring':
	print '***** Top10 proteins using PMF *****'
	print top10
	print 'You can use Bowtie scoring to improve your detection by using the argument --use-scoring'

####################################
# Using the MOWSE score matrix
####################################
def Bowtie(match, picklist):
	# Matching the mass of the spectrum with the mass in the database:
	matches = pd.read_csv(match, parse_dates=True)
	mass = picklist.split('\n')

	# Starting the Molecular Weighted Search (MOWSE):
	# Making the matrix where row is size of protein and collumn is the size of fragment

	protlen = df['Protein Mass']
	peplen = df['Peptide Mass']
	together = zip(protlen, peplen)

	# Making the matrix depending on the length of the protein (in KD) and the peptide in D
	cols = int(round(max(peplen)/100))
	rows = int(round(max(protlen)/10000))

	# Placing each data point in the right location
	def place(prot, pep):
		rows = int(max(prot))/10000 + 1
		cols = int(max(pep))/100 + 1
		Mowse_matrix = np.zeros((rows,cols), dtype = 'float')	
		for i in range(len(prot)):
			r = int(prot[i])/10000
			c = int(pep[i])/100
			Mowse_matrix[r,c] += 1
		return Mowse_matrix
	# Filling up the matrix with our data
	Mowse_matrix = place(protlen, peplen)

	# Heatmap plot of the matrix
	#from matplotlib import pyplot as plt
	#heatmap = plt.pcolor(Mowse_matrix[0:50,0:30])
	#plt.ylabel('Protein length')
	#plt.xlabel('Peptide length')
	#plt.show()
	
	# Normalizing the cols of the matrix by their sum
	for j in range(0,cols):
		Mowse_matrix[:,j] = Mowse_matrix[:,j]/float(sum(Mowse_matrix[:,j]))

	# All nan number to zero
	where_are_NaNs = np.isnan(Mowse_matrix)
	Mowse_matrix[where_are_NaNs] = 0

	# Dividing the rows of the matrix by the maximum value
	for i in range(0,rows):
		maxi = np.max(Mowse_matrix[i,:])	
		Mowse_matrix[i,:] = Mowse_matrix[i,:]/float(maxi)

	where_are_NaNs = np.isnan(Mowse_matrix)
	Mowse_matrix[where_are_NaNs] = 0

	# Computing the score for each entry
	prot_test = matches['Protein Mass']
	pep_test = matches['Peptide Mass']

	divisor = []
	for i in range(len(prot_test)):
		r = int(prot_test[i])/10000
		c = int(pep_test[i])/100 
		divisor.append(Mowse_matrix[r,c])

	matches['Divisor'] = pd.Series(divisor)

	grouped = matches.groupby('Protein')
	protein = []
	result = []
	for protname, group in grouped:
		prot_mass = group['Protein Mass']
		prot_mass = np.unique(prot_mass)[0]
		scores = 50000/(np.prod(group['Divisor'])*prot_mass)
		protein.append(protname)
		result.append(scores)

	results = pd.DataFrame({'Protein': protein, 'Scores': result})
	top10 = results.sort_values('Scores', ascending= False )[:10]
	return top10
if argv[1] == '--use-scoring':
	print '***** Top10 proteins using MOWSE score *****'
	print Bowtie('matches.csv', picklist)