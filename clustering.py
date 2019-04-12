"""
Author - Luke Quinn
Capabilties -	1) Creating numpy pickle file to imporve load speeds
				2) General purpose PCA
	
Note this requires you to have the exon.csv file present
"""
import numpy as np
import csv
import pickle
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time

import phenograph
from sklearn.decomposition import TruncatedSVD
from scipy import sparse
from sklearn.manifold import TSNE

"""
	input (data) - dxn matrix of features x examples where this is a numpy array
	output - new_data, variances, eigenvectors   note all values are real
"""
def PCA(data, FullPCA=False):
	
	if not FullPCA:
		data -= data.mean(axis=1, keepdims=True)
		
		#u, s, vh = np.linalg.svd(data.T)
		s = np.linalg.svd(data, compute_uv=False)	

		return None, np.real(s**2), None
	else:	
		u, s, vh = np.linalg.svd(data.T)
		eigenvectors = vh.T
		variances = s ** 2
		new_data = eigenvectors.T.dot(data)
		
		return np.real(new_data), np.real(variances), np.real(eigenvectors)
	

def BuildNumpyArray(save=False):
	
	data = np.zeros((50281, 15928),dtype=np.float64)
	with open('human_MTG_2018-06-14_exon-matrix.csv') as csv_file:
		csv_reader = csv.reader(csv_file,delimiter=',', quoting=csv.QUOTE_NONE)
		matCount = 0
		for i,row in enumerate(csv_reader):	
			if i > 0:
				geneData = np.array(row[1:], dtype=np.dtype(np.float64))
				data[matCount] = geneData
				matCount += 1
	if save:
		print("Saving matrix")
		np.save("matrixPickle", data)

	print("Data read but not reduced here")
	return data


def BuildNumpyArrayReducedDim(save=False):
	data = np.zeros((48276, 15928),dtype=np.float64)
	with open('human_MTG_2018-06-14_exon-matrix.csv') as csv_file:
		csv_reader = csv.reader(csv_file,delimiter=',', quoting=csv.QUOTE_NONE)
		matCount = 0
		for i,row in enumerate(csv_reader):	
			if i > 0:
				geneData = np.array(row[1:], dtype=np.dtype(np.float64))
				x = np.sum(geneData)
				y = np.count_nonzero(geneData)
				if not(x == 0 or y >= geneData.shape[0]-3): 
					data[matCount] = geneData
					matCount += 1
		
	if save:
		print("Saving matrix")
		np.save("matrixPickle", data)

	print("Data read")
	return data

def QualityGeneExpressionQualityControlRows():
	badRows = []

	with open('human_MTG_2018-06-14_exon-matrix.csv') as csv_file:
		csv_reader = csv.reader(csv_file,delimiter=',', quoting=csv.QUOTE_NONE)
		for i,row in enumerate(csv_reader):
			if i > 0:
				geneData = np.array(row[1:], dtype=np.dtype(np.int16))
				x = np.sum(geneData)
				y = np.count_nonzero(geneData)
				if x == 0 or y >= geneData.shape[0]-3: 
					badRows.append(i)

	print(str(len(badRows)) + " this many bad rows")
	return badRows 

def GenerateQualityControlIdLists():

	sampleId = []
	sampleColumnIndex = []

	with open('human_MTG_2018-06-14_samples-columns.csv') as csv_file:
		csv_reader = csv.reader(csv_file,delimiter=',')
		line_count = 0
		for row in csv_reader:
			if line_count > 0:
				if row[19] >= 500000 and row[25] > 50 and row[28] > 40:
					sampleId.append(row[1])
					sampleColumnIndex.append(line_count)
			line_count += 1

	return sampleId,sampleColumnIndex


def DisplayPCAOnExonMatrix(data):
	
	#if not data:
	#	data = np.load('./matrixPickle.npy')

	new_data, variances, eigenvectors = PCA(data)

	y = []
	total = 0
	for point in variances:
		y.append(point)
		total += point
		if total > 50:
			print("It took " + str(len(y)) + " to account for " + str(total) + " of the variance")
			break
		
	x = [i for i in range(len(y))]

	plt.figure()
	plt.plot(x,y)
	plt.xlabel('Dimension')
	plt.ylabel('Captured Variance')

	plt.savefig('./VariancePlot.png')

	np.save('Variances', variances)

	"""
	plt.figure()
	plt.plot(new_data[0,:], new_data[1:], 'x')
	plt.title('PCA modded data')	
	plt.savefig('./2dPCAPlot.png')
	"""

def Testsetup():
	y = np.random.random_integers(1,100, 100)
	plt.figure()
	plt.stem(y.ravel())
	plt.xlabel('Linear')
	plt.ylabel('Random')
	plt.savefig('./Temp.png')

def DrawScatterPlot(data, name, labelX='X', labelY='Y', colour_seq=None):
	plt.figure()
	plt.scatter(data[0,:], data[1,:], c=colour_seq)
	plt.title(name)
	plt.xlabel(labelX)
	plt.ylabel(labelY)
	name = "./" + name.replace(" ", "") + str(time.time()) + ".png"
	plt.savefig(name)

def VariableSubset(data):
	# cut cells who have less than 1000 genes
	
	keep_cells = np.count_nonzero(data, axis=0) > 1000
	data = data[:,keep_cells]
	#cut genes with low variabilty
	keep_genes = np.std(data, axis=1) > 0
	data = data[keep_genes]

	data /= np.sum(data, axis=0)  # convert to percents
	data *= 1e6  # convert percents to cpm
	data += 1   # add 1 for log transform
	data = np.log(data)  # log transform
	print("The variable subset shape is :", data.shape)
	return data


def SVDWrapper(data, dim=50):
	data = sparse.csr_matrix(data.T)
	svd = TruncatedSVD(n_components=dim, n_iter=7, random_state=42)
	svd.fit(data)
	return svd.transform(data)

def TSNEWrapper(data):
	return TSNE(n_components=2).fit_transform(data)

def TSNEPipeline(data):

	print("Starting SVD at: ", time.time())
	data = SVDWrapper(data, dim=20) # note this returns data.T as a sparse TrucatedSVD
	print("SVD complete at: " , time.time())
	print("result shape is ", data.shape)
	print("Starting TSNE ... ")
	Y = TSNEWrapper(data)
	print("TSNE finished at ", time.time())
	print(Y.shape)
	return Y.T

"""
def SubsetBasedOnClass(data):

	keep_cols_gluta = []
	keep_cols_ = []
	with open('human_MTG_2018-06-14_samples-columns.csv') as csv_file:
		csv_reader = csv.reader(csv_file,delimiter=',')
		line_count = 0
		for row in csv_reader:
			if line_count > 0:
				clusterClass = 
					sampleId.append(row[1])
					sampleColumnIndex.append(line_count)
			line_count += 1

	return sampleId,sampleColumnIndex
"""

def SubsetsBasedOnType(data):

	keep_cols_exc = []
	keep_cols_inb = []
	keep_cols_non_nerual = []
		
	with open('human_MTG_2018-06-14_samples-columns.csv') as csv_file:
		csv_reader = csv.reader(csv_file,delimiter=',')
		line_count = 0
		for row in csv_reader:
			if "Inh" in row[-1]:
				keep_cols_inb.append(line_count)
			elif "Exc" in row[-1]:
				keep_cols_exc.append(line_count)
			else:
				keep_cols_non_nerual.append(line_count)
			line_count += 1
	
	return (keep_cols_exc, keep_cols_inb, keep_cols_non_nerual)

def PhenoClusterSubsets(data, subsets_indices):

	print("Clustering 1")
	exc_colours, graph, Q = phenograph.cluster(data[subsets_indices[0]])
	print("Clustering 2")
	inb_colours, graph, Q = phenograph.cluster(data[subsets_indices[1]])
	print("Clustering 3")
	non_neural_colours, graph, Q = phenograph.cluster(data[subsets_indices[2]])

	print(np.max(exc_colours))
	print(np.max(inb_colours))
	print(np.max(non_neural_colours))

	exc_colours += 1
	
	inb_colours += (np.max(exc_colours) + 1)
	non_neural_colours += (np.max(inb_colours) + 1) 
	return (exc_colours, inb_colours, non_neural_colours)

def merge_arrays_inorder(idx, poll_array):

	idx1 = 0
	idx2 = 0

	orderedIdx = []
	orderedPoll = []

	while idx1 + idx2 < len(idx[0]) + len(idx[1]):
		if idx1 >= len(idx[0]):
			orderedIdx.append(idx2)
			orderedPoll.append(poll_array[1][idx2])
			idx2 += 1
		elif idx2 >= len(idx[1]):
			orderedIdx.append(idx1)
			orderedPoll.append(poll_array[0][idx1])
			idx1 += 1
		elif idx[0][idx1] < idx[1][idx2]:
			orderedIdx.append(idx1)
			orderedPoll.append(poll_array[0][idx1])
			idx1 += 1
		else:
			orderedIdx.append(idx2)
			orderedPoll.append(poll_array[1][idx2])
			idx2 += 1

	return orderedIdx, orderedPoll

def ColourPlot(data, name, subsets_indices, colours):
	
	print("Building Colour Sequence")	
	# 3 way merge sort
	idxT, colourT = merge_arrays_inorder((subsets_indices[0], subsets_indices[1]), (colours[0], colours[1]))
	inorder, colour_seq = merge_arrays_inorder((idxT, subsets_indices[2]), (colourT, colours[2]))
	colour_seq = colour_seq[:-1] # chop the last one off
		


	print("Ploting Data")	
	DrawScatterPlot(data, name, colour_seq=colour_seq)
	

def Main():
	
	#data = BuildNumpyArray(True)  # load array
	data = np.load('matrixPickle.npy') # load saved data
	data = VariableSubset(data) # make subeset
	
	print("Spliting Based on type")	
	subsets_indices = SubsetsBasedOnType(data)
	
	print("PhenoClustering")
	colours = PhenoClusterSubsets(data, subsets_indices)

	print("Start TSNE")
	projection = TSNEPipeline(data)
	
	ColourPlot(projection,"Colour Me crazy", subsets_indices, colours)

Main()
