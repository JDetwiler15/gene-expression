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


def PreprocessData(data):
	
	new_data = np.log((data / 1000000) + 1)  # log_2
	# new_data -= new_data.mean(axis=1, keepdims=True)  #  centers the data

	return new_data

def DrawScatterPlot(data, name, labelX='X', labelY='Y'):
	plt.figure()
	plt.scatter(data[0,:], data[1,:])
	plt.title(name)
	plt.xlabel(labelX)
	plt.ylabel(labelY)
	name = "./" + name.replace(" ", "") + str(time.time()) + ".png"
	plt.savefig(name)



""" DEPRICATED 
def TsneSubset(data, d=1000):
	sample = data[0:d,0:d]
	sample = sample[~np.all(sample == 0, axis=1)]  # remove rows
	print("Starting TSNE on matrix with shape = " + str(sample.shape))
	Y = tsne.tsne(sample, 2, 20, 20.0)
	DrawScatterPlot(Y, "TSNE Scatter Plot D=20")

def TsneWrapper(sample):
	smaple = sample[~np.all(sample == 0, axis=1)]
	print("Starting TSNE on matrix with shape = " + str(sample.shape))
	Y = tsne.tsne(sample, 2, 20, 20.0)
	DrawScatterPlot(Y, "TSNE Scatter Plot D=20@" + str(time.time()))
"""	

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
	print(data)
	print("The variable subset shape is :", data.shape)
	return data


def SVDWrapper(data):
	data = sparse.csr_matrix(data.T)
	svd = TruncatedSVD(n_components=50, n_iter=7, random_state=42)
	svd.fit(data)
	return svd.transform(data)

def TSNEWrapper(data):
	return TSNE(n_components=2).fit_transform(data)

def Main():
	
	#data = BuildNumpyArray(True)  # load array
	t0 = time.time()
	data = np.load('matrixPickle.npy')
	data = VariableSubset(data) # make subeset
	data = PreprocessData(data)

	print("Starting SVD at: ", time.time())
	data = SVDWrapper(data) # note this returns data.T as a sparse TrucatedSVD
	print("SVD complete at: " , time.time())
	print("result shape is ", data.shape)
	print("Starting TSNE ... ")
	Y = TSNEWrapper(data)
	print("TSNE finished at ", time.time())
	print(Y)
	print(Y.shape)
	print("Ploting data")
	DrawScatterPlot(Y.T, "TSNETestPlot", labelX='TSNE X', labelY='TSNE Y')
	tf = time.time()-t0
	print("Process took" +str( tf) + "secs" )


	# tsne it
	# pheno cluster

Main()
