import numpy as np
import csv
import pickle
#import matplotlib.pyplot as plt


"""
	input (data) - dxn matrix of features x examples where this is a numpy array
	output - new_data, variances, eigenvectors   note all values are real
"""
def PCA(data):
	
	data -= data.mean(axis=1, keepdims=True)
	
	#u, s, vh = np.linalg.svd(data.T)
	s = np.linalg.svd(data, compute_uv=False)

	return None, s, None

	"""
	eigenvectors = vh.T
	variances = s ** 2
	new_data = eigenvectors.T.dot(data)
	
	return np.real(new_data), np.real(variances), np.real(eigenvectors)
	"""

def buildNumpyArrayQuick():
	data = np.zeros((48276, 15928))
	with open('human_MTG_2018-06-14_exon-matrix.csv') as csv_file:
		csv_reader = csv.reader(csv_file,delimiter=',', quoting=csv.QUOTE_NONE)
		matCount = 0
		for i,row in enumerate(csv_reader):	
			if i > 0:
				geneData = np.array(row[1:], dtype=np.dtype(np.int16))
				x = np.sum(geneData)
				y = np.count_nonzero(geneData)
				if not(x == 0 or y >= geneData.shape[0]-3): 
					data[matCount] = geneData
					matCount += 1
					if matCount % 1000 == 0:
						print(str(matCount) + "-- out of 48272")	
	
		print(matCount)
	print("Saving matrix")

	np.save("matrixPickle", data)
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


def displayPCAOnExonMatrix():
	
	data = np.load('./matrixPickle.npy')
	new_data, variances, eigenvectors = PCA(data)

	plt.figure()
	plt.stem(variances.ravel())
	plt.xlabel('Dimension')
	plt.ylabel('Captured Variance')

	plt.savefig('./VariancePlot.png')

	"""
	plt.figure()
	plt.plot(new_data[0,:], new_data[1:], 'x')
	plt.title('PCA modded data')	
	plt.savefig('./2dPCAPlot.png')
	"""


def main():
	displayPCAOnExonMatrix()


main()
