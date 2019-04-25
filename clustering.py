"""
Author - Luke Quinn
Capabilties -   1) Creating numpy pickle file to imporve load speeds
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

# see https://github.com/jacoblevine/PhenoGraph
import phenograph
from sklearn.decomposition import TruncatedSVD
from scipy import sparse
from sklearn.manifold import TSNE




def PCAWrapper(data):
    dataT ,variences = SVDWrapper(data, dim=2) # note this returns data.T as a sparse TrucatedSVD
    Xaxis = "PC-1, Var=" + str(variences[0])
    Yaxis = "PC-2, Var=" + str(variences[1])

    DrawScatterPlot(dataT.T, "PCA of Gene data", labelX=Xaxis, labelY=Yaxis)
    

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




def get_genes_to_delete():
    with open('human_MTG_2018-06-14_genes-rows.csv') as file:
        entrez_ids_to_delete = []
        for l in  csv.reader(file, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
            if 'X' in l[1] or 'Y' in l[1]:
                entrez_id = l[2]
                entrez_ids_to_delete.append(entrez_id)

        return entrez_ids_to_delete

def get_samples_to_delete():
    with open('human_MTG_2018-06-14_samples-columns.csv') as file:
        samples_to_delete = []
        next(file)
        sample_name = 0
        total_reads = 19
        percent_exon = 20
        percent_intron = 21
        percent_aligned_reads_total = 28
        reads_unique = 25
        for l in  csv.reader(file, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
            # the total value doesn't change anything so may not be necessary
            total_aligned_to_exon_intron = (float(l[percent_exon]) + float(l[percent_intron])) * float(l[total_reads]) / 100.0
            if float(l[reads_unique]) <= 50 or float(l[percent_aligned_reads_total]) <= 40 or total_aligned_to_exon_intron <= 500000:
                samples_to_delete.append(l[sample_name])

        return samples_to_delete

def BuildFilteredNumpyArray(save=False):
    genesToDelete = get_genes_to_delete()
    samplesToDelete = get_samples_to_delete()

    data = np.zeros((50281 - len(genesToDelete), 15928),dtype=np.float64)
    with open('human_MTG_2018-06-14_exon-matrix.csv') as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',', quoting=csv.QUOTE_NONE)
        matCount = 0
        samples_header = None
        for i,row in enumerate(csv_reader):
            if i == 0:
                samples_header = row
            if i > 0 and row[0].strip('\"') not in genesToDelete:
                geneData = np.array(row[1:], dtype=np.dtype(np.float64))
                data[matCount] = geneData
                matCount += 1
        columnIndicesToDelete = []

        for i,name in enumerate(samples_header):
            if name.strip('\"') in samplesToDelete:
                columnIndicesToDelete.append(i)
        data = np.delete(data, columnIndicesToDelete, axis=1)
        #remove 0 rows
        data = data[~np.all(data == 0, axis=1)]


    if save:
        print("Saving matrix")
        np.save("matrixPickle", data)

    print("Data read but not reduced here")
    return data, columnIndicesToDelete


def DisplayPCAOnExonMatrix(data):
    
    #if not data:
    #   data = np.load('./matrixPickle.npy')

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
    plt.scatter(data[0,:], data[1,:], c=colour_seq, s=1)
    plt.title(name)
    plt.xlabel(labelX)
    plt.ylabel(labelY)
    name = "./" + name.replace(" ", "") + str(time.time()) + ".png"
    plt.savefig(name, transparent=True)

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
    return svd.transform(data), svd.explained_variance_ratio_

def TSNEWrapper(data):
    return TSNE(n_components=2).fit_transform(data)

def TSNEPipeline(data):

    print("Starting SVD at: ", time.time())
    data = SVDWrapper(data, dim=20)[0] # note this returns data.T as a sparse TrucatedSVD
    print("SVD complete at: " , time.time())
    print("result shape is ", data.shape)
    print("Starting TSNE ... ")
    Y = TSNEWrapper(data)
    print("TSNE finished at ", time.time())
    print(Y.shape)
    return Y.T


def SubsetsBasedOnType(data, ignoreIdx):

    keep_cols_exc = []
    keep_cols_inb = []
    keep_cols_non_nerual = []
        
    with open('human_MTG_2018-06-14_samples-columns.csv') as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count in ignoreIdx:
                line_count += 1
                continue
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

def AgglogClusterSubsets(data, subsets_indices):

    print("Clustering Exc ")
    cluster1 = agg(n_clusters=24).fit(data[subsets_indices[0]])
    print("Clustering Inb ")
    cluster2 = agg(n_clusters=45).fit(data[subsets_indices[1]])
    print("Clustering Non-Neural ")
    cluster3 = agg(n_clusters=6).fit(data[subsets_indices[2]])

    exc_colours = cluster1.labels_
    inb_colours = cluster2.labels_
    non_neural_colours = cluster3.labels_

    print(np.max(exc_colours))
    print(np.max(inb_colours))
    print(np.max(non_neural_colours))

    # builds a colour sequences that increases in number
    inb_colours += (np.max(exc_colours) + 1)
    non_neural_colours += (np.max(inb_colours) + 1)
    print(np.max(non_neural_colours))
    return (exc_colours, inb_colours, non_neural_colours)

def ColourPlot(data, name, subsets_indices, colours, labelX1="", labelY1=""):
    
    print("Building Colour Sequence")   
    # 3 way merge sort
    idxT, colourT = merge_arrays_inorder((subsets_indices[0], subsets_indices[1]), (colours[0], colours[1]))
    inorder, colour_seq = merge_arrays_inorder((idxT, subsets_indices[2]), (colourT, colours[2]))
    colour_seq = colour_seq[:-1] # chop the last one off

    print("Ploting Data")   
    DrawScatterPlot(data, name, labelX=labelX1, labelY=labelY1, colour_seq=colour_seq)
    
def ColourPlotAdv(data, name, subsets_indices, colours):

    idxT, colourT = merge_arrays_inorder((subsets_indices[0], subsets_indices[1]), (colours[0], colours[1]))
    inorder, colour_seq = merge_arrays_inorder((idxT, subsets_indices[2]), (colourT, colours[2]))
    colour_seq = colour_seq[:-1] # chop the last one off

    fig, ax = plt.subplots(1,1, figsize=(6,6))

    cmap = plt.cm.gist_ncar
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = (.5, .5, .5, 1.0)

    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace(0,1,75)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    scat = ax.scatter(data[0,:], data[1,:], c=colour_seq,s=1, cmap=cmap, norm=norm)
    ax.set_title("Better colors?")
    name = "./" + name.replace(" ", "") + str(time.time()) + ".png"
    plt.savefig(name, transparent=True)

def Main():

    #data = BuildNumpyArray(True)  # load array
    data, removedIdx = BuildFilteredNumpyArray()
    #data = np.load('matrixPickle.npy') # load saved data   

    data = VariableSubset(data) # make subeset

    print("Spliting Based on type")
    subsets_indices = SubsetsBasedOnType(data, removedIdx)
    
    #print("PhenoClustering")
    #colours = AgglogClusterSubsets(data, subsets_indices)
    colours = PhenoClusterSubsets(data, subsets_indices)

    #PCAWrapper(data)

    print("Start TSNE")
    projection = TSNEPipeline(data)
    #DrawScatterPlot(projection,"Gene data mapping using T-SNE Projection", labelX="T-SNE X", labelY="T-SNE Y")
    ColourPlot(projection,"Gene data TSNE projection coupled with Phenograph clustering", subsets_indices, colours, labelX1="TSNE Dim 1", labelY1="TSNE Dim 2")


if __name__ == "__main__":
    Main()
