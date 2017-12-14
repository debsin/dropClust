import numpy as np
from scipy import io
from sklearn.neighbors import NearestNeighbors, LSHForest
from igraph import Graph, EdgeSeq
import random


random.seed(100)

print "Reading sparce matrix..."
matrix = io.mmread("whole_matrix")
print "Converting matrix to dense format..."
a = np.array(matrix.todense())
print a.shape

with open("subsamples_louv_68K") as f:
  idx = f.read().splitlines()
  
with open("pc_gene_ids_sub.csv") as f:
  idx_gene = f.read().splitlines()

idx = np.array(map(int, idx[1:])) - 1
idx_gene = np.array(map(int, idx_gene[1:])) - 1

print "Initialize LSH..."
lshf = LSHForest(random_state=42, n_neighbors=5)
print "fit LSH..."
lshf.fit(a[idx[:,None],idx_gene])

distances, indices = lshf.kneighbors(a[:,idx_gene], n_neighbors=5)

np.savetxt("neigh.txt",indices)
