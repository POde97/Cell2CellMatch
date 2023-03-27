import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import scipy
import scanpy as sc
from networkx.algorithms import bipartite as bpr
import pandas as pd 
import networkx as nx 
import numpy as np
from itertools import count
from sknetwork.clustering import Louvain, modularity, bimodularity
import igraph as ig
import leidenalg as la


def Optimazed_Build_G_with_Treshold(matrix1,Treshold):
    G=nx.algorithms.bipartite.matrix.from_biadjacency_matrix(scipy.sparse.csr_matrix(matrix1))
    to_remove = [(a,b) for a, b, attrs in G.edges(data=True) if attrs["weight"] < Treshold]
    G.remove_edges_from(to_remove)
    
    return G


def merge_two_dicts(x, y):
    z = x.copy()   # start with keys and values of x
    z.update(y)    # modifies z with keys and values of y
    return z



#K-partite Graph
def KpG(list_of_sample1,list_of_sample2,sameData,T,OptT):
    list_of_sample_Seurat1 = list_of_sample1
    list_of_sample_Seurat2 = list_of_sample2
    
    T = 2.0
    G = nx.Graph()
    idx1=[i for i in range(len(list_of_sample_Seurat1))]
    idx2 = [i for i in range(len(list_of_sample_Seurat1))]

    for i in range(len(list_of_sample_Seurat1)):#
        if sameData==True:
            list_of_sample_Seurat2.remove(list_of_sample_Seurat1[i])
            idx2.remove(idx1[i])
        for j in range(len(list_of_sample_Seurat2)):#
            if list_of_sample_Seurat1[i] != list_of_sample_Seurat2[j]:
                matrix1 = pd.read_csv("matrixPvalue/"+
                                          list_of_sample_Seurat2[j]+"-"+list_of_sample_Seurat1[i]
                                             +"-dimz-50.csv",index_col=[0])
                #Rename
                matrix1.columns = matrix1.columns.str.replace('.', '-')
                matrix1.index = matrix1.index.str.replace('.', '-')
                matrix1.columns =  matrix1.columns.astype(str) + "-" +list_of_sample_Seurat1[i]
                matrix1.index = matrix1.index +"-" + list_of_sample_Seurat2[j]
                T = 2.0
                G_temp =Optimazed_Build_G_with_Treshold(matrix1,T)
                sample2_name = list(matrix1.index)
                sample1_name= list(matrix1.columns)
                names = sample2_name+ sample1_name
                #print(len(names))
                mapping = dict(zip(G_temp, names))
                G_temp = nx.relabel_nodes(G_temp, mapping)
                
                
                
                if OptT ==True:
                    t = [len(c) for c in sorted(nx.connected_components(G_temp),key = len,reverse = True)]
                    while  len([i for i in t if i>1]) == 1:
                        T = T + 1.0
                        
                        to_remove = [(a,b) for a, b, attrs in G_temp.edges(data=True) if attrs["weight"] < T ]
                        G_temp.remove_edges_from(to_remove)
                        t = [len(c) for c in sorted(nx.connected_components(G_temp),key = len,reverse = True)]
                
                
                T = T-1.0
                if T < 2.0:
                    T = 2.0
                #print(T)
                #G = KpG(list_of_sample_Seurat1,list_of_sample_Seurat2,T=T,sameData=sameDat)
                #print([len(c) for c in sorted(nx.connected_components(G),key = len,reverse = True)])
                G_temp = Optimazed_Build_G_with_Treshold(matrix1,T)
                G_temp = nx.relabel_nodes(G_temp, mapping)
                largest_cc = max(nx.connected_components(G_temp),key = len)
                Go = G_temp.subgraph(largest_cc)
                G_temp = Go
                
                
                
                
                #print(len(G_temp.nodes()))
                #print(idx1[i],idx2[j])
                
                Seurat2_nodes = {n for n, d in G_temp.nodes(data=True) if d["bipartite"] == 0}
                Seurat1_nodes = set(G_temp) - Seurat2_nodes
                
                
                dict_temp1 = dict(zip(list(Seurat2_nodes),
                                          [idx2[j] for kk in range(len(list(Seurat2_nodes)))]))
                dict_temp2 = dict(zip(list(Seurat1_nodes),
                                      [idx1[i] for kk in range(len(list(Seurat1_nodes)))]))
                dict_temp = merge_two_dicts(dict_temp1,dict_temp2)
                
                G.add_nodes_from(G_temp.nodes(data=True))
                nx.set_node_attributes(G, dict_temp , name="bipartite")
                G.add_edges_from(G_temp.edges(data = True))
    return G,T


def BipLouvain(G,bip0,bip1):
    biadjacency =  nx.algorithms.bipartite.matrix.biadjacency_matrix(G, 
                                                                     row_order = list(bip0), 
                                                                     column_order=list(bip1))
                
    louvain = Louvain(resolution= 1.0,random_state = 42)
    louvain.fit(biadjacency,force_bipartite=True)
    labels_row = louvain.labels_row_
    labels_col = louvain.labels_col_
    lbz = list(labels_row)+list(labels_col)
    d = {'node': bip0+bip1,'clusterOverall' : lbz}
    
    return d

def BipLeiden(G):
    h = ig.Graph.from_networkx(G)
    h.vs["type"] = h.vs["bipartite"] 
    optimiser = la.Optimiser()
    p_01, p_0, p_1 = la.CPMVertexPartition.Bipartite(h,
                                                     resolution_parameter_01=(3./(2*len(G.edges()))) ,
                                                     weights='weight',types='bipartite',
                                                     degree_as_node_size = True)
    diff = optimiser.optimise_partition_multiplex([p_01, p_0, p_1],layer_weights=[1, -1, -1],n_iterations=-1)
    d = {'node': h.vs['_nx_name'],'clusterOverall' : p_01.membership} 
    
    return d

def Leiden(G):
    h = ig.Graph.from_networkx(G)
    part = la.find_partition(h, la.ModularityVertexPartition,weights="weight",n_iterations=-1)
                #MemebershipTot
    d = {'node': h.vs['_nx_name'],'clusterOverall' : part.membership}
    
    return d
